#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Bioproject-first FASTQ sorter with concurrency and project-level QC/notes.

Metadata (TSV):
  - If --header is set, the first row is skipped.
  - Fixed 1-based columns:
      col2  = Run/Accession
      col19 = BioProject
      col28 = LibraryLayout  (ignored for PE/SE decision)
      col30 = Platform

Behavior:
  1) For each accession, move FASTQs from --fq-dir into <OUT>/<BioProject>/00_fq/
     Expected names ONLY:
       PAIRED: {acc}_1.fastq(.gz), {acc}_2.fastq(.gz)
       SINGLE: {acc}.fastq(.gz)
     If an accession from metadata has no files in --fq-dir -> write it to
       <OUT>/<BioProject>/missing_sra.list
  2) Per BioProject:
       - Determine PE vs SE by actual files under 00_fq/ (metadata ignored).
       - If an accession has BOTH pair (_1/_2) AND single (*.fastq(.gz)): mark as anomaly
         * append accession to <BioProject>/sra_3_fq.note
         * move the "third" single file to <BioProject>/sra_3_fq/
       - Create one read-type marker file:
         * pe.reads   (all non-anomalous are PE)
         * se.reads   (all non-anomalous are SE)
         * pe_se.reads (mixed among non-anomalous)
  3) Platform notes:
       - Collect platforms (from metadata col30) per BioProject
       - Map to buckets: illumina | bgi | roche454 | iontorrent | unknown
       - Create ONE note file named like: "<joined_by_underscore>.platform.note"
         e.g., "illumina_bgi.platform.note"

CLI:
  -m/--metadata  TSV path
  -f/--fq-dir    FASTQ dir (source)
  -o/--out-root  Output root (default: .)
  -t/--threads   Parallel workers for moving files
  --header       Skip first row in metadata
  --dry-run      Do not actually move/touch files; just print actions
"""

import argparse
import csv
import os
import shutil
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

from tqdm import tqdm


# ---------- Platform mapping ----------

def classify_platform_bucket(s: str) -> str:
    s = (s or "").lower()
    if any(k in s for k in ("illumina", "hiseq", "miseq", "novaseq", "nextseq")):
        return "illumina"
    if any(k in s for k in ("bgi", "bgiseq", "mgiseq", "mgitech", "dnbseq")):
        return "bgi"
    if ("roche" in s) or ("454" in s):
        return "roche454"
    if ("ion" in s) or ("torrent" in s):
        return "iontorrent"
    return "unknown"


# ---------- Files expected per accession ----------

def expected_paths(fq_dir: Path, acc: str) -> List[Path]:
    """Return the 6 fixed-form candidates (if exist)."""
    cands = [
        fq_dir / f"{acc}_1.fastq",
        fq_dir / f"{acc}_1.fastq.gz",
        fq_dir / f"{acc}_2.fastq",
        fq_dir / f"{acc}_2.fastq.gz",
        fq_dir / f"{acc}.fastq",
        fq_dir / f"{acc}.fastq.gz",
    ]
    return [p for p in cands if p.exists() and p.is_file()]


def move_file(src: Path, dst_dir: Path, dry_run: bool) -> bool:
    dst_dir.mkdir(parents=True, exist_ok=True)
    dst = dst_dir / src.name
    if dry_run:
        print(f"DRY-RUN: mv -n '{src}' '{dst_dir}/'")
        return True
    if dst.exists():
        return False
    shutil.move(str(src), str(dst))
    return True


# ---------- Read metadata ----------

def iter_rows(meta: Path, skip_header: bool) -> Iterable[Tuple[str, str, str]]:
    """
    Yield (run, bioproject, platform_raw). Skip rows missing run or bioproject.
    col2=run (idx1), col19=bioproject (idx18), col30=platform (idx29)
    """
    with meta.open("r", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        if skip_header:
            next(reader, None)
        for row in reader:
            if len(row) < 30:
                continue
            run = row[1].strip() if len(row) > 1 else ""
            bp  = row[18].strip() if len(row) > 18 else ""
            pf  = row[29].strip() if len(row) > 29 else ""  # may be empty -> unknown
            if not run or not bp:
                continue
            yield run, bp, pf


# ---------- Per-project PE/SE/anomaly analysis ----------

def analyze_project(proj_dir: Path, dry_run: bool) -> Tuple[int, int, int]:
    """
    Inspect <proj_dir>/00_fq, decide PE/SE per accession, handle anomalies.
    Returns (n_pe, n_se, n_anom) counting non-anomalous PE/SE and anomaly count.
    Side effects:
      - moves single file of 3-file anomaly to <proj_dir>/sra_3_fq/
      - writes <proj_dir>/sra_3_fq.note (append)
      - writes one of {pe.reads,se.reads,pe_se.reads}
    """
    fqdir = proj_dir / "00_fq"
    if not fqdir.exists():
        return (0, 0, 0)

    # collect files by accession
    seen: Dict[str, Dict[str, Path]] = defaultdict(dict)  # acc -> {"r1":Path, "r2":Path, "se":Path}
    for p in fqdir.glob("*.fastq*"):
        name = p.name
        if name.endswith(".fastq.gz"):
            stem = name[:-9]  # strip .fastq.gz
        elif name.endswith(".fastq"):
            stem = name[:-6]
        else:
            continue
        if stem.endswith("_1"):
            acc = stem[:-2]
            seen[acc]["r1"] = p
        elif stem.endswith("_2"):
            acc = stem[:-2]
            seen[acc]["r2"] = p
        else:
            acc = stem
            seen[acc]["se"] = p

    n_pe = n_se = n_anom = 0
    anom_dir = proj_dir / "sra_3_fq"
    anom_note = proj_dir / "sra_3_fq.note"
    # clean previous notes/markers to avoid stale state
    for marker in ("pe.reads", "se.reads", "pe_se.reads"):
        m = proj_dir / marker
        if m.exists() and not dry_run:
            m.unlink()

    for acc, parts in seen.items():
        has_pair = ("r1" in parts) and ("r2" in parts)
        has_se   = ("se" in parts)
        if has_pair and not has_se:
            n_pe += 1
        elif has_se and not has_pair:
            n_se += 1
        elif has_pair and has_se:
            # anomaly: move the single file to sra_3_fq/, note accession
            n_anom += 1
            if not dry_run:
                anom_dir.mkdir(exist_ok=True, parents=True)
                with anom_note.open("a") as f:
                    f.write(acc + "\n")
                # move single (the "third") file
                try:
                    se_path = parts["se"]
                    dst = anom_dir / se_path.name
                    if not dst.exists():
                        shutil.move(str(se_path), str(dst))
                except Exception as e:
                    print(f"Warning: failed moving third file for {acc}: {e}", file=sys.stderr)
            else:
                print(f"DRY-RUN: anomaly {acc} -> would write {anom_note.name} and move single to sra_3_fq/")

    # write project-level read-type marker
    marker_name = "pe_se.reads"
    if n_anom == 0:
        if n_pe > 0 and n_se == 0:
            marker_name = "pe.reads"
        elif n_se > 0 and n_pe == 0:
            marker_name = "se.reads"
        elif n_pe == 0 and n_se == 0:
            # no files; keep default mixed to avoid misleading
            marker_name = "pe_se.reads"
    if not dry_run:
        (proj_dir / marker_name).touch()

    return (n_pe, n_se, n_anom)


# ---------- Main ----------

def main():
    ap = argparse.ArgumentParser(
        description="Bioproject-first FASTQ sorter with concurrency, PE/SE audit, anomaly handling, and platform notes."
    )
    ap.add_argument("-m", "--metadata", required=True, help="Metadata TSV path")
    ap.add_argument("-f", "--fq-dir", default="fq", help="Source FASTQ directory")
    ap.add_argument("-o", "--out-root", default=".", help="Output root (BioProject folders here)")
    ap.add_argument("-t", "--threads", type=int, default=min(16, os.cpu_count() or 8), help="Parallel workers")
    ap.add_argument("--header", action="store_true", help="Skip the first row in metadata")
    ap.add_argument("--dry-run", action="store_true", help="Print actions without changing files")
    args = ap.parse_args()

    meta = Path(args.metadata)
    src = Path(args.fq_dir)
    out = Path(args.out_root)

    if not meta.exists():
        sys.exit(f"Error: metadata not found: {meta}")
    if not src.exists():
        sys.exit(f"Error: FASTQ dir not found: {src}")

    out.mkdir(parents=True, exist_ok=True)

    # Build accession -> (bioproject, platform_raw) map and per-project platform buckets
    acc2bp: Dict[str, str] = {}
    bp2plats: Dict[str, set] = defaultdict(set)

    rows = list(iter_rows(meta, skip_header=args.header))
    for run, bp, pf in rows:
        acc2bp[run] = bp
        bp2plats[bp].add(classify_platform_bucket(pf))

    # Stage 1: plan moves / missing lists
    move_tasks: List[Tuple[Path, Path]] = []   # (src_path, dst_dir)
    missing_by_bp: Dict[str, List[str]] = defaultdict(list)

    with tqdm(total=len(acc2bp), desc="Indexing", unit="acc") as pbar:
        for acc, bp in acc2bp.items():
            files = expected_paths(src, acc)
            proj_fq_dir = out / bp / "00_fq"
            if not files:
                # remember missing accession for this project
                missing_by_bp[bp].append(acc)
            else:
                for p in files:
                    move_tasks.append((p, proj_fq_dir))
            pbar.update(1)

    # Stage 2: move in parallel
    if move_tasks:
        with ThreadPoolExecutor(max_workers=max(1, args.threads)) as ex:
            futs = [ex.submit(move_file, s, d, args.dry_run) for (s, d) in move_tasks]
            for _ in tqdm(as_completed(futs), total=len(futs), desc="Moving", unit="file"):
                pass

    # Write missing_sra.list per project
    for bp, miss in missing_by_bp.items():
        proj_dir = out / bp
        if miss:
            path = proj_dir / "missing_sra.list"
            if args.dry_run:
                print(f"DRY-RUN: would write {path} with {len(miss)} accessions")
            else:
                proj_dir.mkdir(parents=True, exist_ok=True)
                with path.open("w") as f:
                    f.write("\n".join(sorted(set(miss))) + "\n")

    # Stage 3: per-project analysis (PE/SE/anomaly + platform note)
    for bp in tqdm(sorted(set(acc2bp.values())), desc="Per-project QC", unit="bp"):
        proj_dir = out / bp

        # PE/SE/anomaly
        n_pe, n_se, n_anom = analyze_project(proj_dir, args.dry_run)

        # Platform note (from metadata collected)
        plats = bp2plats.get(bp, set()) or {"unknown"}
        name = "_".join(sorted(plats)) + ".platform.note"
        note = proj_dir / name
        if args.dry_run:
            print(f"DRY-RUN: would touch {note} (platforms: {sorted(plats)})")
        else:
            proj_dir.mkdir(exist_ok=True, parents=True)
            note.touch()

    print("Done.")


if __name__ == "__main__":
    main()
