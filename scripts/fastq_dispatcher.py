#!/usr/bin/env python3

# date: 2024-08-30 latest update
# contact: yanpengch@qq.com
#  
"""
Bioproject-first FASTQ sorter with concurrency, PE/SE + platform split.

Metadata (TSV):
  - If --header is set, the first row is skipped.
  - Fixed 1-based columns:
      col2  = Run/Accession
      col19 = BioProject
      col28 = LibraryLayout  (ignored for PE/SE decision)
      col30 = Platform

Behavior:
  1) For each accession, move FASTQs from --fq-dir into <OUT>/<BioProject>_<readtype>_<platform>/00_fq/
     Expected names ONLY:
       PAIRED: {acc}_1.fastq(.gz), {acc}_2.fastq(.gz)
       SINGLE: {acc}.fastq(.gz)
     If an accession from metadata has no files in --fq-dir -> write it to
       <OUT>/<BioProject>_<readtype>_<platform>/missing_sra.list
  2) Per BioProject_<readtype>_<platform>:
       - Detect "3-file anomaly": an accession has BOTH pair (_1/_2) AND single (*.fastq(.gz))
         * append accession to sra_3_fq.note
         * move the "third" single file to sra_3_fq/
       - No marker files are created.
  3) Print a summary line per project dir:
       [QC] <dir>: PE=<n_pe>, SE=<n_se>, ANOM=<n_anom>

Platform buckets (case-insensitive):
  - illumina: illumina, hiseq, miseq, novaseq, nextseq
  - bgi:      bgi, bgiseq, mgiseq, mgitech, dnbseq
  - roche454: roche, 454
  - iontorrent: ion, torrent
  - unknown: everything else

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
    """Return the 12 fixed-form candidates (if exist)."""
    cands = [
        fq_dir / f"{acc}_1.fastq",
        fq_dir / f"{acc}_1.fastq.gz",
        fq_dir / f"{acc}_1.fq",
        fq_dir / f"{acc}_1.fq.gz",
        fq_dir / f"{acc}_2.fastq",
        fq_dir / f"{acc}_2.fastq.gz",
        fq_dir / f"{acc}_2.fq",
        fq_dir / f"{acc}_2.fq.gz",
        fq_dir / f"{acc}.fastq",
        fq_dir / f"{acc}.fastq.gz",
        fq_dir / f"{acc}.fq",
        fq_dir / f"{acc}.fq.gz",
    ]
    return [p for p in cands if p.exists() and p.is_file()]


def is_pe_by_files(files: List[Path]) -> bool:
    """Infer PE if any of the files carries _1 or _2."""
    names = [p.name for p in files]
    return any(n.endswith("_1.fastq") or n.endswith("_1.fastq.gz") or 
               n.endswith("_1.fq") or n.endswith("_1.fq.gz") or
               n.endswith("_2.fastq") or n.endswith("_2.fastq.gz") or 
               n.endswith("_2.fq") or n.endswith("_2.fq.gz") for n in names)


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
            pf  = row[29].strip() if len(row) > 29 else ""
            if not run or not bp:
                continue
            yield run, bp, pf


# ---------- Per-project anomaly analysis (no markers) ----------

def analyze_project(proj_dir: Path, dry_run: bool) -> Tuple[int, int, int]:
    """
    Scan 00_fq and:
      - Count PE / SE (for logging only)
      - Handle 3-file anomalies: move the single file to sra_3_fq/ and append note
    No marker files are created.
    """
    fqdir = proj_dir / "00_fq"
    if not fqdir.exists():
        return (0, 0, 0)

    seen: Dict[str, Dict[str, Path]] = defaultdict(dict)
    fq_files = list(fqdir.glob("*.fastq*")) + list(fqdir.glob("*.fq*"))
    for p in fq_files:
        name = p.name
        if name.endswith(".fastq.gz"):
            stem = name[:-9]
        elif name.endswith(".fastq"):
            stem = name[:-6]
        elif name.endswith(".fq.gz"):
            stem = name[:-6]
        elif name.endswith(".fq"):
            stem = name[:-3]    
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

    for acc, parts in seen.items():
        has_pair = ("r1" in parts) and ("r2" in parts)
        has_se   = ("se" in parts)
        if has_pair and not has_se:
            n_pe += 1
        elif has_se and not has_pair:
            n_se += 1
        elif has_pair and has_se:
            n_anom += 1
            if not dry_run:
                anom_dir.mkdir(exist_ok=True, parents=True)
                with anom_note.open("a") as f:
                    f.write(acc + "\n")
                try:
                    se_path = parts["se"]
                    dst = anom_dir / se_path.name
                    if not dst.exists():
                        shutil.move(str(se_path), str(dst))
                except Exception as e:
                    print(f"Warning: failed moving third file for {acc}: {e}", file=sys.stderr)
            else:
                print(f"DRY-RUN: anomaly {acc} -> would note and move")

    return (n_pe, n_se, n_anom)


# ---------- Main ----------

def main():
    ap = argparse.ArgumentParser(
        description="Bioproject-first FASTQ sorter with PE/SE + platform split."
    )
    ap.add_argument("-m", "--metadata", required=True, help="Metadata TSV path")
    ap.add_argument("-f", "--fq-dir", default="fq", help="Source FASTQ directory")
    ap.add_argument("-o", "--out-root", default=".", help="Output root")
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

    acc2bp: Dict[str, str] = {}
    acc2plat: Dict[str, str] = {}

    rows = list(iter_rows(meta, skip_header=args.header))
    for run, bp, pf in rows:
        acc2bp[run] = bp
        acc2plat[run] = classify_platform_bucket(pf)

    move_tasks: List[Tuple[Path, Path]] = []
    missing_by_proj: Dict[str, List[str]] = defaultdict(list)

    with tqdm(total=len(acc2bp), desc="Indexing", unit="acc") as pbar:
        for acc, bp in acc2bp.items():
            files = expected_paths(src, acc)
            plat = acc2plat.get(acc, "unknown")
            if not files:
                # Unknown readtype at this moment; keep previous behavior:
                # record missing under both pe and se groups for this platform.
                for rt in ["pe", "se"]:
                    proj_key = f"{bp}_{rt}_{plat}"
                    missing_by_proj[proj_key].append(acc)
            else:
                readtype = "pe" if is_pe_by_files(files) else "se"
                proj_dir = out / f"{bp}_{readtype}_{plat}" / "00_fq"
                for p in files:
                    move_tasks.append((p, proj_dir))
            pbar.update(1)

    if move_tasks:
        with ThreadPoolExecutor(max_workers=max(1, args.threads)) as ex:
            futs = [ex.submit(move_file, s, d, args.dry_run) for (s, d) in move_tasks]
            for _ in tqdm(as_completed(futs), total=len(futs), desc="Moving", unit="file"):
                pass

    # Write missing lists
    for proj_key, miss in missing_by_proj.items():
        proj_dir = out / proj_key
        if miss:
            path = proj_dir / "missing_sra.list"
            if args.dry_run:
                print(f"DRY-RUN: would write {path} with {len(miss)} accessions")
            else:
                proj_dir.mkdir(parents=True, exist_ok=True)
                with path.open("w") as f:
                    f.write("\n".join(sorted(set(miss))) + "\n")

    # Analyze anomalies per existing project dir (no markers, no platform notes)
    proj_dirs = sorted({ (out / f"{bp}_{rt}_{plat}")
                         for bp in set(acc2bp.values())
                         for rt in ["pe", "se"]
                         for plat in ["illumina","bgi","roche454","iontorrent","unknown"]
                         if (out / f"{bp}_{rt}_{plat}").exists() })

    for proj_dir in tqdm(proj_dirs, desc="Per-project QC", unit="proj"):
        n_pe, n_se, n_anom = analyze_project(proj_dir, args.dry_run)
        # Optional: you can print summary; nothing is written to disk.
        print(f"[QC] {proj_dir.name}: PE={n_pe}, SE={n_se}, ANOM={n_anom}")

    print("Done.")


if __name__ == "__main__":
    main()
