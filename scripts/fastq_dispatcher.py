#!/usr/bin/env python3
# date: 2024-08-30 latest update (patched: 2026-01-23)
# contact: yanpengch@qq.com
"""
FASTQ dispatcher (drive by FASTQs present), Bioproject-first, PE/SE + platform split.

Patch 2026-01-23:
  (A) Avoid empty project dirs:
      - Drive by FASTQs present in --fq_dir (target ~30k), not all metadata runs (~300k).
      - Missing summarized per BioProject in ONE report file (no mkdir for missing-only projects).

  (B) Keep platform_raw and output unknown platforms list.

  (C) Platform directory name uses *platform raw info* but canonicalizes common ones:
      - ion_torrent / iontorrent -> iontorrent
      - pacbio_smrt / pacbio ccs / smrt / ccs / pacbio -> pacbio
      - ls454 / roche454 / 454 -> roche454
      - oxford_nanopore / ont / nanopore -> nanopore
    Then:
      - lowercase
      - remove spaces
      - sanitize to [a-z0-9._-] (others -> '_')
"""

import argparse
import csv
import os
import re
import shutil
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, Iterable, List, Tuple, Optional

try:
    from tqdm import tqdm
except ImportError:  # tqdm is cosmetic; keep the dispatcher usable without it.
    def tqdm(iterable, **_kwargs):
        return iterable


# ---------------- Layout parsing ----------------

def parse_layout(layout_raw: str) -> Optional[str]:
    """
    Parse layout field to 'pe' or 'se' if possible.
    Returns None if unknown/unparseable.
    """
    s = (layout_raw or "").strip().lower()
    if not s:
        return None
    if "paired" in s or s in ("pe", "pair"):
        return "pe"
    if "single" in s or s in ("se", "single_end", "singleend"):
        return "se"
    return None


# ---------------- Platform normalization for directory name ----------------

_SAN_RE = re.compile(r"[^a-z0-9._-]+")  # anything not safe -> '_'
_WS_RE = re.compile(r"\s+")

def normalize_platform_dir(platform_raw: str) -> str:
    """
    Directory tag for platform:
      1) canonicalize common platforms you care about:
         - iontorrent, pacbio, roche454, nanopore
      2) otherwise keep platform raw (but normalized):
         - lowercase, no spaces, sanitize to [a-z0-9._-]
    """
    s0 = (platform_raw or "").strip()
    if not s0:
        return "unknown"

    # keep a "matching string" with spaces/underscores/hyphens removed
    s = s0.lower()
    s_compact = re.sub(r"[\s_\-]+", "", s)


    # --- canonical buckets (order matters) ---
    # illumina : hiseq/novaseq/nextseq/miseq -> illumina
    if any(k in s_compact for k in ("illumina", "hiseq", "miseq", "novaseq", "nextseq")):
        return "illumina"

    # bgi (optional but common)
    if any(k in s_compact for k in ("bgi", "bgiseq", "mgiseq", "mgitech", "dnbseq")):
        return "bgi"

    # nanopore
    if ("oxfordnanopore" in s_compact) or ("nanopore" in s_compact) or (s_compact == "ont") or ("promethion" in s_compact) or ("minion" in s_compact):
        return "nanopore"

    # pacbio
    if ("pacbio" in s_compact) or ("smrt" in s_compact) or ("ccs" in s_compact) or ("hifi" in s_compact) or ("revi" in s_compact):
        return "pacbio"

    # roche454 / ls454
    if ("roche" in s_compact) or ("ls454" in s_compact) or (s_compact == "454") or ("gsflx" in s_compact) or ("titanium" in s_compact):
        return "roche454"

    # iontorrent
    if ("iontorrent" in s_compact) or ("iontorrent" in s.replace(" ", "")) or ("ion" in s_compact and "torrent" in s_compact) or ("pgm" in s_compact) or ("proton" in s_compact) or ("s5" in s_compact):
        return "iontorrent"

    # --- fallback: normalize raw string ---
    s = _WS_RE.sub("", s)       # remove all whitespace
    s = _SAN_RE.sub("_", s)     # sanitize
    s = s.strip("_")
    return s if s else "unknown"


# ---------------- FASTQ scanning (drive by files, not metadata) ----------------

FASTQ_RE = re.compile(r"""
    ^
    (?P<acc>[^/]+?)              # accession
    (?:
        _(?P<mate>[12])          # optional _1/_2
    )?
    \.(?P<ext>fastq|fq)
    (?P<gz>\.gz)?
    $
""", re.VERBOSE | re.IGNORECASE)

def scan_fastq_dir(fq_dir: Path) -> Dict[str, Dict[str, List[Path]]]:
    """
    Return mapping:
      acc -> { 'r1': [...], 'r2': [...], 'se': [...] }
    Only files matching fixed naming patterns are considered.
    """
    acc_files: Dict[str, Dict[str, List[Path]]] = defaultdict(lambda: {"r1": [], "r2": [], "se": []})

    for p in fq_dir.iterdir():
        if not p.is_file():
            continue
        m = FASTQ_RE.match(p.name)
        if not m:
            continue
        acc = m.group("acc")
        mate = m.group("mate")
        if mate == "1":
            acc_files[acc]["r1"].append(p)
        elif mate == "2":
            acc_files[acc]["r2"].append(p)
        else:
            acc_files[acc]["se"].append(p)

    return acc_files


# ---------------- Moving helpers ----------------

def dispatch_file(src: Path, dst_dir: Path, action: str, dry_run: bool) -> Tuple[Path, str]:
    """Move, copy, or symlink one FASTQ without overwriting a destination."""
    dst = dst_dir / src.name
    if dry_run:
        return dst, f"would_{action}"
    if dst.exists():
        return dst, "exists"
    dst_dir.mkdir(parents=True, exist_ok=True)
    if action == "move":
        shutil.move(str(src), str(dst))
    elif action == "copy":
        shutil.copy2(src, dst)
    else:
        dst.symlink_to(src.resolve())
    return dst, action

def append_note(note_path: Path, acc: str, dry_run: bool) -> None:
    if dry_run:
        print(f"DRY-RUN: would append '{acc}' to {note_path}")
        return
    note_path.parent.mkdir(parents=True, exist_ok=True)
    with note_path.open("a") as f:
        f.write(acc + "\n")


# ---------------- Main ----------------

def main():
    ap = argparse.ArgumentParser(
        description="FASTQ dispatcher driven by files present: BioProject, layout, and platform split.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=r"""
Example:
  fastq_dispatcher.py -m meta.tsv -f fq -o out -t 16 --header --action move
"""
    )
    ap.add_argument("-m", "--metadata", required=True, metavar="TSV", help="Metadata TSV path")
    ap.add_argument("-f", "--fq-dir", "--fq_dir", dest="fq_dir", required=True, metavar="DIR",
                    help="Source FASTQ directory (flat)")
    ap.add_argument("-o", "--out-dir", "--out_dir", "--out-root", dest="out_dir", default=".",
                    metavar="DIR", help="Output root")
    ap.add_argument("-t", "--threads", type=int, default=min(16, os.cpu_count() or 4), metavar="INT",
                    help="Concurrent file operations")
    ap.add_argument("--header", "--header-1strow", "--header_1strow", dest="header_1strow",
                    action="store_true", help="Skip the metadata header row")

    # column controls (1-based) with defaults
    ap.add_argument("--run-col", "--run_col", dest="run_col", type=int, default=1, metavar="INT",
                    help="1-based column index for run accession")
    ap.add_argument("--bioproject-col", "--bioproject_col", dest="bioproject_col", type=int, default=22, metavar="INT",
                    help="1-based column index for BioProject")
    ap.add_argument("--layout-col", "--layout_col", dest="layout_col", type=int, default=16, metavar="INT",
                    help="1-based column index for layout; 0 disables metadata layout")
    ap.add_argument("--platform-col", "--platform_col", dest="platform_col", type=int, default=19, metavar="INT",
                    help="1-based column index for platform; 0 disables metadata platform")

    ap.add_argument("--action", choices=("move", "copy", "symlink"), default="move",
                    help="How FASTQs are dispatched; move changes the source directory")
    ap.add_argument("--dry-run", action="store_true", help="Report actions without changing files")
    args = ap.parse_args()

    meta = Path(args.metadata)
    fq_dir = Path(args.fq_dir)
    out = Path(args.out_dir)

    if not meta.is_file():
        sys.exit(f"Error: metadata not found: {meta}")
    if not fq_dir.is_dir():
        sys.exit(f"Error: FASTQ dir not found: {fq_dir}")
    if args.threads < 1:
        ap.error("--threads must be >= 1")
    if args.run_col < 1 or args.bioproject_col < 1 or args.layout_col < 0 or args.platform_col < 0:
        ap.error("column indices must be 1-based positive integers; layout/platform may be 0")
    if not args.dry_run:
        out.mkdir(parents=True, exist_ok=True)

    # 1) Scan FASTQs -> target accessions
    print(f"[1/4] Scanning FASTQ directory: {fq_dir}")
    acc_files = scan_fastq_dir(fq_dir)
    target_accs = set(acc_files.keys())
    print(f"      Found accessions with FASTQs: {len(target_accs)}")

    if not target_accs:
        print("No FASTQ accessions found. Exit.")
        return

    # 2) Read metadata stream: only keep rows whose run is in target_accs
    print(f"[2/4] Reading metadata (only for FASTQ-present runs): {meta}")
    acc2bp: Dict[str, str] = {}
    acc2layout: Dict[str, Optional[str]] = {}
    acc2plat_raw: Dict[str, str] = {}
    acc2plat_dir: Dict[str, str] = {}

    # Total runs per bioproject in metadata (for missing stats)
    bp_total_in_meta: Dict[str, int] = defaultdict(int)

    run_i = args.run_col - 1
    bp_i = args.bioproject_col - 1
    lay_i = (args.layout_col - 1) if args.layout_col and args.layout_col > 0 else None
    pf_i  = (args.platform_col - 1) if args.platform_col and args.platform_col > 0 else None

    with meta.open("r", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        if args.header_1strow:
            next(reader, None)

        for row in tqdm(reader, desc="      metadata rows", unit="row"):
            if len(row) <= max(run_i, bp_i):
                continue
            run = (row[run_i] or "").strip()
            bp  = (row[bp_i] or "").strip()
            if not run or not bp:
                continue

            bp_total_in_meta[bp] += 1

            if run not in target_accs:
                continue

            layout_raw = (row[lay_i].strip() if (lay_i is not None and len(row) > lay_i and row[lay_i]) else "")
            platform_raw = (row[pf_i].strip() if (pf_i is not None and len(row) > pf_i and row[pf_i]) else "")

            acc2bp[run] = bp
            acc2layout[run] = parse_layout(layout_raw)
            acc2plat_raw[run] = platform_raw
            acc2plat_dir[run] = normalize_platform_dir(platform_raw)

    # FASTQ accessions missing in metadata
    accs_without_meta = sorted(target_accs - set(acc2bp.keys()))
    if accs_without_meta:
        for acc in accs_without_meta:
            acc2bp[acc] = "NA"
            acc2layout[acc] = None
            acc2plat_raw[acc] = ""
            acc2plat_dir[acc] = "unknown"

    # 3) Build move tasks (only FASTQ-present accs)
    print(f"[3/4] Dispatching FASTQs into project dirs (no empty dirs)")
    move_tasks: List[Tuple[Path, Path]] = []
    bp_present_fastq: Dict[str, int] = defaultdict(int)

    # record unknown platforms (raw) for later manual mapping
    unknown_platform_rows: List[Tuple[str, str, str]] = []

    for acc in tqdm(sorted(target_accs), desc="      indexing accessions", unit="acc"):
        bp = acc2bp.get(acc, "NA")
        layout_hint = acc2layout.get(acc)  # pe/se/None
        plat_raw = acc2plat_raw.get(acc, "")
        plat_dir = acc2plat_dir.get(acc, "unknown")

        if plat_dir == "unknown" and plat_raw:
            praw = plat_raw.replace("\t", " ").replace("\n", " ").strip()
            unknown_platform_rows.append((acc, bp, praw))

        files = acc_files[acc]
        r1 = files["r1"]
        r2 = files["r2"]
        se = files["se"]

        has_pair = bool(r1) and bool(r2)
        has_se = bool(se)

        # Decide readtype
        if has_pair:
            readtype = "pe"
        elif (r1 or r2) and layout_hint == "pe":
            readtype = "pe"
        elif layout_hint == "se":
            readtype = "se"
        else:
            readtype = "pe" if (r1 or r2) else "se"

        proj_base = out / f"{bp}_{readtype}_{plat_dir}"
        fq_out = proj_base / "00_fq"

        bp_present_fastq[bp] += 1

        # 3-file anomaly: keep pair, discard se
        if has_pair and has_se:
            append_note(proj_base / "sra_3_fq.note", acc, args.dry_run)
            discard_dir = proj_base / "sra_3_fq_discarded"
            for p in r1 + r2:
                move_tasks.append((p, fq_out))
            for p in se:
                move_tasks.append((p, discard_dir))
        else:
            for p in (r1 + r2 + se):
                move_tasks.append((p, fq_out))

    # Execute moves
    collisions: List[Path] = []
    if move_tasks:
        with ThreadPoolExecutor(max_workers=max(1, args.threads)) as ex:
            futs = [ex.submit(dispatch_file, s, d, args.action, args.dry_run) for (s, d) in move_tasks]
            for future in tqdm(as_completed(futs), total=len(futs), desc=f"      {args.action} files", unit="file"):
                destination, status = future.result()
                if status == "exists":
                    collisions.append(destination)

    # 4) Reports
    print(f"[4/4] Writing reports")

    report_bp = out / "bioproject_missing.summary.tsv"
    if args.dry_run:
        print(f"DRY-RUN: would write {report_bp}")
    else:
        with report_bp.open("w") as f:
            f.write("bioproject\tpresent_fastq_runs\ttotal_runs_in_metadata\tmissing_runs_estimate\n")
            for bp in sorted(bp_present_fastq.keys()):
                present = bp_present_fastq[bp]
                total = bp_total_in_meta.get(bp, 0)
                missing = max(total - present, 0) if total else "NA"
                f.write(f"{bp}\t{present}\t{total}\t{missing}\n")

    report_unknown = out / "platform_unknown.tsv"
    if unknown_platform_rows:
        if args.dry_run:
            print(f"DRY-RUN: would write {report_unknown} ({len(unknown_platform_rows)} rows)")
        else:
            with report_unknown.open("w") as f:
                f.write("run\tbioproject\tplatform_raw\n")
                for run, bp, praw in unknown_platform_rows:
                    f.write(f"{run}\t{bp}\t{praw}\n")

    if collisions:
        report_collision = out / "destination_collisions.tsv"
        if args.dry_run:
            print(f"DRY-RUN: {len(collisions)} existing destinations would be skipped")
        else:
            report_collision.write_text(
                "destination\n" + "".join(f"{path}\n" for path in sorted(collisions)),
                encoding="utf-8",
            )
        print(f"WARNING: skipped {len(collisions)} existing destinations", file=sys.stderr)

    print("Done.")


if __name__ == "__main__":
    main()
