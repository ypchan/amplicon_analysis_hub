"""
fastq_sorter.py
Efficient FASTQ sorter for fixed filename patterns, with parallel moves & tqdm progress.

Supported filename forms ONLY (case-sensitive):
  PAIRED:  {accession}_1.fastq , {accession}_1.fastq.gz ,
           {accession}_2.fastq , {accession}_2.fastq.gz
  SINGLE:  {accession}.fastq   , {accession}.fastq.gz

Metadata: TSV without header (default 1-based columns: Run=2, BioProject=19, LibraryLayout=29, Platform=30)
"""

import argparse
import csv
import os
import shutil
import sys
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import List, Tuple

from tqdm import tqdm  # pip install tqdm


def classify_platform(p: str) -> str:
    """Map platform string to one of the standard categories."""
    p = (p or "").lower()
    if "illumina" in p:
        return "Illumina"
    if "roche" in p or "454" in p:
        return "Roche_454"
    if "ion" in p or "torrent" in p:
        return "Ion_Torrent"
    return "Unknown"


def expected_files(fq_dir: Path, accession: str, layout: str) -> List[Path]:
    """
    Return the list of expected file paths based on accession and layout.
    Only the 6 allowed filename forms are checked.
    """
    if layout == "PAIRED":
        candidates = [
            fq_dir / f"{accession}_1.fastq",
            fq_dir / f"{accession}_1.fastq.gz",
            fq_dir / f"{accession}_2.fastq",
            fq_dir / f"{accession}_2.fastq.gz",
        ]
    else:  # SINGLE
        candidates = [
            fq_dir / f"{accession}.fastq",
            fq_dir / f"{accession}.fastq.gz",
        ]
    # Keep only existing regular files
    return [p for p in candidates if p.exists() and p.is_file()]


def read_metadata_rows(path: Path, cols: Tuple[int, int, int, int]):
    """
    Yield tuples (run, bioproject, layout, platform) from a no-header TSV.
    Column indices are 1-based.
    """
    run_i, bp_i, ll_i, pf_i = cols
    with path.open("r", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if len(row) < max(cols):
                continue
            run = row[run_i - 1].strip()
            if not run:
                continue
            bp = (row[bp_i - 1].strip() or "NA")
            ll = row[ll_i - 1].strip().upper()
            pf = row[pf_i - 1].strip()
            yield run, bp, ll, pf


def move_one(src: Path, dest_dir: Path, dry_run: bool) -> bool:
    """
    Move a file into the destination directory.
    Returns True if moved (or dry-run), False if skipped (e.g., file exists).
    """
    dest_dir.mkdir(parents=True, exist_ok=True)
    dst = dest_dir / src.name
    if dry_run:
        print(f"DRY-RUN: mv -n '{src}' '{dest_dir}/'")
        return True
    if dst.exists():  # do not overwrite
        return False
    shutil.move(str(src), str(dst))
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-m", "--metadata", required=True, help="Metadata TSV (no header)")
    ap.add_argument("-f", "--fq-dir", default="fq", help="FASTQ directory (no recursion)")
    ap.add_argument("-r", "--report", default="fq_sorting_report.tsv", help="Output report TSV")
    ap.add_argument("-o", "--out-root", default=".", help="Output root (creates PAIRED/ and SINGLE/)")
    ap.add_argument("--dry-run", action="store_true", help="Show actions without moving files")
    ap.add_argument("-t", "--threads", type=int, default=min(16, (os.cpu_count() or 8)), help="Parallel threads")
    # 1-based column indices
    ap.add_argument("--col-run", type=int, default=2, help="1-based column index for Run/Accession")
    ap.add_argument("--col-bioproject", type=int, default=19, help="1-based column index for BioProject")
    ap.add_argument("--col-layout", type=int, default=29, help="1-based column index for LibraryLayout")
    ap.add_argument("--col-platform", type=int, default=30, help="1-based column index for Platform")
    args = ap.parse_args()

    meta = Path(args.metadata)
    fq_dir = Path(args.fq_dir)
    out_root = Path(args.out_root)
    if not meta.exists():
        sys.exit(f"Error: metadata not found: {meta}")
    if not fq_dir.exists():
        sys.exit(f"Error: FASTQ dir not found: {fq_dir}")

    # Create top-level output directories
    (out_root / "PAIRED").mkdir(parents=True, exist_ok=True)
    (out_root / "SINGLE").mkdir(parents=True, exist_ok=True)

    # Phase 1: Build move task list by checking for expected files
    tasks: List[Tuple[Path, Path, str, str, str]] = []  # (src, dest_dir, layout, bioproject, platformCategory)
    missing_accessions = []
    rows = list(read_metadata_rows(meta, (args.col_run, args.col_bioproject, args.col_layout, args.col_platform)))

    # Deduplicate by run; keep the last occurrence
    last_seen = {}
    for run, bp, ll, pf in rows:
        last_seen[run] = (run, bp, ll, pf)

    unique_rows = list(last_seen.values())
    with tqdm(total=len(unique_rows), desc="Indexing", unit="rec") as pbar:
        for run, bp, ll, pf in unique_rows:
            if ll not in {"SINGLE", "PAIRED"}:
                pbar.update(1)
                continue
            plat_dir = classify_platform(pf)
            dest_dir = out_root / ll / plat_dir / bp / "00_fq"
            files = expected_files(fq_dir, run, ll)
            if not files:
                missing_accessions.append((run, bp, ll))
            else:
                for f in files:
                    tasks.append((f, dest_dir, ll, bp, plat_dir))
            pbar.update(1)

    if not tasks:
        # Still write an empty report and exit
        with open(args.report, "w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
            writer.writerow(["LibLayout", "BioProject", "Platform", "FASTQ_Count"])
        print("Nothing to move. Check accession names and files present.")
        if missing_accessions:
            print(f"Accessions with no matching files: {len(missing_accessions)}")
        sys.exit(0)

    # Phase 2: Parallel moves
    counts = Counter()
    moved_total = 0
    with ThreadPoolExecutor(max_workers=max(1, args.threads)) as ex:
        futs = [ex.submit(move_one, src, dst, args.dry_run) for (src, dst, _, _, _) in tasks]
        for (src, dst, ll, bp, plat), fut in tqdm(zip(tasks, as_completed(futs)),
                                                  total=len(tasks), desc="Moving", unit="file"):
            ok = fut.result()
            if ok:
                counts[(ll, bp, plat)] += 1
                moved_total += 1

    # Phase 3: Write report
    with open(args.report, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["LibLayout", "BioProject", "Platform", "FASTQ_Count"])
        for (ll, bp, plat), cnt in sorted(counts.items()):
            writer.writerow([ll, bp, plat, cnt])

    print(f"Completed. Moved files: {moved_total}/{len(tasks)}. Report: {args.report}")
    if args.dry_run:
        print("Note: dry-run mode; no files were moved.")


if __name__ == "__main__":
    main()
