#!/usr/bin/env python3
"""
is_16s_amplicon - I/O-optimized, concurrent 16S amplicon checker with BLAST.

Now: --threads means threads PER JOB. Total threads ≈ threads * concurrent.

Features:
  - Pure pipeline (no temp files): seqkit head | seqkit fq2fa | blastn -query -
  - Multi-sample concurrency: process many samples in parallel
  - Streamed output: each sample prints as soon as it finishes
  - STDIN or argv inputs; stdout formats: table / tsv / csv
  - Optional --output writes streamed results to a file (tsv/csv)

Examples:
  ls fq/*.gz | python3 is_16s_amplicon.py - --threads 4 --concurrent 3 --nreads 500
  # Total threads ≈ 4 * 3 = 12
"""

import argparse
import csv
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, Any, List, Tuple

HEADER = ["sample_id", "bac_hits", "arch_hits", "total_hits", "total_percent", "is_16S"]


# ------------ I/O-optimized single-sample runner ------------
def run_blast_stream(
    fq_path: Path,
    db: str,
    nreads: int,
    identity: float,
    threads_per_job: int,
) -> Dict[str, Any]:
    """
    Pipeline:
      seqkit head (sample N reads) ->
      seqkit fq2fa (FASTQ->FASTA) ->
      blastn -query - (read from stdin), outfmt 6
    Parse BLAST lines on the fly, count unique qseqid passing identity.
    """
    sample_name = fq_path.name

    # p1: sample N reads
    p1 = subprocess.Popen(
        ["seqkit", "head", "-j", str(max(1, threads_per_job)), "-n", str(nreads), str(fq_path)],
        stdout=subprocess.PIPE,
    )
    # p2: fastq -> fasta
    p2 = subprocess.Popen(
        ["seqkit", "fq2fa", "-j", str(max(1, threads_per_job)), "-"],
        stdin=p1.stdout,
        stdout=subprocess.PIPE,
    )
    if p1.stdout is not None:
        p1.stdout.close()

    # p3: blastn reads fasta from stdin
    p3 = subprocess.Popen(
        [
            "blastn", "-query", "-", "-db", db,
            "-evalue", "1e-5",
            "-outfmt", "6 qseqid sseqid pident length qlen slen",
            "-num_threads", str(max(1, threads_per_job)),
            "-max_target_seqs", "5",     # keep 5 to avoid BLAST warning
            "-max_hsps", "1",
            "-task", "megablast",
            "-word_size", "28",
            "-dust", "no",
        ],
        stdin=p2.stdout,
        stdout=subprocess.PIPE,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    if p2.stdout is not None:
        p2.stdout.close()

    # Stream-parse BLAST output
    hits: Dict[str, str] = {}
    assert p3.stdout is not None
    for line in p3.stdout:
        parts = line.rstrip("\n").split("\t")   # qseqid sseqid pident length qlen slen
        if len(parts) < 3:
            continue
        try:
            pid = float(parts[2])
        except ValueError:
            continue
        if pid >= identity:
            qid, sid = parts[0], parts[1]
            if qid not in hits:                # count each read once
                hits[qid] = sid

    p3.stdout.close()

    # Ensure all subprocesses ended successfully
    for proc in (p1, p2, p3):
        proc.wait()
        if proc.returncode not in (0, None):
            raise RuntimeError(f"Subprocess failed: {proc.args} (rc={proc.returncode})")

    total_hits = len(hits)
    bac_hits = sum(1 for sid in hits.values() if str(sid).startswith("bacteria__"))
    arch_hits = sum(1 for sid in hits.values() if str(sid).startswith("archaea__"))
    percent_hits = 100.0 * total_hits / max(1, nreads)
    is_amplicon = "YES" if percent_hits >= 50.0 else "NO"

    return {
        "sample_id": sample_name,
        "bac_hits": bac_hits,
        "arch_hits": arch_hits,
        "total_hits": total_hits,
        "total_percent": f"{percent_hits:.1f}",
        "is_16S": is_amplicon,
    }


# ---------------- Input & formatting helpers ----------------
def gather_inputs(argv_inputs: List[str]) -> List[str]:
    """Allow mixing '-' (stdin) and explicit file paths."""
    files: List[str] = []
    for tok in argv_inputs:
        if tok == "-":
            for line in sys.stdin:
                line = line.strip()
                if line:
                    files.append(line)
        else:
            files.append(tok)
    return files


def make_table_formatter(inputs: List[str]) -> Tuple[str, dict]:
    """Estimate widths for a neat streaming table."""
    widths = {h: len(h) for h in HEADER}
    if inputs:
        max_sample_len = max(len(Path(f).name) for f in inputs)
        widths["sample_id"] = max(widths["sample_id"], max_sample_len)
    fmt = "  ".join("{%s:%ds}" % (h, widths[h]) for h in HEADER)
    return fmt, widths


# ------------------------------- Main -------------------------------
def main():
    ap = argparse.ArgumentParser(description="I/O-optimized, concurrent 16S amplicon checker (BLAST).")
    ap.add_argument("inputs", nargs="+", help="FASTQ paths; use '-' to read from stdin (can be mixed)")
    ap.add_argument("-d", "--db", default="/share/cn1_fs/database/dada2_gtdb_ref/arch_bac_nr_16s",
                    help="BLAST database prefix")
    ap.add_argument("-n", "--nreads", type=int, default=1000,
                    help="Number of reads to sample per file (default: 1000)")
    ap.add_argument("-p", "--identity", type=float, default=60.0,
                    help="Identity cutoff percent (default: 60)")
    ap.add_argument("-t", "--threads", type=int, default=4,
                    help="Threads PER JOB (per sample). Total threads ≈ threads * concurrent")
    ap.add_argument("-c", "--concurrent", type=int, default=1,
                    help="How many samples to process concurrently (default: 1)")
    ap.add_argument("--format", choices=["table", "tsv", "csv"], default="table",
                    help="STDOUT format (default: table)")
    ap.add_argument("-o", "--output", default=None,
                    help="Also write streamed results to this file (tsv/csv)")
    ap.add_argument("--out-format", choices=["tsv", "csv"], default=None,
                    help="Force output file format (defaults by extension)")
    args = ap.parse_args()

    files = gather_inputs(args.inputs)
    if not files:
        sys.exit("No input files.")

    # Concurrency semantics: threads_per_job = args.threads (per job)
    concurrent = max(1, args.concurrent)
    threads_per_job = max(1, args.threads)  # per job threads; total ≈ threads_per_job * concurrent

    # stdout writer (streaming)
    out_writer = None
    fmt = None
    if args.format == "tsv":
        out_writer = csv.DictWriter(sys.stdout, fieldnames=HEADER, delimiter="\t", lineterminator="\n")
        out_writer.writeheader(); sys.stdout.flush()
    elif args.format == "csv":
        out_writer = csv.DictWriter(sys.stdout, fieldnames=HEADER, lineterminator="\n")
        out_writer.writeheader(); sys.stdout.flush()
    else:  # table
        fmt, _ = make_table_formatter(files)
        print(fmt.format(**{h: h for h in HEADER}), flush=True)

    # optional file writer (streaming)
    file_writer = None
    f_handle = None
    if args.output:
        out_path = Path(args.output)
        file_fmt = args.out_format if args.out_format else ("tsv" if out_path.suffix.lower() == ".tsv" else "csv")
        out_path.parent.mkdir(parents=True, exist_ok=True)
        f_handle = out_path.open("w", newline="")
        if file_fmt == "tsv":
            file_writer = csv.DictWriter(f_handle, fieldnames=HEADER, delimiter="\t", lineterminator="\n")
        else:
            file_writer = csv.DictWriter(f_handle, fieldnames=HEADER, lineterminator="\n")
        file_writer.writeheader()

    # run concurrently; print each result as soon as it finishes
    try:
        with ThreadPoolExecutor(max_workers=concurrent) as ex:
            fut2file = {
                ex.submit(
                    run_blast_stream,
                    Path(fq),
                    args.db,
                    args.nreads,
                    args.identity,
                    threads_per_job,
                ): fq
                for fq in files
            }

            for fut in as_completed(fut2file):
                res = fut.result()

                # stdout
                if args.format == "table":
                    print(fmt.format(**{h: str(res[h]) for h in HEADER}), flush=True)
                else:
                    out_writer.writerow(res); sys.stdout.flush()

                # optional file
                if file_writer:
                    file_writer.writerow(res)
                    f_handle.flush()
    finally:
        if f_handle:
            f_handle.close()


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(line_buffering=True)
    except Exception:
        pass
    main()
