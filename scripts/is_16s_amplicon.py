#!/usr/bin/env python3
"""Screen FASTQ files against the bundled bacterial/archaeal 16S BLAST DB."""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any, Iterable, TextIO

VERSION = "2.0.0"
HEADER = ["sample_id", "bac_hits", "arch_hits", "total_hits", "total_percent", "is_16S"]
DEFAULT_DB = Path(__file__).resolve().parent.parent / "data" / "arc_bac_16s_blastDB" / "arch_bac_16s_ref_90"


def open_text(path: Path) -> TextIO:
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def sample_fasta(path: Path, limit: int) -> tuple[str, int]:
    """Read at most *limit* FASTQ records and return FASTA plus actual count."""
    chunks: list[str] = []
    count = 0
    with open_text(path) as handle:
        while count < limit:
            header = handle.readline()
            if not header:
                break
            sequence = handle.readline().strip()
            plus = handle.readline()
            quality = handle.readline().strip()
            if not sequence or not plus or not quality:
                raise ValueError(f"truncated FASTQ record at read {count + 1}")
            if not header.startswith("@") or not plus.startswith("+"):
                raise ValueError(f"invalid FASTQ record at read {count + 1}")
            if len(sequence) != len(quality):
                raise ValueError(f"sequence/quality lengths differ at read {count + 1}")
            # Synthetic IDs guarantee uniqueness even in concatenated public
            # FASTQs whose original read identifiers are duplicated.
            read_id = f"read_{count + 1}"
            chunks.append(f">{read_id}\n{sequence}\n")
            count += 1
    return "".join(chunks), count


def database_exists(prefix: Path) -> bool:
    return any(
        candidate.exists()
        for candidate in (
            Path(f"{prefix}.nhr"), Path(f"{prefix}.ndb"), Path(f"{prefix}.00.nhr")
        )
    )


def run_screen(
    fastq: Path,
    db: Path,
    nreads: int,
    identity: float,
    query_coverage: float,
    hit_threshold: float,
    threads: int,
) -> dict[str, Any]:
    fasta, sampled = sample_fasta(fastq, nreads)
    if sampled == 0:
        return {
            "sample_id": fastq.name,
            "bac_hits": 0,
            "arch_hits": 0,
            "total_hits": 0,
            "total_percent": "0.000",
            "is_16S": "NO",
        }
    command = [
        "blastn", "-query", "-", "-db", str(db), "-task", "blastn",
        "-word_size", "11", "-evalue", "1e-10", "-dust", "no",
        "-num_threads", str(threads), "-max_target_seqs", "5", "-max_hsps", "1",
        "-outfmt", "6 qseqid sseqid pident qcovhsp bitscore",
    ]
    completed = subprocess.run(
        command,
        input=fasta,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if completed.returncode:
        raise RuntimeError(f"blastn failed for {fastq}: {completed.stderr.strip()}")

    hits: dict[str, str] = {}
    for line in completed.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 4:
            continue
        try:
            passed = float(fields[2]) >= identity and float(fields[3]) >= query_coverage
        except ValueError:
            continue
        if passed:
            hits.setdefault(fields[0], fields[1])
    bacterial = sum(subject.lower().startswith("bacteria_") for subject in hits.values())
    archaeal = sum(subject.lower().startswith("archaea_") for subject in hits.values())
    percent = 100.0 * len(hits) / sampled
    return {
        "sample_id": fastq.name,
        "bac_hits": bacterial,
        "arch_hits": archaeal,
        "total_hits": len(hits),
        "total_percent": f"{percent:.3f}",
        "is_16S": "YES" if percent >= hit_threshold else "NO",
    }


def gather_inputs(tokens: Iterable[str]) -> list[Path]:
    inputs: list[Path] = []
    for token in tokens:
        values = (line.strip() for line in sys.stdin) if token == "-" else (token,)
        inputs.extend(Path(value) for value in values if value)
    return list(dict.fromkeys(inputs))


def make_writer(handle: TextIO, output_format: str) -> csv.DictWriter:
    return csv.DictWriter(
        handle,
        fieldnames=HEADER,
        delimiter="\t" if output_format == "tsv" else ",",
        lineterminator="\n",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=(
            "Examples:\n"
            "  is_16s_amplicon.py sample.fastq.gz\n"
            "  printf 'a.fastq.gz\\nb.fastq.gz\\n' | is_16s_amplicon.py - -c 4 --format tsv\n"
            "\nInterpretation: this is a conservative content screen, not taxonomic assignment. "
            "For low-diversity or highly divergent environmental samples, inspect the hit table "
            "before excluding data."
        ),
    )
    parser.add_argument("inputs", nargs="+", help="FASTQ paths; '-' reads paths from stdin")
    parser.add_argument("-d", "--db", type=Path, default=DEFAULT_DB, help="BLAST database prefix")
    parser.add_argument("-n", "--nreads", type=int, default=100, help="Maximum reads sampled per FASTQ")
    parser.add_argument("-p", "--identity", type=float, default=70.0, help="Minimum nucleotide identity percent")
    parser.add_argument("-q", "--query-coverage", type=float, default=70.0, help="Minimum query coverage percent")
    parser.add_argument("--hit-threshold", type=float, default=50.0, help="Percent of sampled reads with hits required for YES")
    parser.add_argument("-t", "--threads", type=int, default=1, help="BLAST threads per FASTQ")
    parser.add_argument("-c", "--concurrent", type=int, default=1, help="FASTQ files processed concurrently")
    parser.add_argument("--format", choices=("table", "tsv", "csv"), default="table", help="Standard-output format")
    parser.add_argument("-o", "--output", type=Path, help="Optional TSV/CSV output path")
    parser.add_argument("--out-format", choices=("tsv", "csv"), help="Output-file format; inferred from suffix otherwise")
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    args = parser.parse_args()
    if args.nreads < 1 or args.threads < 1 or args.concurrent < 1:
        parser.error("--nreads, --threads, and --concurrent must be >= 1")
    for name in ("identity", "query_coverage", "hit_threshold"):
        if not 0 <= getattr(args, name) <= 100:
            parser.error(f"--{name.replace('_', '-')} must be between 0 and 100")
    return args


def main() -> int:
    args = parse_args()
    if shutil.which("blastn") is None:
        print("error: blastn is not installed or not in PATH", file=sys.stderr)
        return 127
    if not database_exists(args.db):
        print(f"error: BLAST database not found for prefix {args.db}", file=sys.stderr)
        return 2
    files = gather_inputs(args.inputs)
    if not files:
        print("error: no input paths were provided", file=sys.stderr)
        return 2
    missing = [str(path) for path in files if not path.is_file()]
    if missing:
        print("error: FASTQ file(s) not found: " + ", ".join(missing), file=sys.stderr)
        return 2

    worker = lambda path: run_screen(
        path, args.db, args.nreads, args.identity, args.query_coverage,
        args.hit_threshold, args.threads
    )
    with ThreadPoolExecutor(max_workers=min(args.concurrent, len(files))) as executor:
        results = list(executor.map(worker, files))

    file_handle: TextIO | None = None
    try:
        file_writer = None
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            file_handle = args.output.open("w", encoding="utf-8", newline="")
            file_format = args.out_format or ("tsv" if args.output.suffix.lower() == ".tsv" else "csv")
            file_writer = make_writer(file_handle, file_format)
            file_writer.writeheader()

        if args.format in {"tsv", "csv"}:
            stdout_writer = make_writer(sys.stdout, args.format)
            stdout_writer.writeheader()
            for result in results:
                stdout_writer.writerow(result)
                if file_writer:
                    file_writer.writerow(result)
        else:
            widths = {field: max(len(field), *(len(str(row[field])) for row in results)) for field in HEADER}
            template = "  ".join(f"{{{field}:<{widths[field]}}}" for field in HEADER)
            print(template.format(**{field: field for field in HEADER}))
            for result in results:
                print(template.format(**result))
                if file_writer:
                    file_writer.writerow(result)
    finally:
        if file_handle:
            file_handle.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
