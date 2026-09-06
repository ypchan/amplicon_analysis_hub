#!/usr/bin/env python3
"""Safely normalize FASTQ suffixes and optionally gzip uncompressed inputs."""

from __future__ import annotations

import argparse
import gzip
import os
import shutil
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

VERSION = "2.0.0"


def transformed_name(path: Path, old_suffix: str, new_suffix: str, trim: tuple[str, ...]) -> str:
    name = path.name[: -len(old_suffix)] if old_suffix else path.name
    for token in trim:
        name = name.replace(token, "")
    if not name:
        raise ValueError(f"normalization creates an empty sample name for {path.name}")
    return name + new_suffix


def transfer(source: Path, destination: Path, dry_run: bool) -> str:
    if source.resolve() == destination.resolve():
        return f"unchanged\t{source}"
    if dry_run:
        return f"would_normalize\t{source}\t{destination}"
    if destination.exists():
        raise FileExistsError(f"destination already exists: {destination}")
    source_gz = source.name.endswith(".gz")
    destination_gz = destination.name.endswith(".gz")
    if source_gz == destination_gz:
        source.replace(destination)
    else:
        # Write atomically so an interrupted compression never leaves a valid-
        # looking partial destination.
        fd, temp_name = tempfile.mkstemp(prefix=f".{destination.name}.", dir=destination.parent)
        os.close(fd)
        temp_path = Path(temp_name)
        try:
            source_open = gzip.open if source_gz else Path.open
            destination_open = gzip.open if destination_gz else Path.open
            with source_open(source, "rb") as src, destination_open(temp_path, "wb") as dst:
                shutil.copyfileobj(src, dst, length=1024 * 1024)
            temp_path.replace(destination)
            source.unlink()
        except Exception:
            temp_path.unlink(missing_ok=True)
            raise
    return f"normalized\t{source}\t{destination}"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-i", "--input-dir", "--indir", dest="input_dir", required=True, type=Path,
                        help="Flat directory containing FASTQs")
    parser.add_argument("-1", "--reads1-suffix-in", "--reads1_suffix_in", dest="r1_in",
                        help="Current R1 or SE suffix")
    parser.add_argument("-2", "--reads2-suffix-in", "--reads2_suffix_in", dest="r2_in",
                        help="Current R2 suffix; requires --reads1-suffix-in")
    parser.add_argument("-f", "--reads1-suffix-out", "--reads1_suffix_std", dest="r1_out",
                        default="_R1.fastq.gz", help="Normalized R1/SE suffix")
    parser.add_argument("-r", "--reads2-suffix-out", "--reads2_suffix_std", dest="r2_out",
                        default="_R2.fastq.gz", help="Normalized R2 suffix")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Concurrent file operations")
    parser.add_argument("-s", "--trim", "--trim_str", dest="trim", default="",
                        help="Comma-separated literal tokens removed from sample names")
    parser.add_argument("--dry-run", action="store_true", help="Validate and print without changing files")
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    args = parser.parse_args()
    if not args.r1_in:
        parser.error("--reads1-suffix-in is required")
    if args.threads < 1:
        parser.error("--threads must be >= 1")
    if not args.input_dir.is_dir():
        parser.error(f"input directory not found: {args.input_dir}")
    if args.r2_in and args.r1_in == args.r2_in:
        parser.error("R1 and R2 input suffixes must differ")
    return args


def main() -> int:
    args = parse_args()
    trim = tuple(token for token in args.trim.split(",") if token)
    r1_files = sorted(path for path in args.input_dir.iterdir() if path.is_file() and path.name.endswith(args.r1_in))
    if not r1_files:
        print(f"error: no files end with {args.r1_in!r}", file=sys.stderr)
        return 2

    tasks: list[tuple[Path, Path]] = []
    for r1 in r1_files:
        r1_dest = r1.with_name(transformed_name(r1, args.r1_in, args.r1_out, trim))
        tasks.append((r1, r1_dest))
        if args.r2_in:
            sample = r1.name[: -len(args.r1_in)]
            r2 = r1.with_name(sample + args.r2_in)
            if not r2.is_file():
                print(f"error: missing mate for {r1.name}: {r2.name}", file=sys.stderr)
                return 2
            r2_dest = r2.with_name(transformed_name(r2, args.r2_in, args.r2_out, trim))
            tasks.append((r2, r2_dest))

    destinations = [destination for _, destination in tasks]
    if len(destinations) != len(set(destinations)):
        print("error: multiple inputs map to the same output name", file=sys.stderr)
        return 2
    conflicts = [str(destination) for source, destination in tasks
                 if source.resolve() != destination.resolve() and destination.exists()]
    if conflicts:
        print("error: destination(s) already exist: " + ", ".join(conflicts), file=sys.stderr)
        return 2

    try:
        with ThreadPoolExecutor(max_workers=min(args.threads, len(tasks))) as executor:
            for result in executor.map(lambda pair: transfer(*pair, args.dry_run), tasks):
                print(result)
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(f"Finished: {len(tasks)} file(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
