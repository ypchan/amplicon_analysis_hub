#!/usr/bin/env python3
"""Summarize Cutadapt text reports into deterministic TSV tables."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

VERSION = "2.0.0"
INTEGER = r"([\d,]+)"


@dataclass(frozen=True)
class AdapterHit:
    name: str = "NA"
    sequence: str = "NA"
    count: int = 0
    strand: str = "NA"


def _integer(pattern: str, text: str, default: int = 0) -> int:
    match = re.search(pattern, text, flags=re.IGNORECASE)
    return int(match.group(1).replace(",", "")) if match else default


def _best_adapter(text: str, heading: str) -> AdapterHit:
    """Return the most frequently trimmed adapter in a report section."""
    if heading:
        marker = rf"===\s*{heading}:\s*Adapter\s+['\"]?([^'\"=\n]+)['\"]?\s*==="
    else:
        marker = r"===\s*Adapter\s+['\"]?([^'\"=\n]+)['\"]?\s*==="
    starts = list(re.finditer(marker, text, flags=re.IGNORECASE))
    hits: list[AdapterHit] = []
    for index, match in enumerate(starts):
        end = starts[index + 1].start() if index + 1 < len(starts) else len(text)
        block = text[match.end():end]
        sequence_match = re.search(r"Sequence:\s*([A-Z]+)", block, flags=re.IGNORECASE)
        count = _integer(rf"Trimmed:\s*{INTEGER}\s+times", block)
        rc_count = _integer(rf"Reverse-complemented:\s*{INTEGER}", block)
        hits.append(
            AdapterHit(
                name=match.group(1).strip(),
                sequence=sequence_match.group(1).upper() if sequence_match else "NA",
                count=max(count, rc_count),
                strand="-" if rc_count > count else "+",
            )
        )
    return max(hits, key=lambda hit: hit.count, default=AdapterHit())


def parse_report(path: Path, mode: str) -> dict[str, object]:
    text = path.read_text(encoding="utf-8", errors="replace")
    sample = path.name.removesuffix(".cutadapt.log")
    if mode == "PE":
        total = _integer(rf"Total read pairs processed:\s*{INTEGER}", text)
        r1 = _best_adapter(text, "First read")
        r2 = _best_adapter(text, "Second read")
        return {
            "sample": sample,
            "total_reads_or_pairs": total,
            "r1_name": r1.name,
            "r1_sequence": r1.sequence,
            "r1_count": r1.count,
            "r1_percent": round(100 * r1.count / total, 3) if total else 0.0,
            "r1_strand": r1.strand,
            "r2_name": r2.name,
            "r2_sequence": r2.sequence,
            "r2_count": r2.count,
            "r2_percent": round(100 * r2.count / total, 3) if total else 0.0,
            "r2_strand": r2.strand,
            "parse_status": "ok" if total else "missing_total",
        }
    total = _integer(rf"Total reads processed:\s*{INTEGER}", text)
    read = _best_adapter(text, "")
    return {
        "sample": sample,
        "total_reads_or_pairs": total,
        "read_name": read.name,
        "read_sequence": read.sequence,
        "read_count": read.count,
        "read_percent": round(100 * read.count / total, 3) if total else 0.0,
        "read_strand": read.strand,
        "parse_status": "ok" if total else "missing_total",
    }


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def summarize(rows: list[dict[str, object]], mode: str) -> list[dict[str, object]]:
    total = len(rows)
    read_labels = ("r1", "r2") if mode == "PE" else ("read",)
    result: list[dict[str, object]] = []
    for label in read_labels:
        counts = Counter(
            (str(row[f"{label}_name"]), str(row[f"{label}_sequence"]), str(row[f"{label}_strand"]))
            for row in rows
        )
        for (name, sequence, strand), count in sorted(counts.items(), key=lambda item: (-item[1], item[0])):
            result.append(
                {
                    "read": label.upper(),
                    "primer": name,
                    "sequence": sequence,
                    "strand": strand,
                    "samples": count,
                    "total_samples": total,
                    "sample_percent": round(100 * count / total, 3) if total else 0.0,
                }
            )
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-d", "--dir", required=True, type=Path, help="Directory containing *.cutadapt.log")
    parser.add_argument("-m", "--mode", required=True, type=str.upper, choices=("PE", "SE"), help="Read layout")
    parser.add_argument("-o", "--output-dir", type=Path, default=Path("."), help="Directory for summary TSVs")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Concurrent log parsers")
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("--threads must be >= 1")
    if not args.dir.is_dir():
        parser.error(f"--dir is not a directory: {args.dir}")
    return args


def main() -> int:
    args = parse_args()
    reports = sorted(args.dir.glob("*.cutadapt.log"))
    if not reports:
        print(f"error: no *.cutadapt.log files in {args.dir}", file=sys.stderr)
        return 2
    with ThreadPoolExecutor(max_workers=min(args.threads, len(reports))) as executor:
        rows = list(executor.map(lambda path: parse_report(path, args.mode), reports))
    rows.sort(key=lambda row: str(row["sample"]))
    detail_fields = list(rows[0])
    write_tsv(args.output_dir / "cutadapt_details.tsv", rows, detail_fields)
    summary_rows = summarize(rows, args.mode)
    write_tsv(
        args.output_dir / "cutadapt_summary.tsv",
        summary_rows,
        ["read", "primer", "sequence", "strand", "samples", "total_samples", "sample_percent"],
    )
    failed = sum(row["parse_status"] != "ok" for row in rows)
    print(f"Parsed {len(rows)} reports ({failed} incomplete); output: {args.output_dir}")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
