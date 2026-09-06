#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# contact: yanpengch@qq.com
# date: 2025-09-12

import argparse
import csv
import sys
import time
import fileinput
import urllib.parse
import urllib.request
import ssl
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Tuple, Optional

# Try tqdm for a rich progress bar; provide a lightweight fallback if missing
try:
    from tqdm import tqdm as _tqdm  # type: ignore
except Exception:
    _tqdm = None  # fallback will be used


# ENA Portal API endpoint and fields to request
ENA_API = "https://www.ebi.ac.uk/ena/portal/api/filereport"
# Request fields needed to support fastq/sra/submitted preferences
FIELDS = ",".join([
    "run_accession",
    "fastq_ftp", "fastq_md5", "fastq_bytes",
    "sra_ftp", "sra_md5", "sra_bytes",
    "submitted_ftp", "submitted_md5", "submitted_bytes",
])

VERSION = "2.0.0"
UA = "amplicon-analysis-hub/2.0 (ENA Portal API client)"
TIMEOUT = 60
MAX_RETRIES = 5
RETRY_BACKOFF = 1.8  # Exponential backoff factor for retries


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Resolve ENA download links for large SRA accession lists.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("accession_file", help="one SRR/ERR/DRR per line")
    parser.add_argument(
        "--out-tsv",
        default="links.tsv",
        help="Detailed TSV (run_accession,kind,url,md5,bytes)",
    )
    parser.add_argument(
        "--missing",
        default="ENA_missing.tsv",
        help="Accessions with no links found",
    )
    parser.add_argument(
        "--prefer",
        default="fastq",
        choices=["fastq", "sra", "submitted", "any"],
        help="Preferred file type",
    )
    parser.add_argument(
        "--scheme",
        choices=["https", "ftp"],
        default="https",
        help="URL scheme prepended when ENA returns host/path without a scheme",
    )
    parser.add_argument(
        "--batch",
        type=int,
        default=100,
        help="Accessions per API request",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=8,
        help="Concurrent API requests",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    args = parser.parse_args()
    if args.batch < 1 or args.batch > 1000:
        parser.error("--batch must be between 1 and 1000")
    if args.threads < 1:
        parser.error("--threads must be >= 1")
    return args


def chunked(iterable: List[str], n: int):
    """
    Yield successive chunks of size n from a list.
    """
    buf: List[str] = []
    for x in iterable:
        if x:
            buf.append(x)
        if len(buf) == n:
            yield buf
            buf = []
    if buf:
        yield buf


def fetch_batch(
    batch: List[str],
    prefer: str,
    scheme: str,
    retry_hint: bool = True,
) -> Tuple[Dict[str, List[Tuple[str, str, str, str]]], Optional[str]]:
    """
    Query ENA filereport for a batch of accessions and return a mapping:
      { run_accession: [(kind, url, md5, bytes), ...], ... }

    - prefer: 'fastq' | 'sra' | 'submitted' | 'any'
    - scheme: 'http' | 'ftp'
    """
    params = {
        "accession": ",".join(batch),
        "result": "read_run",
        "fields": FIELDS,
        "format": "tsv",
    }
    q = ENA_API + "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(q, headers={"User-Agent": UA})
    ctx = ssl.create_default_context()

    sleep = 0.7
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            with urllib.request.urlopen(req, timeout=TIMEOUT, context=ctx) as r:
                data = r.read().decode("utf-8", "replace").splitlines()
        except Exception as e:
            if attempt >= MAX_RETRIES:
                # On final failure: mark entire batch as missing
                return {acc: [] for acc in batch}, f"[ERROR] {type(e).__name__}: {e}"
            time.sleep(sleep)
            sleep *= RETRY_BACKOFF
            continue

        # Empty or header-only result: fallback to single-accession queries
        if not data or (len(data) == 1 and data[0].strip().startswith("run_accession")):
            if retry_hint and len(batch) > 1:
                merged: Dict[str, List[Tuple[str, str, str, str]]] = {}
                errs: List[str] = []
                for acc in batch:
                    partial, err = fetch_batch([acc], prefer, scheme, retry_hint=False)
                    merged.update(partial)
                    if err:
                        errs.append(err)
                return merged, "; ".join(set(errs)) if errs else None
            else:
                return {acc: [] for acc in batch}, None

        # Parse TSV rows
        out: Dict[str, List[Tuple[str, str, str, str]]] = {}
        reader = csv.DictReader(data, delimiter="\t")

        def split_and_norm(
            row: dict, field_name: str, md5_field: str, bytes_field: str, kind: str
        ) -> List[Tuple[str, str, str, str]]:
            """
            Expand semicolon-separated *_ftp fields and normalize to full URLs.
            """
            bucket: List[Tuple[str, str, str, str]] = []
            raw = (row.get(field_name) or "").strip()
            if not raw:
                return bucket
            md5s = (row.get(md5_field) or "").strip().split(";") if md5_field else []
            byt = (row.get(bytes_field) or "").strip().split(";") if bytes_field else []
            parts = [p.strip() for p in raw.split(";") if p.strip()]
            for i, p in enumerate(parts):
                # ENA often returns host/path without scheme; add it if missing
                if p.startswith(("ftp://", "http://", "https://")):
                    url = p
                else:
                    url = ("https://" if scheme == "https" else "ftp://") + p
                m = md5s[i] if i < len(md5s) else ""
                b = byt[i] if i < len(byt) else ""
                bucket.append((kind, url, m, b))
            return bucket

        for row in reader:
            run = (row.get("run_accession") or "").strip()
            if not run:
                continue

            buckets = {
                "fastq": split_and_norm(row, "fastq_ftp", "fastq_md5", "fastq_bytes", "fastq"),
                "sra": split_and_norm(row, "sra_ftp", "sra_md5", "sra_bytes", "sra"),
                "submitted": split_and_norm(row, "submitted_ftp", "submitted_md5", "submitted_bytes", "submitted"),
            }
            all_items = buckets["fastq"] + buckets["sra"] + buckets["submitted"]
            # "prefer" uses the requested kind when available, then falls back
            # to any available link instead of incorrectly reporting it missing.
            items = buckets[prefer] if prefer != "any" and buckets[prefer] else all_items

            out[run] = items

        return out, None

    # Should not reach here
    return {acc: [] for acc in batch}, "[ERROR] Unknown state"


class _FallbackBar:
    """
    Minimal stderr progress bar used when tqdm is unavailable.
    """
    def __init__(self, total: int, desc: str = "", unit: str = "it"):
        self.total = max(total, 1)
        self.n = 0
        self.desc = desc
        self.unit = unit
        self._render()

    def update(self, inc: int = 1):
        self.n += inc
        if self.n > self.total:
            self.n = self.total
        self._render()

    def _render(self):
        width = 40
        filled = int(width * self.n / self.total)
        bar = "#" * filled + "-" * (width - filled)
        pct = (self.n / self.total) * 100.0
        sys.stderr.write(f"\r{self.desc} [{bar}] {self.n}/{self.total} ({pct:.1f}%) {self.unit}")
        sys.stderr.flush()
        if self.n >= self.total:
            sys.stderr.write("\n")

    def close(self):
        # Ensure newline after completion
        if self.n < self.total:
            sys.stderr.write("\n")


def _make_progress(total: int, desc: str, unit: str):
    """
    Create a tqdm progress bar if available; otherwise use the fallback.
    """
    if _tqdm is not None:
        return _tqdm(total=total, desc=desc, unit=unit)
    return _FallbackBar(total=total, desc=desc, unit=unit)


def main():
    args = parse_args()

    # Read accessions from file
    runs: List[str] = []
    with fileinput.input(files=args.accession_file) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                runs.append(line.split()[0])
    runs = list(dict.fromkeys(runs))
    if not runs:
        print("[ERROR] Empty input file", file=sys.stderr)
        sys.exit(1)

    out_path = Path(args.out_tsv)
    missing_path = Path(args.missing)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    missing_path.parent.mkdir(parents=True, exist_ok=True)

    total = len(runs)
    resolved: Dict[str, List[Tuple[str, str, str, str]]] = {}

    # Build progress bar over total runs
    pbar = _make_progress(total=total, desc="Resolving", unit="runs")

    # Concurrent requests
    batches = list(chunked(runs, args.batch))
    with ThreadPoolExecutor(max_workers=args.threads) as pool:
        futures = {
            pool.submit(fetch_batch, batch, args.prefer, args.scheme): batch
            for batch in batches
        }
        for fut in as_completed(futures):
            batch = futures[fut]
            mapping, err = fut.result()
            if err:
                print(err, file=sys.stderr)
            # ENA can omit an accession from an otherwise successful batch.
            # Record it explicitly so progress and the missing report stay exact.
            for acc in batch:
                resolved[acc] = mapping.get(acc, [])
            pbar.update(len(batch))

    pbar.close()

    # Futures finish out of order; write in the user's accession order to make
    # repeated runs byte-for-byte comparable.
    with open(out_path, "w", encoding="utf-8", newline="") as ftsv:
        writer = csv.writer(ftsv, delimiter="\t", lineterminator="\n")
        writer.writerow(("run_accession", "kind", "url", "md5", "bytes"))
        for acc in runs:
            writer.writerows(
                (acc, kind, url, md5, byt)
                for kind, url, md5, byt in resolved.get(acc, [])
            )

    missing_all = [acc for acc in runs if not resolved.get(acc)]

    # Write missing accessions, if any
    if missing_all:
        with open(missing_path, "w", encoding="utf-8") as f:
            for m in sorted(set(missing_all)):
                f.write(m + "\n")
        print(f"[WARN] {len(missing_all)} accessions missing; see {args.missing}", file=sys.stderr)
    else:
        missing_path.write_text("", encoding="utf-8")

    print(f"Done: {args.out_tsv} (details), {args.missing} (missing)")


if __name__ == "__main__":
    main()
