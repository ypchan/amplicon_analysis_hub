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

UA = "ena-bulk-linker/2025 (non-commercial, ENA API)"
TIMEOUT = 60
MAX_RETRIES = 5
RETRY_BACKOFF = 1.8  # Exponential backoff factor for retries


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Resolve ENA download links for massive SRA accessions."
    )
    parser.add_argument("accession_file", help="one SRR/ERR/DRR per line")
    parser.add_argument(
        "--out-tsv",
        default="links.tsv",
        help="Detailed TSV (run_accession,kind,url,md5,bytes) [default: links.tsv]",
    )
    parser.add_argument(
        "--missing",
        default="ENA_missing.tsv",
        help="Accessions with no links found [default: ENA_missing.tsv]",
    )
    parser.add_argument(
        "--prefer",
        default="fastq",
        choices=["fastq", "sra", "submitted", "any"],
        help="Preferred file type",
    )
    parser.add_argument(
        "--scheme",
        default="ftp",
        choices=["http", "ftp"],
        help="URL scheme to prepend when ENA returns host/path without scheme",
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
        help="Number of concurrent threads (API requests)",
    )
    return parser.parse_args()


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
                    url = ("http://" if scheme == "http" else "ftp://") + p
                m = md5s[i] if i < len(md5s) else ""
                b = byt[i] if i < len(byt) else ""
                bucket.append((kind, url, m, b))
            return bucket

        for row in reader:
            run = (row.get("run_accession") or "").strip()
            if not run:
                continue

            items: List[Tuple[str, str, str, str]] = []
            # Order: fastq -> sra -> submitted
            items += split_and_norm(row, "fastq_ftp", "fastq_md5", "fastq_bytes", "fastq")
            items += split_and_norm(row, "sra_ftp", "sra_md5", "sra_bytes", "sra")
            items += split_and_norm(row, "submitted_ftp", "submitted_md5", "submitted_bytes", "submitted")

            if prefer in ("fastq", "sra", "submitted"):
                items = [x for x in items if x[0] == prefer]

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
                runs.append(line)
    if not runs:
        print("[ERROR] Empty input file", file=sys.stderr)
        sys.exit(1)

    # Initialize output files
    with open(args.out_tsv, "w", encoding="utf-8") as ftsv:
        ftsv.write("run_accession\tkind\turl\tmd5\tbytes\n")
    with open(args.missing, "w", encoding="utf-8"):
        pass  # truncate

    total = len(runs)
    processed = 0
    missing_all: List[str] = []

    # Build progress bar over total runs
    pbar = _make_progress(total=total, desc="Resolving", unit="runs")

    # Concurrent requests
    with ThreadPoolExecutor(max_workers=args.threads) as pool:
        futures = []
        for batch in chunked(runs, args.batch):
            futures.append(pool.submit(fetch_batch, batch, args.prefer, args.scheme))

        for fut in as_completed(futures):
            mapping, err = fut.result()
            if err:
                print(err, file=sys.stderr)

            rows_for_outfile: List[Tuple[str, str, str, str, str]] = []
            for acc, items in mapping.items():
                if not items:
                    missing_all.append(acc)
                    continue
                for kind, url, md5, byt in items:
                    rows_for_outfile.append((acc, kind, url, md5, byt))

            # Append to TSV in the main thread to avoid file write contention
            if rows_for_outfile:
                with open(args.out_tsv, "a", encoding="utf-8") as ftsv:
                    for r in rows_for_outfile:
                        ftsv.write("\t".join(r) + "\n")

            # Update progress by the number of runs completed in this future
            processed += len(mapping)
            pbar.update(len(mapping))

    pbar.close()

    # Write missing accessions, if any
    if missing_all:
        with open(args.missing, "a", encoding="utf-8") as f:
            for m in missing_all:
                f.write(m + "\n")
        print(f"[WARN] {len(missing_all)} accessions missing; see {args.missing}", file=sys.stderr)

    print(f"Done: {args.out_tsv} (details), {args.missing} (missing)")


if __name__ == "__main__":
    main()
