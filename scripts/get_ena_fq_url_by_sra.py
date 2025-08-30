#!/usr/bin/env python3


import argparse
import base64
import csv
import sys
import time
import urllib.parse
import urllib.request
import ssl
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, List, Tuple, Optional

# ENA Portal API endpoint and fields to request
ENA_API = "https://www.ebi.ac.uk/ena/portal/api/filereport"
FIELDS = ",".join([
    "run_accession",
    "fastq_ftp","fastq_md5","fastq_bytes",
    "sra_ftp","sra_md5","sra_bytes",
    "submitted_ftp","submitted_md5","submitted_bytes"
])

UA = "ena-bulk-linker/2025 (non-commercial, ENA API)"
TIMEOUT = 60
MAX_RETRIES = 5
RETRY_BACKOFF = 1.8  # exponential backoff factor

def chunked(iterable, n):
    """Yield successive chunks of size n from iterable."""
    buf = []
    for x in iterable:
        if x:
            buf.append(x.strip())
        if len(buf) == n:
            yield buf
            buf = []
    if buf:
        yield buf

def thunder_encode(url: str) -> str:
    """Convert an http/ftp URL to thunder:// scheme (legacy compatibility)."""
    payload = f"AA{url}ZZ".encode("utf-8")
    return "thunder://" + base64.b64encode(payload).decode("ascii")

def fetch_batch(batch: List[str], prefer: str, scheme: str, thunder: bool, retry_hint: bool = True
               ) -> Tuple[Dict[str, List[Tuple[str,str,str,str]]], Optional[str]]:
    """
    Query ENA filereport for a batch of accessions and return a mapping:
      { run_accession: [(kind, url, md5, bytes), ...], ... }
    - prefer: 'fastq' | 'sra' | 'submitted' | 'any'
    - scheme: 'http' | 'ftp'
    - thunder: if True, convert URLs to thunder://
    """
    params = {
        "accession": ",".join(batch),
        "result": "read_run",
        "fields": FIELDS,
        "format": "tsv"
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
                # On final failure, mark all accessions in this batch as missing.
                return {acc: [] for acc in batch}, f"[ERROR] {type(e).__name__}: {e}"
            time.sleep(sleep)
            sleep *= RETRY_BACKOFF
            continue

        # If only header or empty result is returned, fallback to single-accession queries.
        if not data or (len(data) == 1 and data[0].strip().startswith("run_accession")):
            if retry_hint and len(batch) > 1:
                merged: Dict[str, List[Tuple[str,str,str,str]]] = {}
                errs: List[str] = []
                for acc in batch:
                    partial, err = fetch_batch([acc], prefer, scheme, thunder, retry_hint=False)
                    merged.update(partial)
                    if err:
                        errs.append(err)
                return merged, "; ".join(set(errs)) if errs else None
            else:
                return {acc: [] for acc in batch}, None

        # Parse TSV rows
        out: Dict[str, List[Tuple[str,str,str,str]]] = {}
        reader = csv.DictReader(data, delimiter="\t")
        for row in reader:
            run = (row.get("run_accession") or "").strip()
            if not run:
                continue

            bucket: List[Tuple[str,str,str,str]] = []

            def split_and_norm(field_name: str, md5_field: str, bytes_field: str, kind: str):
                """Expand semicolon-separated ENA *_ftp fields and normalize into full URLs."""
                raw = (row.get(field_name) or "").strip()
                if not raw:
                    return
                md5s = (row.get(md5_field) or "").strip().split(";")
                byt  = (row.get(bytes_field) or "").strip().split(";")
                parts = [p.strip() for p in raw.split(";") if p.strip()]
                for i, p in enumerate(parts):
                    # ENA typically returns host/path without scheme; add scheme if missing.
                    if p.startswith(("ftp://", "http://", "https://")):
                        url = p
                    else:
                        url = ("http://" if scheme == "http" else "ftp://") + p
                    url_out = thunder_encode(url) if thunder else url
                    m = md5s[i] if i < len(md5s) else ""
                    b = byt[i] if i < len(byt) else ""
                    bucket.append((kind, url_out, m, b))

            split_and_norm("fastq_ftp",     "fastq_md5",     "fastq_bytes",     "fastq")
            split_and_norm("sra_ftp",       "sra_md5",       "sra_bytes",       "sra")
            split_and_norm("submitted_ftp", "submitted_md5", "submitted_bytes", "submitted")

            # Sort by preference: fastq -> sra -> submitted
            priority = {"fastq": 0, "sra": 1, "submitted": 2}
            bucket.sort(key=lambda x: priority.get(x[0], 9))

            if prefer in ("fastq", "sra", "submitted"):
                bucket = [x for x in bucket if x[0] == prefer]

            out[run] = bucket

        return out, None

    # Should not reach here
    return {acc: [] for acc in batch}, "[ERROR] Unknown state"

def write_groups_by_run_tsv(
    run_rows: Dict[str, List[Tuple[str,str,str,str]]],
    run_order: List[str],
    group_size: int,
    groups_dir: Path,
    base_name: str = "ena_fq_urls"
) -> List[Path]:
    """
    
    """
    groups_dir.mkdir(parents=True, exist_ok=True)
    paths: List[Path] = []
    if group_size <= 0:
        group_size = 1000

    idx = 1
    current_count = 0
    out_fh = None
    out_path: Optional[Path] = None

    def open_new_file(nonlocal_vars):
        nonlocal idx, current_count, out_fh, out_path
        if out_fh:
            out_fh.close()
        out_path = groups_dir / f"{idx:03d}_{base_name}.tsv"
        out_fh = open(out_path, "w", encoding="utf-8")
        out_fh.write("run_accession\tkind\turl\tmd5\tbytes\n")
        paths.append(out_path)
        idx += 1
        current_count = 0

    # 初始化第一个文件
    open_new_file(locals())

    for run in run_order:
        rows = run_rows.get(run)
        if not rows:
            continue
        # 如果这个 run 放不下，换新文件（保证一个 run 不会被拆分）
        if current_count >= group_size:
            open_new_file(locals())
        for kind, url, md5, byt in rows:
            out_fh.write("\t".join([run, kind, url, md5, byt]) + "\n")
        current_count += 1

    if out_fh:
        out_fh.close()
    return paths

def main():
    parser = argparse.ArgumentParser(
        description="Resolve ENA download links for massive SRA accessions and output Xunlei group task files."
    )
    parser.add_argument("accession_file", help="Text file: one SRR/ERR/DRR per line")
    parser.add_argument("--out-urls", default="urls.txt", help="Output: all URLs in one file")
    parser.add_argument("--out-tsv",  default="links.tsv", help="Output: detailed TSV (run,kind,url,md5,bytes)")
    parser.add_argument("--missing",  default="missing.txt", help="Output: accessions with no links found")
    parser.add_argument("--prefer",   default="fastq", choices=["fastq","sra","submitted","any"], help="Preferred file type")
    parser.add_argument("--scheme",   default="http", choices=["http","ftp"], help="URL scheme for ENA links")
    parser.add_argument("--batch",    type=int, default=200, help="Accessions per API request")
    parser.add_argument("--threads",  type=int, default=8, help="Number of concurrent threads (API requests)")
    parser.add_argument("--group-size", type=int, default=1000, help="SRAs per group TSV (per-file run count)")
    parser.add_argument("--groups-dir", default="ena_fq_url", help="Directory to store grouped TSV files")  # 默认修改
    parser.add_argument("--thunder-encode", action="store_true",
                        help="Encode links as thunder:// instead of plain http/ftp")
    args = parser.parse_args()

    # Read accessions
    runs: List[str] = []
    with open(args.accession_file, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                runs.append(line)
    if not runs:
        print("[ERROR] Empty input file", file=sys.stderr)
        sys.exit(1)

    # Prepare outputs
    Path(args.out_urls).write_text("", encoding="utf-8")
    Path(args.out_tsv).write_text("run_accession\tkind\turl\tmd5\tbytes\n", encoding="utf-8")
    Path(args.missing).write_text("", encoding="utf-8")

    total = len(runs)
    processed = 0
    missing_all: List[str] = []
    all_urls: List[str] = []

    
    all_rows_by_run: Dict[str, List[Tuple[str,str,str,str]]] = {}

    # Use threads (not 'workers') for concurrency
    with ThreadPoolExecutor(max_workers=args.threads) as pool:
        futures = []
        for batch in chunked(runs, args.batch):
            futures.append(pool.submit(fetch_batch, batch, args.prefer, args.scheme, args.thunder_encode))

        for fut in as_completed(futures):
            mapping, err = fut.result()
            if err:
                print(err, file=sys.stderr)

            rows_for_outfile: List[Tuple[str,str,str,str,str]] = []
            url_lines: List[str] = []
            for acc, items in mapping.items():
                if not items:
                    missing_all.append(acc)
                    continue
                # 汇总同一个 run 的所有行，便于后续“按 run 计数分组”
                bucket = all_rows_by_run.setdefault(acc, [])
                for kind, url, md5, byt in items:
                    rows_for_outfile.append((acc, kind, url, md5, byt))
                    url_lines.append(url)
                    bucket.append((kind, url, md5, byt))

            # 追写到总 TSV 与 URL 列表（保持原有产物）
            with open(args.out_tsv, "a", encoding="utf-8") as ftsv:
                for r in rows_for_outfile:
                    ftsv.write("\t".join(r) + "\n")
            with open(args.out_urls, "a", encoding="utf-8") as furl:
                for u in url_lines:
                    furl.write(u + "\n")

            all_urls.extend(url_lines)
            processed += len(mapping)
            if processed % 10000 == 0 or processed == total:
                print(f"[{processed}/{total}] resolved", file=sys.stderr)

    # Write missing accessions
    if missing_all:
        with open(args.missing, "a", encoding="utf-8") as f:
            for m in missing_all:
                f.write(m + "\n")
        print(f"[WARN] {len(missing_all)} accessions missing; see {args.missing}", file=sys.stderr)

    # ===
    group_paths = write_groups_by_run_tsv(
        run_rows=all_rows_by_run,
        run_order=runs,              
        group_size=args.group_size,
        groups_dir=Path(args.groups_dir),
        base_name="ena_fq_urls"
    )
    print(f"Generated {len(group_paths)} grouped TSV files under '{args.groups_dir}' "
          f"(e.g., 001_ena_fq_urls.tsv).")

    print(f"Done: {args.out_urls} (all links), {args.out_tsv} (details), {args.missing} (missing)")

if __name__ == "__main__":
    main()
