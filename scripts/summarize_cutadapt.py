#!/usr/bin/env python3

"""
summarize_cutadapt.py -- parse cutadapt log files and summarize primer use
date: 2025-7-30
contact: yanpengch@qq.com
"""

import re
import pandas as pd
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

def parse_cutadapt(log: str, mode: str) -> list:
    """Parse a single Cutadapt log file and extract primer information."""
    sample_name = Path(log).stem.replace('.cutadapt', '')

    with open(log, 'rt') as infh:
        txt = infh.read()

    if mode == 'PE':
        r1_name = r1_sequence = r1_strand = 'NA'
        r2_name = r2_sequence = r2_strand = 'NA'
        r1_count = r2_count = 0
        r1_percent = r2_percent = 0

        try:
            r_pairs_count = int(re.search(r"Total read pairs processed:\s+([\d,]+)", txt).group(1).replace(',', ''))
            r1_count, r1_percent = map(lambda x: float(x.strip('%')) if '%' in x else int(x),
                                       re.search(r"Read 1 with adapter:\s+(\d+)\s+\(([\d.]+%)\)", txt).groups())
            r2_count, r2_percent = map(lambda x: float(x.strip('%')) if '%' in x else int(x),
                                       re.search(r"Read 2 with adapter:\s+(\d+)\s+\(([\d.]+%)\)", txt).groups())
        except Exception:
            return [sample_name, r1_name, r1_sequence, r1_count, r1_percent, r1_strand,
                    r2_name, r2_sequence, r2_count, r2_percent, r2_strand]

        if r1_percent >= 0.5:
            r1_matches = re.findall(
                r"=== First read: Adapter (\S+).*?Sequence: ([A-Z]+);.*?Trimmed: ([\d,]+) times;.*?Reverse-complemented: ([\d,]+)",
                txt, re.S)
            r1_best = max(r1_matches, key=lambda x: max(int(x[2].replace(',', '')), int(x[3].replace(',', ''))), default=None)
            if r1_best:
                r1_name, r1_sequence, r1_trim, r1_rc = r1_best
                r1_trim, r1_rc = int(r1_trim.replace(',', '')), int(r1_rc.replace(',', ''))
                r1_count = max(r1_trim, r1_rc)
                r1_strand = '+' if r1_trim >= r1_rc else '-'

        if r2_percent >= 0.5:
            r2_matches = re.findall(
                r"=== Second read: Adapter (\S+).*?Sequence: ([A-Z]+);.*?Trimmed: ([\d,]+) times;.*?Reverse-complemented: ([\d,]+)",
                txt, re.S)
            r2_best = max(r2_matches, key=lambda x: max(int(x[2].replace(',', '')), int(x[3].replace(',', ''))), default=None)
            if r2_best:
                r2_name, r2_sequence, r2_trim, r2_rc = r2_best
                r2_trim, r2_rc = int(r2_trim.replace(',', '')), int(r2_rc.replace(',', ''))
                r2_count = max(r2_trim, r2_rc)
                r2_strand = '+' if r2_trim >= r2_rc else '-'

        return [
            sample_name, r1_name, r1_sequence, r1_count, r1_count / r_pairs_count * 100, r1_strand,
            r2_name, r2_sequence, r2_count, r2_count / r_pairs_count * 100, r2_strand
        ]

    else:  # SE mode
        r_name = r_sequence = r_strand = 'NA'
        r_count = 0
        r_percent = 0

        try:
            r_total_count = int(re.search(r"Total reads processed:\s+([\d,]+)", txt).group(1).replace(',', ''))
            r_count, r_percent = map(lambda x: float(x.strip('%')) if '%' in x else int(x),
                                     re.search(r"Reads with adapter:\s+(\d+)\s+\(([\d.]+%)\)", txt).groups())
        except Exception:
            return [sample_name, r_name, r_sequence, r_count, r_percent, r_strand]

        if r_percent >= 0.5:
            r_matches = re.findall(
                r"=== Adapter (\S+).*?Sequence: ([A-Z]+);.*?Trimmed: ([\d,]+) times;.*?Reverse-complemented: ([\d,]+)",
                txt, re.S)
            r_best = max(r_matches, key=lambda x: max(int(x[2].replace(',', '')), int(x[3].replace(',', ''))), default=None)
            if r_best:
                r_name, r_sequence, r_trim, r_rc = r_best
                r_trim, r_rc = int(r_trim.replace(',', '')), int(r_rc.replace(',', ''))
                r_count = max(r_trim, r_rc)
                r_strand = '+' if r_trim >= r_rc else '-'

        return [sample_name, r_name, r_sequence, r_count, r_count / r_total_count * 100, r_strand]


def batch_parse_cutadapt(directory: str, mode: str, threads: int = 4):
    """Parse all Cutadapt logs in a directory using multithreading and summarize primer usage."""
    log_files = list(Path(directory).glob("*.cutadapt.log"))
    total_samples = len(log_files)
    results = []

    # Multithreading for parsing
    with ThreadPoolExecutor(max_workers=threads) as executor:
        future_to_file = {executor.submit(parse_cutadapt, str(f), mode.upper()): f for f in log_files}
        for future in as_completed(future_to_file):
            result = future.result()
            results.append(result)

    # Construct detailed table
    if mode.upper() == 'PE':
        columns = ['sample_name', 'r1_name', 'r1_sequence', 'r1_count', 'r1_percent', 'r1_strand',
                   'r2_name', 'r2_sequence', 'r2_count', 'r2_percent', 'r2_strand']
    else:
        columns = ['sample_name', 'r_name', 'r_sequence', 'r_count', 'r_percent', 'r_strand']

    df = pd.DataFrame(results, columns=columns)
    df.to_csv('cutadapt_details.csv', sep='\t', index=False)

    # Primer usage summary
    def summarize(df, group_cols, primer_type):
        summary = df.groupby(group_cols).size().reset_index(name='n')
        summary.insert(0, 'primer_type', primer_type)
        summary['count'] = summary['n'].apply(lambda x: f"{x}/{total_samples}")
        summary['percent'] = summary['n'].apply(lambda x: f"{int(round(x / total_samples * 100))}%")
        summary = summary.drop(columns='n')
        return summary

    if mode.upper() == 'PE':
        r1_summary = summarize(df, ['r1_name', 'r1_sequence', 'r1_strand'], 'R1')
        r1_summary.columns = ['primer_type', 'name', 'sequence', 'strand', 'count', 'percent']

        r2_summary = summarize(df, ['r2_name', 'r2_sequence', 'r2_strand'], 'R2')
        r2_summary.columns = ['primer_type', 'name', 'sequence', 'strand', 'count', 'percent']

        summary_df = pd.concat([r1_summary, r2_summary], ignore_index=True)

    else:
        r_summary = summarize(df, ['r_name', 'r_sequence', 'r_strand'], 'SE')
        r_summary.columns = ['primer_type', 'name', 'sequence', 'strand', 'count', 'percent']
        summary_df = r_summary

    summary_df.to_csv('cutadapt_summary.csv', sep='\t', index=False)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-d', '--dir', type=str, required=True, help='Directory containing *.cutadapt.log files')
    parser.add_argument('-m', '--mode', type=str, choices=['PE', 'SE'], required=True, help='Sequencing mode: PE or SE')
    parser.add_argument('-t', '--threads', type=int, default=4, help='Number of threads to use (default: 4)')
    args = parser.parse_args()

    batch_parse_cutadapt(args.dir, args.mode, threads=args.threads)

    print("    ✅ Cutadapt parsing complete.")
    print("        📄 Details: cutadapt_details.csv")
    print("        📊 Summary: cutadapt_summary.csv")
