"""
unify_fq_suffix.py -- Unify FASTQ filenames, compress if needed, and standardize suffix naming.

date: 2025-07-09
contact: yanpengch@qq.com

Usage:
    unify_fq_suffix.py -i input_dir -1 _R1.fastq.gz -2 _R2.fastq.gz -f _R1.fq.gz -r _R2.fq.gz -t 8
    unify_fq_suffix.py -i input_dir -1 _R1.fastq.gz -2 _R2.fastq.gz -t 8
"""

import os
import sys
import glob
import time
import shutil
import argparse
import datetime
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

def compress_and_move(src_path, dest_path):
    if os.path.abspath(src_path) == os.path.abspath(dest_path):
        return "skipped"

    if os.path.exists(dest_path):
        return "skipped"

    if src_path.endswith(".gz"):
        shutil.move(src_path, dest_path)
    else:
        with open(dest_path, "wb") as out_f:
            subprocess.run(["gzip", "-c", src_path], stdout=out_f, check=True)
        os.remove(src_path)
    return "done"


def process_pair(r1_file, r2_file, r1_out, r2_out):
    status1 = compress_and_move(r1_file, r1_out)
    status2 = compress_and_move(r2_file, r2_out)
    return f"    ✔ {os.path.basename(r1_out)} / {os.path.basename(r2_out)} : {status2}"


def process_single(se_file, se_out):
    status = compress_and_move(se_file, se_out)
    return f"    ✔ {os.path.basename(se_out)} : {status}"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--indir", required=True, help="Input directory containing FASTQ files, required")
    parser.add_argument("-1", "--reads1_suffix_in", help="Suffix of input R1 or SE reads, optional")
    parser.add_argument("-2", "--reads2_suffix_in", help="Suffix of input R2 reads, optional")
    parser.add_argument("-f", "--reads1_suffix_std", default='_R1.fq.gz', help="Standard suffix for R1/SE reads (if -1, default: _R1.fq.gz)")
    parser.add_argument("-r", "--reads2_suffix_std", default="_R2.fq.gz", help="Standard suffix for R2reads (if -2, default: _R2.fq.gz)")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads (default: 4)")
    parser.add_argument("-s", "--trim_str", default="", help="Comma-separated patterns to remove from filenames")
    args = parser.parse_args()

    if not args.reads1_suffix_in and not args.reads2_suffix_in:
        parser.error("At least one of -1 or -2 must be specified.")

    start_time = time.time()

    indir = args.indir
    r1_suffix = args.reads1_suffix_in
    r2_suffix = args.reads2_suffix_in
    r1_suffix_std = args.reads1_suffix_std
    r2_suffix_std = args.reads2_suffix_std
    trim_patterns = [p for p in args.trim_str.split(",") if p]

    pattern = ""
    file_count = 0
    if r1_suffix and r2_suffix:
        pattern = os.path.join(indir, f"*{r1_suffix}")
        file_count = sum(1 for _ in glob.iglob(pattern, recursive=True)) * 2
    elif r1_suffix:
        pattern = os.path.join(indir,  f"*{r1_suffix}")
        file_count = sum(1 for _ in glob.iglob(pattern, recursive=True))
    elif r2_suffix:
        pattern = os.path.join(indir, f"*{r2_suffix}")
        file_count = sum(1 for _ in glob.iglob(pattern, recursive=True))

    threads = min(args.threads, max(1, file_count))
    tasks = []

    with ThreadPoolExecutor(max_workers=threads) as executor:
        if r1_suffix and r2_suffix:
            for r1_file in glob.iglob(os.path.join(indir, f"*{r1_suffix}"), recursive=True):
                r2_file = r1_file.replace(r1_suffix, r2_suffix)
                name_base = os.path.basename(r1_file).rsplit(r1_suffix, 1)[0]
                for pat in trim_patterns:
                    name_base = name_base.replace(pat, "")
                dirpath = os.path.dirname(r1_file)
                r1_out = os.path.join(dirpath, f"{name_base}{r1_suffix_std}")
                r2_out = os.path.join(dirpath, f"{name_base}{r2_suffix_std}")
                tasks.append(executor.submit(process_pair, r1_file, r2_file, r1_out, r2_out))

        elif r1_suffix:
            for se_file in glob.iglob(os.path.join(indir,  f"*{r1_suffix}"), recursive=True):
                name_base = os.path.basename(se_file).rsplit(r1_suffix, 1)[0]
                for pat in trim_patterns:
                    name_base = name_base.replace(pat, "")
                dirpath = os.path.dirname(se_file)
                se_out = os.path.join(dirpath, f"{name_base}{r1_suffix_std}")
                tasks.append(executor.submit(process_single, se_file, se_out))

        elif r2_suffix:
            for se_file in glob.iglob(os.path.join(indir,  f"*{r2_suffix}"), recursive=True):
                name_base = os.path.basename(se_file).rsplit(r2_suffix, 1)[0]
                for pat in trim_patterns:
                    name_base = name_base.replace(pat, "")
                dirpath = os.path.dirname(se_file)
                se_out = os.path.join(dirpath, f"{name_base}{r2_suffix_std}")
                tasks.append(executor.submit(process_single, se_file, se_out))

        for future in as_completed(tasks):
            try:
                print(future.result())
            except Exception as e:
                print(f"❌ Error during task: {e}", file=sys.stderr)

    elapsed = time.time() - start_time
    print(f"🎉 Finished. Elapsed time: {datetime.timedelta(seconds=int(elapsed))}")
