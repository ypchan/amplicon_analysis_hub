#!/usr/bin/env bash

set -e

usage() {
cat <<EOF
check_amplicon.sh - Check if FASTQ is 16S Bacteria/Archaea amplicon using BLAST

date: 2025-08-01
contact: yanpengch@qq.com

Usage:
  $0 -i <input.fq.gz> -d <blastdb> [options]

Required:
  -i, --input     Input FASTQ (.fq.gz)

Optional:
  -d, --db        BLASTN database prefix [default: /share/cn1_fs/database/dada2_gtdb_ref/arch_bac_nr_16s]
  -n, --nreads    Number of reads to sample [default: 1000]
  -p, --identity  Identity cutoff (%) [default: 60]
  -t, --threads   Threads for BLASTN [default: 4]
  -h, --help      Show this help message
EOF
exit 1
}

# Defaults
nreads=1000
identity=60
threads=4
db="/share/cn1_fs/database/dada2_gtdb_ref/arch_bac_nr_16s"

# Parse args
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input) input="$2"; shift 2 ;;
    -d|--db) db="$2"; shift 2 ;;
    -n|--nreads) nreads="$2"; shift 2 ;;
    -p|--identity) identity="$2"; shift 2 ;;
    -t|--threads) threads="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "❌ Unknown argument: $1"; usage ;;
  esac
done

[[ -z "$input" ]] && usage
[[ ! -f "$input" ]] && echo "❌ Input file not found: $input" && exit 1

sample_name=$(basename "$input")

# Create temporary working directory
tmpdir=$(mktemp -d -t amplicon_check_XXXXXX)

# Ensure cleanup on exit
cleanup() {
  [[ -d "$tmpdir" ]] && rm -rf "$tmpdir"
}
trap cleanup EXIT

# Step 1: Sample N reads
seqkit head -j "$threads" -n "$nreads" "$input" -o "$tmpdir/sample.fq.gz"
seqkit fq2fa -j "$threads" "$tmpdir/sample.fq.gz" -o "$tmpdir/sample.fa"

# Step 2: Run BLASTN (optimized)
blastn -query "$tmpdir/sample.fa" -db "$db" \
  -out "$tmpdir/blastn.tsv" \
  -evalue 1e-5 \
  -outfmt "6 qseqid sseqid pident length qlen slen" \
  -num_threads "$threads" \
  -max_target_seqs 5 \
  -dust no \
  -word_size 20

# Step 3: Filter by identity
awk -v id="$identity" '$3 >= id' "$tmpdir/blastn.tsv" > "$tmpdir/filtered.tsv"

# Step 4: Compute alignment ratio
awk -v OFS='\t' '{printf "%s\t%s\t%.1f\t%d\t%d\t%d\t%.1f\n", $1, $2, $3, $4, $5, $6, $5/$6*100}' "$tmpdir/filtered.tsv" \
  > "$tmpdir/blastn_with_ratio.tsv"

# Step 5: Count hits
total_hits=$(cut -f1 "$tmpdir/blastn_with_ratio.tsv" | sort -u | wc -l)
bacteria_hits=$(awk '$2 ~ /^bacteria__/ {print $1}' "$tmpdir/blastn_with_ratio.tsv" | sort -u | wc -l)
archaea_hits=$(awk '$2 ~ /^archaea__/ {print $1}' "$tmpdir/blastn_with_ratio.tsv" | sort -u | wc -l)

# Step 6: Decision
percent_hits=$(awk -v t="$total_hits" -v n="$nreads" 'BEGIN {printf "%.1f", t/n*100}')
is_amplicon=$(awk -v p="$percent_hits" 'BEGIN {print (p >= 50) ? "YES" : "NO"}')

# Step 7: Output summary (saved in working dir)
summary="$tmpdir/check_summary.tsv"
mkdir -p "$(dirname "$summary")"

if [[ ! -f "$summary" ]]; then
  printf "%-20s %-10s %-10s %-12s %-15s %-12s\n" "sample_id" "bac_hits" "arch_hits" "total_hits" "total_percent" "is_16S" > "$summary"
fi
printf "%-20s %-10d %-10d %-12d %-15.1f %-12s\n" "$sample_name" "$bacteria_hits" "$archaea_hits" "$total_hits" "$percent_hits" "$is_amplicon" >> "$summary"

# Display
cat "$summary"
