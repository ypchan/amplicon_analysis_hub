#!/usr/bin/env bash

# date: 2025-08-01
# contact:yanpengch@qq.com

# ──────────────── Usage Function ────────────────
usage() {
  cat <<EOF
Usage: amplicon_reads_lost_check.sh -i <input_file> [-o <output_file>]

Analyze DADA2 track table and calculate:
  - Non-chimera read percentage per sample
  - Merged read percentage per sample
  - Count how many samples have ≥50% nonchim or merged reads
  - If >50% samples have both <50%, suggest switching to SE analysis

Options:
  -i, --input      Input summary table (TSV format with header)
  -o, --output     Output TSV file [default: reads_lost_ratio.tsv]
  -h, --help       Show this help message

Example:
  amplicon_reads_lost_check.sh -i track.summary.tsv -o reads_lost_ratio.details.tsv
EOF
  exit 1
}

# ──────────────── Parse Arguments ────────────────
input_file=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input) input_file="$2"; shift 2 ;;
    -o|--output) output_file="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "❌ Unknown option: $1" >&2; usage ;;
  esac
done

# ──────────────── Validate Input ────────────────
if [[ -z "$input_file" || ! -f "$input_file" ]]; then
  echo "❌ Input file missing or does not exist." >&2
  usage
fi

if [[ -z "$output_file" ]];then
	outdir="$(dirname $input_file)"
	output_file="${outdir}/reads_lost_ratio.tsv"
else
	outdir="$(dirname $output_file)"
fi

# ──────────────── File Setup────────────────
echo -e "sample\tinput\tmerged\tnonchim\tmerged_pct\tnonchim_pct" > "$output_file"

summary_out="$(echo $output_file | sed 's/.tsv//;s/.details//')"".summary.tsv"

# ──────────────── Initialize Counters ────────────────
total=0
nonchim_high=0
merged_high=0
both_low=0

# Analyze line by line (skip header)
tail -n +2 "$input_file" | awk -v OFS='\t' -v out="$output_file" '
{
  sample = $1
  input = $2 + 0
  merged = $(NF-1) + 0
  nonchim = $NF + 0

  merged_pct = input > 0 ? (merged / input) * 100 : 0
  nonchim_pct = input > 0 ? (nonchim / input) * 100 : 0

  printf "%s\t%d\t%d\t%d\t%.2f\t%.2f\n", sample, input, merged, nonchim, merged_pct, nonchim_pct >> out

  total++
  if (nonchim_pct >= 50) nonchim_high++
  if (merged_pct >= 50) merged_high++
  if (nonchim_pct < 50 && merged_pct < 50) both_low++
}
END {
  printf "Sample Count                : %d\n", total
  printf "  nonchim reads left ≥ 50%%  : %d\n", nonchim_high
  printf "  merged reads ≥ 50%%        : %d\n", merged_high
  printf "  reads retained < 50%%      : %d\n", both_low

  if (both_low > total / 4)
    printf "\n⚠️  Suggestion: More than 25%% of samples have low merged and nonchim rates. Switch to SE analysis may improve results.\n\n"
}' | tee "$summary_out"

# ──────────────── Suggestion Note ────────────────
if grep -q 'Suggestion' "$summary_out";then
	touch "$outdir/suggestion.pe2se.note"
else
	touch "$outdir/suggestion.is_pe.note"
fi
echo "✅  Results saved to: $output_file"
echo "📄 Summary saved to: $summary_out"
echo ""

exit 0
