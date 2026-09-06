#!/usr/bin/env bash

# Summarize read retention from a DADA2 track.summary.tsv by column name.

set -Eeuo pipefail
INPUT=""
OUTPUT=""
RETAINED_THRESHOLD=50
SAMPLE_FRACTION=25

usage() {
  cat <<'EOF'
Usage: amplicon_reads_lost_check.sh -i track.summary.tsv [options]

Options:
  -i, --input FILE            DADA2 track table (required)
  -o, --output FILE           Per-sample TSV (default: beside input as reads_lost_ratio.tsv)
      --retained-threshold N  Low-retention cutoff percent (default: 50)
      --sample-fraction N     Recommend review when this percent of PE samples is low
                              (default: 25)
  -h, --help                  Show this help

Columns are located by header names: input, nonchim, and optional merged. This
makes the command valid for both PE and SE track tables.
EOF
}

parsed="$(getopt -o i:o:h -l input:,output:,retained-threshold:,sample-fraction:,help -- "$@")" || { usage >&2; exit 2; }
eval "set -- $parsed"
while true; do
  case "$1" in
    -i|--input) INPUT="$2"; shift 2 ;;
    -o|--output) OUTPUT="$2"; shift 2 ;;
    --retained-threshold) RETAINED_THRESHOLD="$2"; shift 2 ;;
    --sample-fraction) SAMPLE_FRACTION="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    --) shift; break ;;
  esac
done
[[ -f "$INPUT" ]] || { echo "Input table not found: $INPUT" >&2; exit 2; }
for value in "$RETAINED_THRESHOLD" "$SAMPLE_FRACTION"; do
  [[ "$value" =~ ^([0-9]+([.][0-9]*)?|[.][0-9]+)$ ]] || { echo "Thresholds must be numeric" >&2; exit 2; }
  awk -v value="$value" 'BEGIN{exit !(value>=0 && value<=100)}' || { echo "Thresholds must be in [0,100]" >&2; exit 2; }
done

if [[ -z "$OUTPUT" ]]; then OUTPUT="$(dirname -- "$INPUT")/reads_lost_ratio.tsv"; fi
mkdir -p -- "$(dirname -- "$OUTPUT")"
summary="${OUTPUT%.tsv}.summary.tsv"
note_dir="$(dirname -- "$OUTPUT")"

awk -F '\t' -v OFS='\t' -v threshold="$RETAINED_THRESHOLD" -v fraction="$SAMPLE_FRACTION" \
  -v detail="$OUTPUT" -v summary="$summary" '
  NR==1 {
    for(i=1;i<=NF;i++){ name=$i; gsub(/^"|"$/, "", name); col[name]=i }
    if(!("input" in col) || !("nonchim" in col)){ print "Missing input/nonchim column" > "/dev/stderr"; exit 2 }
    has_merged=("merged" in col)
    print "sample","input","merged","nonchim","merged_pct","nonchim_pct" > detail
    next
  }
  {
    sample=$1; input=$(col["input"])+0; nonchim=$(col["nonchim"])+0
    if(has_merged){ merged=$(col["merged"])+0; merged_pct=input>0 ? 100*merged/input : 0 }
    else { merged="NA"; merged_pct="NA" }
    nonchim_pct=input>0 ? 100*nonchim/input : 0
    printf "%s\t%d\t%s\t%d\t%s\t%.3f\n", sample,input,merged,nonchim,(has_merged?sprintf("%.3f",merged_pct):"NA"),nonchim_pct >> detail
    total++
    if(nonchim_pct < threshold) low_nonchim++
    if(has_merged && merged_pct < threshold) low_merged++
    if(has_merged && merged_pct < threshold && nonchim_pct < threshold) low_both++
  }
  END {
    if(total==0){ print "No sample rows" > "/dev/stderr"; exit 2 }
    print "metric","value" > summary
    print "samples",total >> summary
    print "layout",(has_merged?"PE":"SE") >> summary
    print "retained_threshold_percent",threshold >> summary
    print "low_nonchim_samples",low_nonchim+0 >> summary
    print "low_merged_samples",(has_merged?low_merged+0:"NA") >> summary
    print "low_both_samples",(has_merged?low_both+0:"NA") >> summary
    low_fraction=has_merged ? 100*low_both/total : 0
    print "low_both_percent",(has_merged?sprintf("%.3f",low_fraction):"NA") >> summary
    print "review_pe_as_se",(has_merged && low_fraction>=fraction?"YES":"NO") >> summary
  }
' "$INPUT"

if awk -F '\t' '$1=="review_pe_as_se" && $2=="YES"{found=1} END{exit !found}' "$summary"; then
  rm -f -- "$note_dir/suggestion.is_pe.note"
  printf 'At least %s%% of samples have merged and nonchim retention below %s%%.\n' "$SAMPLE_FRACTION" "$RETAINED_THRESHOLD" > "$note_dir/suggestion.pe2se.note"
else
  rm -f -- "$note_dir/suggestion.pe2se.note"
  printf 'Retention review did not trigger the PE-to-SE threshold.\n' > "$note_dir/suggestion.is_pe.note"
fi
printf 'Details: %s\nSummary: %s\n' "$OUTPUT" "$summary"
