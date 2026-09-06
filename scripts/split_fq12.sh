#!/usr/bin/env bash

# Split an interleaved/mixed FASTQ using a delimited field in the read ID.

set -Eeuo pipefail

INPUT=""
R1_OUT="reads_R1.fastq.gz"
R2_OUT="reads_R2.fastq.gz"
DELIMITER='[.]'
FIELD=3
R1_VALUE=1
R2_VALUE=2
FORCE=false

usage() {
  cat <<'EOF'
Usage: split_fq12.sh -i FASTQ [options]

Options:
  -i, --input FILE       Mixed/interleaved FASTQ; .gz is detected by suffix (required)
  -1, --r1-out FILE      R1 output (default: reads_R1.fastq.gz)
  -2, --r2-out FILE      R2 output (default: reads_R2.fastq.gz)
      --delimiter REGEX  awk split regex for first header token (default: [.])
      --field INT        1-based split field holding the mate label (default: 3)
      --r1-value STR     Field value identifying R1 (default: 1)
      --r2-value STR     Field value identifying R2 (default: 2)
      --force            Replace existing output files
  -h, --help             Show this help

The input must contain complete four-line FASTQ records. Unknown mate labels are
counted and skipped. Output compression is selected independently by each .gz suffix.
EOF
}

parsed="$(getopt -o i:1:2:h -l input:,r1-out:,r2-out:,delimiter:,field:,r1-value:,r2-value:,force,help -- "$@")" || { usage >&2; exit 2; }
eval "set -- $parsed"
while true; do
  case "$1" in
    -i|--input) INPUT="$2"; shift 2 ;;
    -1|--r1-out) R1_OUT="$2"; shift 2 ;;
    -2|--r2-out) R2_OUT="$2"; shift 2 ;;
    --delimiter) DELIMITER="$2"; shift 2 ;;
    --field) FIELD="$2"; shift 2 ;;
    --r1-value) R1_VALUE="$2"; shift 2 ;;
    --r2-value) R2_VALUE="$2"; shift 2 ;;
    --force) FORCE=true; shift ;;
    -h|--help) usage; exit 0 ;;
    --) shift; break ;;
  esac
done
[[ $# -eq 0 ]] || { echo "Unexpected arguments: $*" >&2; exit 2; }
[[ -f "$INPUT" ]] || { echo "Input FASTQ not found: $INPUT" >&2; exit 2; }
[[ "$FIELD" =~ ^[1-9][0-9]*$ ]] || { echo "--field must be >= 1" >&2; exit 2; }
[[ "$R1_OUT" != "$R2_OUT" ]] || { echo "R1 and R2 outputs must differ" >&2; exit 2; }
[[ "$INPUT" != "$R1_OUT" && "$INPUT" != "$R2_OUT" ]] || { echo "Outputs must differ from input" >&2; exit 2; }
if [[ "$FORCE" == false && ( -e "$R1_OUT" || -e "$R2_OUT" ) ]]; then
  echo "Output exists; pass --force to replace it" >&2
  exit 2
fi

temp_dir="$(mktemp -d "${TMPDIR:-/tmp}/amplicon-split.XXXXXX")"
trap 'rm -rf -- "$temp_dir"' EXIT
temp_r1="$temp_dir/r1.fastq"
temp_r2="$temp_dir/r2.fastq"
counts="$temp_dir/counts.tsv"

if [[ "$INPUT" == *.gz ]]; then reader=(gzip -dc -- "$INPUT"); else reader=(cat -- "$INPUT"); fi
"${reader[@]}" | awk -v r1="$temp_r1" -v r2="$temp_r2" -v counts="$counts" \
  -v delim="$DELIMITER" -v field="$FIELD" -v one="$R1_VALUE" -v two="$R2_VALUE" '
  BEGIN { n1=0; n2=0; unknown=0; invalid=0 }
  {
    header=$0
    if ((getline sequence) <= 0 || (getline plus) <= 0 || (getline quality) <= 0) { invalid=1; exit }
    if (substr(header,1,1)!="@" || substr(plus,1,1)!="+" || length(sequence)!=length(quality)) { invalid=1; exit }
    split(header, tokens, /[[:space:]]+/)
    n=split(tokens[1], parts, delim)
    label=(field <= n ? parts[field] : "")
    record=header ORS sequence ORS plus ORS quality
    if (label==one) { print record >> r1; n1++ }
    else if (label==two) { print record >> r2; n2++ }
    else { unknown++ }
  }
  END { close(r1); close(r2); print "r1\t" n1 > counts; print "r2\t" n2 >> counts; print "unknown\t" unknown >> counts; print "invalid\t" invalid >> counts }
'

invalid="$(awk -F '\t' '$1=="invalid"{print $2}' "$counts")"
[[ "$invalid" == 0 ]] || { echo "Invalid or truncated FASTQ record" >&2; exit 1; }
[[ -e "$temp_r1" ]] || : > "$temp_r1"
[[ -e "$temp_r2" ]] || : > "$temp_r2"

write_output() {
  local source="$1" destination="$2" temp_output
  temp_output="$destination.tmp.$$"
  mkdir -p -- "$(dirname -- "$destination")"
  if [[ "$destination" == *.gz ]]; then gzip -c -- "$source" > "$temp_output"; else cp -- "$source" "$temp_output"; fi
  mv -f -- "$temp_output" "$destination"
}
write_output "$temp_r1" "$R1_OUT"
write_output "$temp_r2" "$R2_OUT"
cat "$counts"
printf 'Wrote: %s and %s\n' "$R1_OUT" "$R2_OUT"
