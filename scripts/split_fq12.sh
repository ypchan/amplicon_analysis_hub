#!/usr/bin/env bash
# split_by_segment.sh — split mixed FASTQ into R1/R2 by header segment, with optional gzip output

set -euo pipefail

# Function to show usage message
usage() {
  cat <<EOF >&2
Usage: $(basename "$0") -i INPUT [-1 R1_OUT] [-2 R2_OUT] [-h]

  -i FILE    path to mixed FASTQ (can be .gz)
  -1 FILE    output path for R1 reads (default: reads_R1.fastq or .gz if you include .gz)
  -2 FILE    output path for R2 reads (default: reads_R2.fastq or .gz if you include .gz)
  -h         display this help and exit

If R1_OUT or R2_OUT ends with .gz, the script will gzip-compress that output.
EOF
}

# Default output filenames
r1_out="reads_R1.fastq"
r2_out="reads_R2.fastq"

# Parse command-line options
infile=""
while getopts "i:1:2:h" opt; do
  case "$opt" in
    i) infile="$OPTARG" ;;
    1) r1_out="$OPTARG" ;;
    2) r2_out="$OPTARG" ;;
    h) usage; exit 0 ;;
    *) usage; exit 1 ;;
  esac
done

# Verify input file was provided
if [[ -z "$infile" ]]; then
  echo "Error: input file is required." >&2
  usage
  exit 1
fi

# Prepare (empty) output files or set up gzip pipes
# If file ends with .gz, we'll pipe through gzip; otherwise truncate/create normally.
if [[ "$r1_out" == *.gz ]]; then
  : > /dev/null  # just ensure script doesn't error
  r1_pipe="gzip > \"$r1_out\""
  r1_mode="pipe"
else
  : > "$r1_out"
  r1_pipe="$r1_out"
  r1_mode="file"
fi

if [[ "$r2_out" == *.gz ]]; then
  : > /dev/null
  r2_pipe="gzip > \"$r2_out\""
  r2_mode="pipe"
else
  : > "$r2_out"
  r2_pipe="$r2_out"
  r2_mode="file"
fi

# Choose read command based on input extension
if [[ "$infile" == *.gz ]]; then
  read_cmd="gzip -dc \"$infile\""
else
  read_cmd="cat \"$infile\""
fi

# Process FASTQ records (4 lines each) and write to appropriate outputs
eval "$read_cmd" | awk -v R1="$r1_pipe" -v M1="$r1_mode" \
                     -v R2="$r2_pipe" -v M2="$r2_mode" '
  BEGIN {
    # nothing to do
  }
  NR % 4 == 1 {
    header = $0
    # split first field by ".", take third element as segment index
    split($1, parts, "\\.")
    segment = parts[3]
    # read the rest of the FASTQ record
    getline seq
    getline plus
    getline qual

    # choose output based on segment
    if (segment == "1") {
      if (M1 == "pipe") {
        print header ORS seq ORS plus ORS qual | R1
      } else {
        print header ORS seq ORS plus ORS qual >> R1
      }
    }
    else if (segment == "2") {
      if (M2 == "pipe") {
        print header ORS seq ORS plus ORS qual | R2
      } else {
        print header ORS seq ORS plus ORS qual >> R2
      }
    }
  }
  END {
    # close gzip pipes if used
    if (M1 == "pipe") close(R1)
    if (M2 == "pipe") close(R2)
  }
'

echo "Done splitting:"
echo "  R1 → $r1_out"
echo "  R2 → $r2_out"
