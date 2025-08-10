#!/usr/bin/env bash

# fastq_sorter.sh
# Sort FASTQ files by library layout, BioProject, and platform.


usage() {
  cat <<EOF
Usage: fq_sorter.sh -m <metadata.tsv> [-f <fastq_dir>] [-r <report.tsv>]

Options:
  -m, --metadata   Path to metadata TSV file (required)
  -f, --fq-dir     Directory containing FASTQ files (default: ./fq)
  -r, --report     Output classification report file (default: fq_sorting_report.tsv)
  -h, --help       Show this help message

Description:
  Sorts FASTQ files into PAIRED or SINGLE directories based on library layout, then by BioProject,
  then by sequencing platform (Illumina, Roche_454, Ion_Torrent). Generates a summary report.
EOF
  exit 1
}

fq_dir="fq"
report_file="fq_sorting_report.tsv"
metadata=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    -m|--metadata) metadata="$2"; shift 2;;
    -f|--fq-dir) fq_dir="$2"; shift 2;;
    -r|--report) report_file="$2"; shift 2;;
    -h|--help) usage;;
    *) echo "Unknown option: $1"; usage;;
  esac
done

[[ -z "$metadata" ]] && echo "Error: metadata file is required" && usage

mkdir -p PAIRED SINGLE

> "$report_file"
echo -e "LibLayout\tBioProject\tPlatform\tFASTQ_Count" >> "$report_file"

while IFS=$'\t' read -r _ accession _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ bioproject _ _ _ _ _ _ _ _ lib_layout _ platform _ _ _; do
    ll=$(echo "$lib_layout" | tr '[:lower:]' '[:upper:]')
    case "$ll" in
        SINGLE) layout_dir="SINGLE" ;;
        PAIRED) layout_dir="PAIRED" ;;
        *) continue ;;
    esac

    plat_lower=$(echo "$platform" | tr '[:upper:]' '[:lower:]')
    if [[ "$plat_lower" == *illumina* ]]; then
        plat_dir="Illumina"
    elif [[ "$plat_lower" == *roche* || "$plat_lower" == *454* ]]; then
        plat_dir="Roche_454"
    elif [[ "$plat_lower" == *ion* || "$plat_lower" == *torrent* ]]; then
        plat_dir="Ion_Torrent"
    else
        plat_dir="Other"
    fi

    dest_dir="$layout_dir/$bioproject/$plat_dir"
    mkdir -p "$dest_dir"

    count=0
    # Match files according to layout to avoid mixing single and paired patterns
    # - PAIRED: accession_1/2 or accession_R1/R2 (fastq/fq, optional .gz)
    # - SINGLE: accession (fastq/fq, optional .gz)
    shopt -s nullglob nocaseglob
    if [[ "$layout_dir" == "PAIRED" ]]; then
      patterns=(
        "$fq_dir/${accession}_1.fastq"        "$fq_dir/${accession}_1.fastq.gz"
        "$fq_dir/${accession}_2.fastq"        "$fq_dir/${accession}_2.fastq.gz"
      )
    else
      patterns=(
        "$fq_dir/${accession}.fastq"          "$fq_dir/${accession}.fastq.gz"
      )
    fi

    for pat in "${patterns[@]}"; do
      for fq in $pat; do
        [[ -e "$fq" ]] || continue
        mv "$fq" "$dest_dir"/
        ((count++))
      done
    done
    shopt -u nocaseglob nullglob

    echo -e "${layout_dir}\t${bioproject}\t${plat_dir}\t${count}" >> "$report_file"
done < "$metadata"

echo "fq files sorting finished. Report saved to $report_file"
