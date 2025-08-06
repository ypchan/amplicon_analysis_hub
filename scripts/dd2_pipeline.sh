#!/usr/bin/env bash

# ───────────────────────────────────────────────
# 🧪 DADA2 Amplicon Processing Pipeline  (Full Auto Version)
# Date: 2025-8-5
# Author: yanpengch@qq.com
# Description: From raw FASTQ to ASV table using fastp, cutadapt and DADA2
# ───────────────────────────────────────────────

# ──────────────── Usage Function ────────────────
usage() {
  cat <<EOF
dd2_pipeline.sh — Automated DADA2 amplicon analysis pipeline

📋 Workflow Overview:
    1️⃣  fastp         - Quality filtering and read trimming
    2️⃣  cutadapt      - Primer detection and trimming
    3️⃣  DADA2         - Denoising, ASV inference, chimera removal

🔧 Required Parameters:
    -i, --input_dir     <DIR>       Input directory containing raw FASTQ files
    -1, --r1_suffix     <STR>       Suffix for R1 reads (e.g., _R1.fq.gz)

⚙️ Optional Parameters:
    -2, --r2_suffix     <STR>       Suffix for R2 reads (required if mode=PE)
    -m, --mode          <SE|PE>     Sequencing mode: Single-End or Paired-End (default: PE)
    -p, --platform      <STR>       Sequencing platform: illumina, 454, iontorrent (default: illumina)
    -t, --threads       <INT>       Number of threads to use (default: 4)
    -r, --step          <INT>       Run specific step only:
                                      1=fastp, 2=cutadapt, 3=dada2, 0=all (default: 0)
    -s, --submit        <true|false> Submit DADA2 step via SLURM (default: false)
    -h, --help                      Show this help message and exit

📤 Example:
    dd2_pipeline.sh -i ./raw_data -1 _R1.fq.gz -2 _R2.fq.gz -m PE -t 8 -r 0

EOF
  exit 1
}

# ────────────── Parse Command-Line Arguments ──────────────
OPTIONS=i:1:2:t:m:p:r:sh
LONGOPTS=input_dir:,r1_suffix:,r2_suffix:,threads:,mode:,platform:,step:,submit,help

PARSED=$(getopt --options=$OPTIONS --longoptions=$LONGOPTS --name "$0" -- "$@") || exit 2
eval set -- "$PARSED"

# ────────────── Default Parameters ──────────────
threads=4
mode="PE"
platform="illumina"
step=0
submit=false
primer_file="/home/data/t170527/database/16s_primer.tsv"

# ────────────── Read User Parameters ──────────────
while true; do
  case "$1" in
    -i|--input_dir) input_dir="$2"; shift 2 ;;
    -1|--r1_suffix) r1_suffix="$2"; shift 2 ;;
    -2|--r2_suffix) r2_suffix="$2"; shift 2 ;;
    -t|--threads) threads="$2"; shift 2 ;;
    -m|--mode) mode="$2"; shift 2 ;;
    -p|--platform) platform="$2"; shift 2 ;;
    -r|--step) step="$2"; shift 2 ;;
    -s|--submit) submit=true; shift ;;
    -h|--help) usage ;;
    --) shift; break ;;
    *) echo "❌ Invalid option: $1" >&2; exit 3 ;;
  esac
done

# ────────────── Validate Required Inputs ──────────────
[[ -z "$input_dir" || -z "$r1_suffix" ]] && echo "❌ Required parameters missing." >&2 && usage


# ─────────────── Logging Functions ─────────────
log() { echo "$(date '+[%F %T]') $*" }

elapsed() {
  local s=$1; local e=$(date +%s)
  printf "Elapsed time: %02d:%02d:%02d\n" $(( (e-s)/3600 )) $(( ((e-s)%3600)/60 )) $(( (e-s)%60 ))
}

# ────────────── Log Configuration Summary ──────────────
log "🔧 Pipeline Configuration:"
echo "  Input directory    : $input_dir"
echo "  R1 suffix          : $r1_suffix"
echo "  R2 suffix          : $r2_suffix"
echo "  Mode               : $mode"
echo "  Threads            : $threads"
echo "  Platform           : $platform"
echo "  Step               : $step"
echo "  Submit via SLURM   : $submit"


# ────────────── Step 1: FASTQ Stats and fastp ──────────────
if [[ "$step" == 0 || "$step" == 1 ]]; then
  log "📊 Step 1: FASTQ statistics and fastp QC"
  start_t=$(date +%s)
  mkdir -p 01_fastp

  # Generate sample list and input files
  sample_list=$(find "$input_dir" -name "*$r1_suffix" -exec basename {} \; | sed "s/$r1_suffix//")
  fq_files=("$input_dir"/*"$r1_suffix")
  [[ "$mode" == "PE" ]] && fq_files+=("$input_dir"/*"$r2_suffix")

  # Run seqkit for statistics
  seqkit stats -j "$threads" "${fq_files[@]}" > 01_fastp/seqkit.stat.tsv
  log "✅ Seqkit finished. $(elapsed $start_t)"

  # Define fastp common options
  fastp_opts="--thread 1 --length_required 100 --n_base_limit 0 --cut_tail \
--qualified_quality_phred 20 --unqualified_percent_limit 20 \
--html /dev/null --json /dev/null"

  # Build rush command per mode
  if [[ "$mode" == "PE" ]]; then
    rush_vars="-v r1=\"$r1_suffix\",r2=\"$r2_suffix\",input_dir=\"$input_dir\""
    fastp_cmd="fastp -i {input_dir}/{1}{r1} -I {input_dir}/{1}{r2} \
-o 01_fastp/{1}{r1} -O 01_fastp/{1}{r2} $fastp_opts &> 01_fastp/{1}.fastp.log"
  else
    rush_vars="-v r1=\"$r1_suffix\",input_dir=\"$input_dir\""
    fastp_cmd="fastp -i {input_dir}/{1}{r1} -o 01_fastp/{1}{r1} $fastp_opts &> 01_fastp/{1}.fastp.log"
  fi

  # Execute fastp in parallel
  echo "$sample_list" | rush -j "$threads" $rush_vars \
    -c --eta --succ-cmd-file fastp.rush.done "$fastp_cmd"

  # Summarize fastp logs
  log "📊 Summarizing fastp results"
  summary="01_fastp/fastp.summary.tsv"
  [[ "$mode" == "PE" ]] && echo -e "sample\tR1_in\tR2_in\tR1_out\tR2_out\tinsert" > "$summary" || echo -e "sample\tR1_in\tR1_out" > "$summary"

  for log_file in 01_fastp/*.fastp.log; do
    sample=$(basename "$log_file" .fastp.log)
    if [[ "$mode" == "PE" ]]; then
      awk '
        /Read1 before filtering:/ {getline; r1in=$3}
        /Read2 before filtering:/ {getline; r2in=$3}
        /Read1 after filtering:/ {getline; r1out=$3}
        /Read2 after filtering:/ {getline; r2out=$3}
        /Insert size peak/ {insert=$NF}
        END {printf "%s\t%s\t%s\t%s\t%s\t%s\n", sample, r1in, r2in, r1out, r2out, insert}
      ' sample="$sample" "$log_file" >> "$summary"
    else
      awk '
        /Read1 before filtering:/ {getline; r1in=$3}
        /Read1 after filtering:/ {getline; r1out=$3}
        END {printf "%s\t%s\t%s\n", sample, r1in, r1out}
      ' sample="$sample" "$log_file" >> "$summary"
    fi
  done
  log "✅ fastp finished. $(elapsed $start_t)"
fi

# ────────────── Step 2: cutadapt ──────────────
if [[ "$step" == 0 || "$step" == 2 ]]; then
  log "✂️  Step 2: Primer trimming with cutadapt"
  start_t=$(date +%s)
  mkdir -p 02_cutadapt

  # Build primer options from primer file
  fwd_primers=$(awk '$1=="forward" {printf "-g %s=^%s ", $2, $3}' "$primer_file")
  rev_primers=$(awk '$1=="reverse" {printf "-G %s=^%s ", $2, $3}' "$primer_file")
  all_primers=$(awk '$1~/^(forward|reverse)$/ {printf "-g %s=^%s ", $2, $3}' "$primer_file")

  # Run cutadapt
  if [[ "$mode" == "PE" ]]; then
    cut_opts="$fwd_primers $rev_primers --revcomp -j 1"
    rush_cmd="cutadapt $cut_opts -o 02_cutadapt/{1}${r1_suffix} -p 02_cutadapt/{1}${r2_suffix} 01_fastp/{1}${r1_suffix} 01_fastp/{1}${r2_suffix} &> 02_cutadapt/{1}.cutadapt.log"
    rush_vars="-v r1=\"$r1_suffix\",r2=\"$r2_suffix\""
  else
    cut_opts="$all_primers --revcomp -j 1"
    rush_cmd="cutadapt $cut_opts -o 02_cutadapt/{1}${r1_suffix} 01_fastp/{1}${r1_suffix} &> 02_cutadapt/{1}.cutadapt.log"
    rush_vars="-v r1=\"$r1_suffix\""
  fi

  echo "$sample_list" | rush -j "$threads" $rush_vars -c --eta --succ-cmd-file cutadapt.rush.done "$rush_cmd"

  # Summarize primer usage
  log "📊 Summarizing primer usage"
  summarize_cutadapt.py -d 02_cutadapt -m "$mode" -t "$threads"

  log "✅ cutadapt finished. $(elapsed $start_t)"
fi

# ────────────── Step 3: DADA2 ──────────────
if [[ "$step" == 0 || "$step" == 3 ]]; then
  log "🧬 Step 3: Running DADA2 (R)"
  mkdir -p 03_dada2
  dd_cmd="dada2.R -i 02_cutadapt --output_dir 03_dada2 --mode $mode --reads1_suffix $r1_suffix --threads $threads --platform $platform"
  [[ "$mode" == "PE" ]] && dd_cmd+=" --reads2_suffix $r2_suffix"

  if [[ "$submit" == true ]]; then
    log "📤 Submitting DADA2 job via SLURM"
    cat > dada2.slurm.sh <<EOF
#!/bin/bash
#SBATCH --job-name=dada2
#SBATCH --partition=cn
#SBATCH --output=%x.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=$threads
#SBATCH --mem=500G
#SBATCH --time=8:00:00

source /home/software/miniconda3/etc/profile.d/conda.sh
conda activate dada2
$dd_cmd
EOF
    sbatch dada2.slurm.sh
    exit 0
  fi

  log "📋 Executing DADA2 locally"
  log "$dd_cmd"
  eval $dd_cmd | tee 03_dada2/dada2_run.log
  [[ $? -ne 0 ]] && log "❌ DADA2 failed. See log for details." && exit 1
  log "✅ DADA2 completed successfully."
fi

exit 0