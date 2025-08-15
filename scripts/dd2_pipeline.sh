#!/usr/bin/env bash
#───────────────────────────────────────────────
# 🧪 DADA2 Amplicon Pipeline (Full Auto Version)
# Author: yanpengch@qq.com
# Date: 2025-08-15 (last update)
# Description: From raw FASTQ to ASV table using fastp, cutadapt and DADA2
#───────────────────────────────────────────────

#─────────────── Default Parameters ─────────────
threads=4
mode="PE"
platform="illumina"
primer_file="/home/chenyanpeng/database/16s_primer.tsv"
partition="cn"
mem_gb=500
walltime="10-00:00:00"
slurm=false
classifier=false

#─────────────── Usage Function ────────────────
usage() {
  cat <<EOF
dd2_pipeline.sh: Process amplicon sequencing data with DADA2.

Steps:
  1. fastp filtering & QC
  2. Primer detection and trimming (cutadapt)
  3. DADA2 denoising & ASV generation

Usage: dd2_pipeline.sh [options]

Required:
  -i, --input_dir DIR         Input directory with raw FASTQ files
  -1, --r1_suffix STR         R1 FASTQ suffix (e.g. _1.fastq.gz)

Optional:
  -2, --r2_suffix STR         R2 suffix (required if mode=PE)
  -t, --threads   INT         Total threads (default: 4)
  -m, --mode      SE|PE       Mode (default: PE)
  -p, --platform  STR         illumina|454|iontorrent (default: illumina)
  --primer_file   FILE        Primer table (default: ${primer_file})
  --classifier                Enable taxonomy classification step (pass through to dada2.R)
  --slurm                     Submit DADA2 via SLURM
  --partition     NAME        SLURM partition (default: ${partition})
  --mem           INT         SLURM memory GB (default: ${mem_gb})
  --request_time  D-HH:MM:SS  SLURM walltime (default: ${walltime})
  -h, --help                  Show help
EOF
  exit 1
}

#─────────────── Logging Functions ─────────────
log()  { echo -e "$(date '+[%F %T]') \033[1;32m$*\033[0m"; }
warn() { echo -e "$(date '+[%F %T]') \033[1;33m$*\033[0m"; }
err()  { echo -e "$(date '+[%F %T]') \033[1;31m$*\033[0m" >&2; }

elapsed() {
  local s=$1; local e=$(date +%s)
  printf "Elapsed time: %02d:%02d:%02d\n" $(( (e-s)/3600 )) $(( ((e-s)%3600)/60 )) $(( (e-s)%60 ))
}

#─────────────── Parse Arguments ────────────────
ARGS=$(getopt -o i:1:2:t:m:p:h -l input_dir:,r1_suffix:,r2_suffix:,threads:,mode:,platform:,primer_file:,slurm,partition:,mem:,request_time:,classifier,help -n "dada2_pipeline.sh" -- "$@") || { err "Try --help for usage."; exit 1; }
eval set -- "$ARGS"
while true; do
  case "$1" in
    -i|--input_dir) input_dir="$2"; shift 2;;
    -1|--r1_suffix) r1_suffix="$2"; shift 2;;
    -2|--r2_suffix) r2_suffix="$2"; shift 2;;
    -t|--threads) threads="$2"; shift 2;;
    -m|--mode) mode="$2"; shift 2;;
    -p|--platform) platform="$2"; shift 2;;
    --primer_file) primer_file="$2"; shift 2;;
    --slurm) slurm=true; shift;;
    --partition) partition="$2"; shift 2;;
    --classifier) classifier=true; shift;;
    --mem) mem_gb="$2"; shift 2;;
    --request_time) walltime="$2"; shift 2;;
    -h|--help) usage;;
    --) shift; break;;
    *) err "Internal error: $1"; exit 1;;
  esac
done

#─────────────── Validate Input ────────────────
[[ -z "${input_dir:-}" ]] && err "ERROR: Missing --input_dir" && usage
[[ -z "${r1_suffix:-}" ]] && err "ERROR: Missing --r1_suffix" && usage
input_dir="${input_dir%/}"

mode=$(echo "$mode" | tr '[:upper:]' '[:lower:]')
platform=$(echo "$platform" | tr '[:upper:]' '[:lower:]')

if [[ "$mode" != "pe" && "$mode" != "se" ]]; then
  err "ERROR: --mode must be SE or PE (case-insensitive)"; usage
fi
if [[ "$mode" == "pe" && -z "${r2_suffix:-}" ]]; then
  err "ERROR: --r2_suffix is required in PE mode"; usage
fi
[[ -f "$primer_file" ]] || { err "not found $primer_file"; usage; }

#─────────────── Pipeline Banner ───────────────
cat <<'EOF'

            DADA2 Amplicon Pipeline
╭──────────────────────────────────────────────╮
│  Raw Reads  → fastp  →  cutadapt  →  dada2   │
╰──────────────────────────────────────────────╯
EOF

#─────────────── Finished Check ────────────────
if [[ -f dd2_finished.note ]]; then
  log "Finished jobs in $(pwd). Nothing to do."
  exit 0
fi

#─────────────── Step 1: FASTQ Statistics ──────
log "Step 1: FASTQ statistics using seqkit"
start_t=$(date +%s)
fqfiles=$(find "$input_dir" -type f \( -name "*$r1_suffix" -o -name "*$r2_suffix" \))
[[ -z "$fqfiles" ]] && { err "No matching files found."; exit 1; }

if [[ -f seqkit.stat.tsv ]]; then
  existing_count=$(($(wc -l < seqkit.stat.tsv) - 1))
  new_count=$(echo "$fqfiles" | wc -l)
  if [[ "$existing_count" -eq "$new_count" ]]; then
    log "seqkit.stat.tsv exists and file count matches. Skipping seqkit stats."
  else
    warn "File count changed. Re-running seqkit stats..."
    seqkit stats -j "$threads" $fqfiles | sed "s|$input_dir/||;s|$r1_suffix||;s|$r2_suffix||" > seqkit.stat.tsv
  fi
else
  seqkit stats -j "$threads" $fqfiles | sed "s|$input_dir/||;s|$r1_suffix||;s|$r2_suffix||" > seqkit.stat.tsv
fi
elapsed $start_t
echo ""

#─────────────── Step 2: fastp QC ──────────────
log "Step 2: QC using fastp"
start_t=$(date +%s)
mkdir -p 01_fastp
sample_list=$(find "$input_dir" -maxdepth 2 -name "*$r1_suffix" -exec basename {} \; | sed "s/$r1_suffix//")

if [[ "$mode" == "pe" ]]; then
  echo "$sample_list" | rush -j "$threads" -v r1="$r1_suffix",r2="$r2_suffix",input_dir="$input_dir" \
    --continue --eta --succ-cmd-file fastp.rush.finished \
    'fastp -i {input_dir}/{1}{r1} -I {input_dir}/{1}{r2} -o 01_fastp/{1}{r1} -O 01_fastp/{1}{r2} --thread 1 --length_required 100 --n_base_limit 0 --cut_tail --qualified_quality_phred 20 --unqualified_percent_limit 20 --html /dev/null --json /dev/null &> 01_fastp/{1}.fastp.log'
else
  echo "$sample_list" | rush -j "$threads" -v r1="$r1_suffix",input_dir="$input_dir" \
    --continue --eta --succ-cmd-file fastp.rush.finished \
    'fastp -i {input_dir}/{1}{r1} -o 01_fastp/{1}{r1} --thread 1 --length_required 100 --n_base_limit 0 --cut_tail --qualified_quality_phred 20 --unqualified_percent_limit 20 --html /dev/null --json /dev/null &> 01_fastp/{1}.fastp.log'
fi

sample_count=$(printf '%s\n' $sample_list | sed '/^$/d' | wc -l)
fastp_finished_count=$(wc -l < fastp.rush.finished)

if (( sample_count != fastp_finished_count )) ; then
  warn "Sample count: $sample_count, fastp finished $fastp_finished_count"
  find 01_fastp -maxdepth 2 -name "*$r1_suffix" -exec basename {} \; | sed "s/$r1_suffix//" | grep -w -v -f - <(echo $sample_list | tr ' ' '\n') | sort -u > fastp.rush.failed.list
  err "ERROR: fastp failed"
  exit 1
fi

# fastp summary
if [[ "$mode" == "pe" ]]; then
  awk_cmd='
    /Read1 before filtering:/ { getline; r1in=$3 }
    /Read2 before filtering:/ { getline; r2in=$3 }
    /Read1 after filtering:/ { getline; r1out=$3 }
    /Read2 after filtering:/ { getline; r2out=$3 }
    /Insert size peak/ { insert=$NF }
    END { print sample, r1in, r2in, r1out, r2out, insert }'
else
  awk_cmd='
    /Read1 before filtering:/ { getline; r1in=$3 }
    /Read1 after filtering:/ { getline; r1out=$3 }
    END { print sample, r1in, r1out }'
fi

> fastp.filter.tsv
for f in 01_fastp/*.fastp.log; do
  sample=$(basename "$f" .fastp.log)
  awk -v sample="$sample" "$awk_cmd" "$f" >> fastp.filter.tsv
done
log "fastp summary -> fastp.filter.tsv"
log "--------------------- fastp finished. $(elapsed $start_t)"
echo ""

#─────────────── Step 3: cutadapt ──────────────
log "Step 3: cutadapt primer trimming"
start_t=$(date +%s)
mkdir -p 02_cutadapt

f_primers=$(grep '^forward' "$primer_file" | while read a b c d; do echo "-g ${b}=^${c}"; done | xargs)
r_primers=$(grep '^reverse' "$primer_file" | while read a b c d; do echo "-G ${b}=^${c}"; done | xargs)
fr_primers=$(grep -e '^forward' -e '^reverse' "$primer_file" | while read a b c d; do echo "-g ${b}=^${c}"; done | xargs)

if [[ "$mode" == "pe" ]]; then
  cutadapt_opts="$f_primers $r_primers --revcomp -j 1"
  echo "$sample_list" | rush -j "$threads" -v r1="$r1_suffix",r2="$r2_suffix",opt="${cutadapt_opts}" \
    --continue --eta --succ-cmd-file cutadapt.rush.finished \
    'cutadapt {opt} -o 02_cutadapt/{1}{r1} -p 02_cutadapt/{1}{r2} 01_fastp/{1}{r1} 01_fastp/{1}{r2} &> 02_cutadapt/{1}.cutadapt.log'
else
  cutadapt_opts="$fr_primers --revcomp -j 1"
  echo "$sample_list" | rush -j "$threads" -v r1="$r1_suffix",opt="${cutadapt_opts}" \
    --continue --eta --succ-cmd-file cutadapt.rush.finished \
    'cutadapt {opt} -o 02_cutadapt/{1}{r1} 01_fastp/{1}{r1} &> 02_cutadapt/{1}.cutadapt.log'
fi

cutadapt_finished_count=$(wc -l < cutadapt.rush.finished)
if (( sample_count != cutadapt_finished_count )) ; then
  find 02_cutadapt -maxdepth 2 -name "*$r1_suffix" -exec basename {} \; | sed "s/$r1_suffix//" | grep -w -v -f - <(echo $sample_list | tr ' ' '\n') | sort -u > cutadapt.rush.failed.list
  err "ERROR: cutadapt failed"
  exit 1
fi

if [[ $mode == "pe" ]]; then
  summarize_cutadapt.py -d 02_cutadapt/ -m PE -t $threads
else
  summarize_cutadapt.py -d 02_cutadapt/ -m SE -t $threads
fi
if [ $? -ne 0 ]; then
  err "ERROR: summarize_cutadapt.py"
  exit 1
fi
log "--------------------- cutadapt finished. $(elapsed $start_t)"
echo ""

#─────────────── Step 4: DADA2 ────────────────
log "Step 4: dada2.R"
start_t=$(date +%s)
mkdir -p 03_dada2
dd_cmd="dada2.R -i 02_cutadapt --output_dir 03_dada2 --mode $mode --reads1_suffix $r1_suffix --threads $threads --platform $platform"
[[ "$mode" == "pe" ]] && dd_cmd+=" --reads2_suffix $r2_suffix"
[[ "$classifier" == true ]] && dd_cmd+=" --classifier /mnt/nfs_ME4084storage03/chenyanpeng/database/gtdb_both_ssu_reps_r226.assignTaxonomy.fna"

cat > dada2.slurm.sh <<EOF
#!/bin/bash
#SBATCH --job-name=dd2
#SBATCH --partition=$partition
#SBATCH --output=%x.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=$threads
#SBATCH --mem=$mem_gb
#SBATCH --time=$walltime

exec 2>&1
source /home/software/miniconda3/etc/profile.d/conda.sh
conda activate dada2
$dd_cmd
EOF

if [[ "$slurm" != true ]]; then
  rm -f dada2.slurm.sh
  if ! eval "$dd_cmd" 2>&1 | tee dd2.log; then
    err "ERROR: dada2.R failed"
    exit 1
  fi
  log "$(elapsed $start_t)"
else
  log "sbatch dada2.slurm.sh"
  sbatch dada2.slurm.sh
fi

log "Step 4: check, should PE → SE?"
if [[ ! -f 03_dada2/track.summary.tsv ]]; then
  err "dada2.R error"
  exit 1
fi
amplicon_reads_lost_check.sh -i 03_dada2/track.summary.tsv

if [[ -f 03_dada2/reads_lost_ratio.summary.tsv ]]; then
  log "amplicon_reads_lost_check.sh finished"
else
  err "amplicon_reads_lost_check.sh error"
  exit 1
fi

if [[ -f 03_dada2/suggestion.pe2se.note ]]; then
  warn "PE → SE suggested, rerunning dada2.R in SE mode"
  dada2.R -i 02_cutadapt --output_dir 03_dada2 --mode SE --reads1_suffix $r1_suffix --threads $threads --platform $platform
fi

if [[ ! -f 03_dada2/seqtab.nochim.rds || ! -f 03_dada2/track.summary.tsv ]]; then
  err "ERROR: dada2 failed"
  exit 1
fi

log "--------------------- dada2 finished. $(elapsed $start_t)"
echo ""

#─────────────── Cleanup ───────────────────────
log "Step cleanup: 00_fq 01_fastp 02_cutadapt"
rm -rf 00_fq 01_fastp 02_cutadapt
log "dd2_pipeline finished."
touch dd2_finished.note
exit 0