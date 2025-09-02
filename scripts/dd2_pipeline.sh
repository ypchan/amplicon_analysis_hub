#!/usr/bin/env bash
#───────────────────────────────────────────────
# DADA2 Amplicon Pipeline
# Date: 2025-09-01 (last update)
# Contact: yanpengch@qq.com
# Description: From raw FASTQ to ASV table using fastp, cutadapt and DADA2
#───────────────────────────────────────────────

#─────────────── Default Parameters ─────────────
THREADS=4
MODE="PE"
PLATFORM="illumina"
BLASTDB_16S="/mnt/nfs_ME4084storage03/chenyanpeng/database/dada2_gtdb_ref/arch_bac_nr_16s"
PRIMER_FILE="/mnt/nfs_ME4084storage03/chenyanpeng/database/16s_primer.tsv"
PARTITION="cn"
MEM_GB=500G
WALLTIME="10-00:00:00"
SLURM=false
CLASSIFIER=false

#─────────────── Usage Function ────────────────
usage() {
  cat <<EOF
dd2_pipeline.sh: Process amplicon sequencing data with DADA2.

Steps:
  0. is 16S amplicon data? 
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
  --primer_file   FILE        Primer table (default: ${PRIMER_FILE})
  --classifier                Enable taxonomy classification step (pass through to dada2.R)
  --slurm                     Submit DADA2 via SLURM
  --partition     NAME        SLURM partition (default: ${PARTITION})
  --mem           INT         SLURM memory GB (default: ${MEM_GB})
  --request_time  D-HH:MM:SS  SLURM walltime (default: ${WALLTIME})
  -h, --help                  Show help

Use:
  # PE
  dd2_pipeline.sh --input_dir 00_fq --r1_suffix _1.fastq.gz --r2_suffix _2.fastq.gz --threads 24 --mode PE --platform illumina
  #SE
  dd2_pipeline.sh --input_dir 00_fq --r1_suffix .fastq.gz --threads 24 --mode SE --platform illumina
  # 454
  dd2_pipeline.sh --input_dir 00_fq --r1_suffix .fastq.gz --threads 24 --mode SE --platform 454
  # iontorrent
  dd2_pipeline.sh --input_dir 00_fq --r1_suffix .fastq.gz --threads 24 --mode SE --platform iontorrent
EOF
  exit 1
}

#─────────────── Logging Functions ─────────────
_ts() { date '+[%F %T]'; }

log()  { printf '%s %s\n' "$(_ts)" "$*"; }
warn() { printf '%s WARN: %s\n' "$(_ts)" "$*"; }
err()  { printf '%s ERROR: %s\n' "$(_ts)" "$*" >&2; }

require_cmd() { command -v "$1" >/dev/null 2>&1 || { err "Missing software: $1"; exit 127; }; }

elapsed() {
  local s=$1; local e=$(date +%s)
  printf "Elapsed time: %02d:%02d:%02d\n" $(( (e-s)/3600 )) $(( ((e-s)%3600)/60 )) $(( (e-s)%60 ))
}

#─────────────── Parse Arguments ────────────────
ARGS=$(getopt -o i:1:2:t:m:p:h -l input_dir:,r1_suffix:,r2_suffix:,threads:,mode:,platform:,primer_file:,slurm,partition:,mem:,request_time:,classifier,help -n "dada2_pipeline.sh" -- "$@") || { err "Try --help for usage."; exit 1; }
eval set -- "$ARGS"
while true; do
  case "$1" in
    -i|--input_dir) INPUT_DIR="$2"; shift 2;;
    -1|--r1_suffix) R1_SUFFIX="$2"; shift 2;;
    -2|--r2_suffix) R2_SUFFIX="$2"; shift 2;;
    -t|--threads) THREADS="$2"; shift 2;;
    -m|--mode) MODE="$2"; shift 2;;
    -p|--platform) PLATFORM="$2"; shift 2;;
    --primer_file) PRIMER_FILE="$2"; shift 2;;
    --slurm) SLURM=true; shift;;
    --partition) PARTITION="$2"; shift 2;;
    --classifier) CLASSIFIER=true; shift;;
    --mem) MEM_GB="$2"; shift 2;;
    --request_time) WALLTIME="$2"; shift 2;;
    -h|--help) usage;;
    --) shift; break;;
    *) err "Internal error: $1"; exit 1;;
  esac
done

#─────────────── Validate Input ────────────────
[[ -z "${INPUT_DIR:-}" ]] && err "Missing --input_dir" && usage
[[ -z "${R1_SUFFIX:-}" ]] && err "Missing --r1_suffix" && usage
input_dir="${INPUT_DIR%/}"

MODE=$(echo "$MODE" | tr '[:upper:]' '[:lower:]')
PLATFORM=$(echo "$PLATFORM" | tr '[:upper:]' '[:lower:]')

if [[ "$MODE" != "pe" && "$MODE" != "se" ]]; then
  err "--mode must be SE or PE "; usage
fi

if [[ "$MODE" == "pe" && -z "${R2_SUFFIX:-}" ]]; then
  err "--r2_suffix is required in PE mode"; usage
fi

if [[ ! -f "$PRIMER_FILE" ]]; then
  err "not found $PRIMER_FILE"; usage
fi

if [[ ! -f "$BLASTDB_16S.ndb" ]]; then
  err "16S DB not found: $BLASTDB_16S"
  exit 1
fi

# ─────────────── Validate softwares ────────────────
for c in fastp cutadapt seqkit rush awk sed gzip is_16s_amplicon.py summarize_cutadapt.py dada2.R; do 
  require_cmd "$c"
done

#─────────────── Pipeline Banner ───────────────
cat <<'EOF'

            DADA2 Amplicon Pipeline
╭────────────────────────────────────────────────────────╮
│  Raw Reads -> is 16s?  → fastp  →  cutadapt  →  dada2  │
╰────────────────────────────────────────────────────────╯
EOF

#─────────────── Finished Check ────────────────
if [[ -f dd2_finished.note ]]; then
  log "Finished jobs in $(pwd). Nothing to do."
  exit 0
fi

#─────────────── Input Check ────────────────
if [[ ! -d "$INPUT_DIR" ]]; then
  err "Input directory not found: $INPUT_DIR"
  exit 1
fi

#─────────────── Step 1: is 16s amplicon data? ──────
log "Step 0: Check if data is 16S amplicon sequencing"
start_t=$(date +%s)
find "$INPUT_DIR" -type f -name "*$R1_SUFFIX" \
  | is_16s_amplicon.py - --db "${BLASTDB_16S}" \
      --nreads 100 --threads 1 --concurrent "$THREADS" --format tsv \
      --output is_16s.tsv 1>/dev/null 2> is_16s.err
[[ -s is_16s.tsv ]] || { err "is_16s.tsv not generated or empty"; exit 1; }

NON_16S_COUNT=$(awk -F '\t' '$6=="NO" {print $1}' is_16s.tsv | wc -l)
echo "    non-16S sample count: $NON_16S_COUNT"

awk -F '\t' '$6=="NO" {print $1}' is_16s.tsv \
  | sed "s/${R1_SUFFIX}//" \
  | while read -r a;do \
      # Remove PE or SE reads
      rm -f "$INPUT_DIR/${a}_1.fastq.gz" \
          "$INPUT_DIR/${a}_2.fastq.gz" \
          "$INPUT_DIR/${a}.fastq.gz" \
          "$INPUT_DIR/${a}${R1_SUFFIX:-}" \
          "$INPUT_DIR/${a}${R2_SUFFIX:-}"
    done
if [[ $MODE == "pe" ]]; then
  fqfiles=$(find "$INPUT_DIR" -type f \( -name "*$R1_SUFFIX" -o -name "*$R2_SUFFIX" \))
else
  fqfiles=$(find "$INPUT_DIR" -type f -name "*$R1_SUFFIX")
fi
elapsed $start_t
[[ -z "$fqfiles" ]] && { err "No matching files found after removing non-16S samples."; exit 0; }
echo ""

#─────────────── Step 1: FASTQ Statistics ──────
log "Step 1: FASTQ statistics using seqkit"
start_t=$(date +%s)

if [[ -f seqkit.stat.tsv ]]; then
  existing_count=$(($(wc -l < seqkit.stat.tsv) - 1))
  new_count=$(echo "$fqfiles" | wc -l)
  if [[ "$existing_count" -eq "$new_count" ]]; then
    log "seqkit.stat.tsv exists and file count matches. Skipping seqkit stats."
  else
    warn "File count changed. Re-running seqkit stats..."
    seqkit stats -j "$THREADS" $fqfiles | sed "s|$INPUT_DIR/||;s|$R1_SUFFIX||;s|$R2_SUFFIX||" > seqkit.stat.tsv
  fi
else
  seqkit stats -j "$THREADS" $fqfiles | sed "s|$INPUT_DIR/||;s|$R1_SUFFIX||;s|$R2_SUFFIX||" > seqkit.stat.tsv
fi
elapsed $start_t
echo ""

#─────────────── Step 2: fastp QC ──────────────
log "Step 2: QC using fastp"
start_t=$(date +%s)
mkdir -p 01_fastp
sample_list=$(find "$INPUT_DIR" -maxdepth 2 -name "*$R1_SUFFIX" -exec basename {} \; | sed "s/$R1_SUFFIX//")

if [[ "$MODE" == "pe" ]]; then
  echo "$sample_list" | rush -j "$THREADS" -v r1="$R1_SUFFIX",r2="$R2_SUFFIX",input_dir="$INPUT_DIR" \
    --continue --eta --succ-cmd-file fastp.rush.finished \
    'fastp -i {input_dir}/{1}{r1} -I {input_dir}/{1}{r2} -o 01_fastp/{1}{r1} -O 01_fastp/{1}{r2} --thread 1 --length_required 100 --n_base_limit 0 --cut_tail --qualified_quality_phred 20 --unqualified_percent_limit 20 --html /dev/null --json /dev/null &> 01_fastp/{1}.fastp.log'
else
  echo "$sample_list" | rush -j "$THREADS" -v r1="$R1_SUFFIX",input_dir="$INPUT_DIR" \
    --continue --eta --succ-cmd-file fastp.rush.finished \
    'fastp -i {input_dir}/{1}{r1} -o 01_fastp/{1}{r1} --thread 1 --length_required 100 --n_base_limit 0 --cut_tail --qualified_quality_phred 20 --unqualified_percent_limit 20 --html /dev/null --json /dev/null &> 01_fastp/{1}.fastp.log'
fi

sample_count=$(printf '%s\n' $sample_list | sed '/^$/d' | wc -l)
fastp_finished_count=$(wc -l < fastp.rush.finished)

if (( sample_count != fastp_finished_count )) ; then
  warn "Sample count: $sample_count, fastp finished $fastp_finished_count"
  find 01_fastp -maxdepth 2 -name "*$R1_SUFFIX" -exec basename {} \; | sed "s/$R1_SUFFIX//" | grep -w -v -f - <(echo $sample_list | tr ' ' '\n') | sort -u > fastp.rush.failed.list
  err "fastp failed"
  exit 1
fi

# fastp summary
if [[ "$MODE" == "pe" ]]; then
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

f_primers=$(grep '^forward' "$PRIMER_FILE" | while read a b c d; do echo "-g ${b}=^${c}"; done | xargs)
r_primers=$(grep '^reverse' "$PRIMER_FILE" | while read a b c d; do echo "-G ${b}=^${c}"; done | xargs)
fr_primers=$(grep -e '^forward' -e '^reverse' "$PRIMER_FILE" | while read a b c d; do echo "-g ${b}=^${c}"; done | xargs)

if [[ "$MODE" == "pe" ]]; then
  cutadapt_opts="$f_primers $r_primers --revcomp -j 1"
  echo "$sample_list" | rush -j "$THREADS" -v r1="$R1_SUFFIX",r2="$R2_SUFFIX",opt="${cutadapt_opts}" \
    --continue --eta --succ-cmd-file cutadapt.rush.finished \
    'cutadapt {opt} -o 02_cutadapt/{1}{r1} -p 02_cutadapt/{1}{r2} 01_fastp/{1}{r1} 01_fastp/{1}{r2} &> 02_cutadapt/{1}.cutadapt.log'
else
  cutadapt_opts="$fr_primers --revcomp -j 1"
  echo "$sample_list" | rush -j "$THREADS" -v r1="$R1_SUFFIX",opt="${cutadapt_opts}" \
    --continue --eta --succ-cmd-file cutadapt.rush.finished \
    'cutadapt {opt} -o 02_cutadapt/{1}{r1} 01_fastp/{1}{r1} &> 02_cutadapt/{1}.cutadapt.log'
fi

cutadapt_finished_count=$(wc -l < cutadapt.rush.finished)
if (( sample_count != cutadapt_finished_count )) ; then
  find 02_cutadapt -maxdepth 2 -name "*$R1_SUFFIX" -exec basename {} \; | sed "s/$R1_SUFFIX//" | grep -w -v -f - <(echo $sample_list | tr ' ' '\n') | sort -u > cutadapt.rush.failed.list
  err "cutadapt failed"
  exit 1
fi

if [[ $MODE == "pe" ]]; then
  summarize_cutadapt.py -d 02_cutadapt/ -m PE -t $THREADS
else
  summarize_cutadapt.py -d 02_cutadapt/ -m SE -t $THREADS
fi
if [ $? -ne 0 ]; then
  err "summarize_cutadapt.py"
  exit 1
fi
log "--------------------- cutadapt finished. $(elapsed $start_t)"
echo ""

#─────────────── Step 4: DADA2 ────────────────
log "Step 4: dada2.R"
start_t=$(date +%s)
mkdir -p 03_dada2
dd_cmd="dada2.R -i 02_cutadapt --output_dir 03_dada2 --mode $MODE --reads1_suffix $R1_SUFFIX --threads $THREADS --platform $PLATFORM"
[[ "$MODE" == "pe" ]] && dd_cmd+=" --reads2_suffix $R2_SUFFIX"
[[ "$CLASSIFIER" == true ]] && dd_cmd+=" --classifier /mnt/nfs_ME4084storage03/chenyanpeng/database/gtdb_both_ssu_reps_r226.assignTaxonomy.fna"

cat > dada2.slurm.sh <<EOF
#!/bin/bash
#SBATCH --job-name=dd2
#SBATCH --partition=$PARTITION
#SBATCH --output=%x.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=$THREADS
#SBATCH --mem=$MEM_GB
#SBATCH --time=$WALLTIME

exec 2>&1
source /home/software/miniconda3/etc/profile.d/conda.sh
conda activate dada2
$dd_cmd
EOF

if [[ "$SLURM" != true ]]; then
  rm -f dada2.slurm.sh
  if ! eval "$dd_cmd" 2>&1 | tee dd2.log; then
    err "dada2.R failed"
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
if [[ $MODE == "pe" ]]; then
  amplicon_reads_lost_check.sh -i 03_dada2/track.summary.tsv
else
  amplicon_reads_lost_check.sh -i 03_dada2/track.summary.tsv &>/dev/null
fi

if [[ -f 03_dada2/reads_lost_ratio.summary.tsv ]]; then
  log "amplicon_reads_lost_check.sh finished"
else
  err "amplicon_reads_lost_check.sh error"
  exit 1
fi

if [[ -f 03_dada2/suggestion.pe2se.note && "$MODE" == "pe" ]]; then
  warn "PE → SE suggested, rerunning dada2.R in SE mode"
  dada2.R -i 02_cutadapt --output_dir 03_dada2 --mode SE --reads1_suffix $R1SUFFIX --threads $THREADS --platform $PLATFORM
fi

if [[ ! -f 03_dada2/seqtab.nochim.rds || ! -f 03_dada2/track.summary.tsv ]]; then
  err "dada2 failed"
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
