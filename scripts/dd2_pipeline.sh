#!/usr/bin/env bash
#───────────────────────────────────────────────
# DADA2 Amplicon Pipeline
# Date: 2025-09-01 (last update)
# Contact: yanpengch@qq.com
# Description: From raw FASTQ to ASV table using fastp, cutadapt and DADA2
#───────────────────────────────────────────────

#─────────────── Default Parameters ─────────────
set -Eeo pipefail

THREADS=4
MODE="PE"
PLATFORM="illumina"
BLASTDB_16S="$(realpath "$(dirname -- "$(realpath "${BASH_SOURCE[0]}")")/../data/arc_bac_16s_blastDB/arch_bac_16s_ref_90")"
PRIMER_FILE="$(realpath "$(dirname -- "$(realpath "${BASH_SOURCE[0]}")")/../data/16s_primer.tsv")"
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
  printf "    Elapsed time: %02d:%02d:%02d\n" $(( (e-s)/3600 )) $(( ((e-s)%3600)/60 )) $(( (e-s)%60 ))
}

#─────────────── Parse Arguments ────────────────
ARGS=$(getopt -o i:1:2:t:m:p:h -l input_dir:,r1_suffix:,r2_suffix:,threads:,mode:,platform:,primer_file:,classifier,help -n "dada2_pipeline.sh" -- "$@") || { err "Try --help for usage."; exit 1; }
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
    --classifier) CLASSIFIER=true; shift;;
    -h|--help) usage;;
    --) shift; break;;
    *) err "Internal error: $1"; exit 1;;
  esac
done

#─────────────── Validate Input ────────────────
[[ -z "${INPUT_DIR:-}" ]] && err "Missing --input_dir" && usage
[[ -z "${R1_SUFFIX:-}" ]] && err "Missing --r1_suffix" && usage
input_dir="${INPUT_DIR%/}"

MODE=$(echo "$MODE" | tr '[:lower:]' '[:upper:]')
PLATFORM=$(echo "$PLATFORM" | tr '[:upper:]' '[:lower:]')

if [[ "$MODE" != "PE" && "$MODE" != "SE" ]]; then
  err "--mode must be SE or PE "; usage
fi

if [[ "$MODE" == "PE" && -z "${R2_SUFFIX:-}" ]]; then
  err "--r2_suffix is required in PE mode"; usage
fi

if [[ ! -f "$PRIMER_FILE" ]]; then
  err "not found $PRIMER_FILE"; usage
fi

if [[ ! -f "$BLASTDB_16S.nhr" ]]; then
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

#─────────────── STEP 0: is 16s amplicon data? ──────
log "STEP 0: Is 16S amplicon sequencing"
start_t=$(date +%s)

R1_FILE=( "$INPUT_DIR"/*"$R1_SUFFIX" )
SAMPLE_COUNT=${#R1_FILE[@]}

: > is_16s.tsv
printf '%s\n' "${R1_FILE[@]}" \
  | is_16s_amplicon.py - --db "$BLASTDB_16S" --threads 1 \
      --concurrent "$THREADS" --format tsv --output is_16s.tsv \
  >/dev/null 2>&1

[[ -s is_16s.tsv ]] || { err "is_16s.tsv not generated or empty"; exit 1; }

NON_16S_COUNT=$(awk -F '\t' '$6=="NO" {print $1}' is_16s.tsv | wc -l)
echo "------------------------------------------"
printf "|    sample  count: %-5s        |\n" "$SAMPLE_COUNT"
printf "|    non-16s count: %-5s        |\n" "$NON_16S_COUNT"
echo "-----------------------------------------|"

# Remove non-16S samples
awk -F '\t' '$6=="NO" {print $1}' is_16s.tsv | sed "s/${R1_SUFFIX}//" | while read -r a;do 
  rm -f -- "$INPUT_DIR/${a}${R1_SUFFIX}" "$INPUT_DIR/${a}${R2_SUFFIX}" 
done

if (( SAMPLE_COUNT == NON_16S_COUNT )); then
  warn "No matching files found after removing non-16S samples."
  rm -rf -- 00_fq 01_fastp 02_cutadapt
  touch dd2_finished.note
  exit 0
fi

elapsed $start_t
echo ""

#─────────────── STEP 1: FASTQ Statistics ──────
log "STEP 1: FASTQ statistics using seqkit"
start_t=$(date +%s)

if [[ -f seqkit.stat.tsv ]]; then
  SEQKIT_COUNT=$(($(wc -l < seqkit.stat.tsv) - 1))
  if [[ "$MODE" == "PE" ]]; then
    SAMPLE_COUNT=$(ls $INPUT_DIR/*$R1_SUFFIX $INPUT_DIR/*$R2_SUFFIX 2>/dev/null | wc -l)
  else
    SAMPLE_COUNT=$(ls $INPUT_DIR/*$R1_SUFFIX 2>/dev/null | wc -l)
  fi

  if (( "$SEQKIT_COUNT" == "$SAMPLE_COUNT" )) ; then
    log "seqkit.stat.tsv exists and file count matches. Skipping seqkit stats."
  else
    warn "File count changed. Re-running seqkit stats..."
    if [[ "$MODE" == "PE" ]]; then
      seqkit stats -j "$THREADS" $INPUT_DIR/*$R1_SUFFIX $INPUT_DIR/*$R2_SUFFIX | sed "s|$INPUT_DIR/||;s|$R1_SUFFIX||;s|$R2_SUFFIX||" > seqkit.stat.tsv
    else
      seqkit stats -j "$THREADS" $INPUT_DIR/*$R1_SUFFIX | sed "s|$INPUT_DIR/||;s|$R1_SUFFIX||" > seqkit.stat.tsv
    fi  
  fi
else
  if [[ "$MODE" == "PE" ]]; then
    seqkit stats -j "$THREADS" $INPUT_DIR/*$R1_SUFFIX $INPUT_DIR/*$R2_SUFFIX | sed "s|$INPUT_DIR/||;s|$R1_SUFFIX||;s|$R2_SUFFIX||" > seqkit.stat.tsv
  else
    seqkit stats -j "$THREADS" $INPUT_DIR/*$R1_SUFFIX | sed "s|$INPUT_DIR/||;s|$R1_SUFFIX||" > seqkit.stat.tsv
  fi
fi
elapsed $start_t
echo ""

#─────────────── STEP 2: fastp QC ──────────────
log "STEP 2: QC using fastp"
start_t=$(date +%s)
mkdir -p 01_fastp
SAMPLE_LIST=$(ls $INPUT_DIR/*$R1_SUFFIX | sed "s|$INPUT_DIR/||;s|$R1_SUFFIX||")

if [[ "$MODE" == "PE" ]]; then
  echo "$SAMPLE_LIST" | rush -j "$THREADS" -v r1="$R1_SUFFIX",r2="$R2_SUFFIX",input_dir="$INPUT_DIR" \
    --continue --eta --succ-cmd-file fastp.rush.finished \
    'fastp -i {input_dir}/{1}{r1} -I {input_dir}/{1}{r2} -o 01_fastp/{1}{r1} -O 01_fastp/{1}{r2} --thread 1 --length_required 100 --n_base_limit 0 --cut_tail --qualified_quality_phred 20 --unqualified_percent_limit 20 --html /dev/null --json /dev/null &> 01_fastp/{1}.fastp.log'
else
  echo "$SAMPLE_LIST" | rush -j "$THREADS" -v r1="$R1_SUFFIX",input_dir="$INPUT_DIR" \
    --continue --eta --succ-cmd-file fastp.rush.finished \
    'fastp -i {input_dir}/{1}{r1} -o 01_fastp/{1}{r1} --thread 1 --length_required 100 --n_base_limit 0 --cut_tail --qualified_quality_phred 20 --unqualified_percent_limit 20 --html /dev/null --json /dev/null &> 01_fastp/{1}.fastp.log'
fi

SAMPLE_COUNT=$(printf '%s\n' $SAMPLE_LIST | sed '/^$/d' | wc -l)
FASTP_COUNT=$(wc -l < fastp.rush.finished)

if [[ $SAMPLE_COUNT -ne $FASTP_COUNT ]]; then
  warn "Sample count: $SAMPLE_COUNT, fastp finished $FASTP_COUNT"
  ls 01_fastp/*$R1_SUFFIX | sed "s|$INPUT_DIR/||;s|$R1_SUFFIX||" | grep -w -v -f - <(echo $SAMPLE_LIST | tr ' ' '\n') | sort -u > fastp.rush.failed.list
  err "fastp failed"
  exit 1
fi

# fastp summary
if [[ "$MODE" == "PE" ]]; then
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

#─────────────── STEP 3: cutadapt ──────────────
log "STEP 3: cutadapt primer trimming"
start_t=$(date +%s)
mkdir -p 02_cutadapt

f_primers=$(grep '^forward' "$PRIMER_FILE" | while read a b c d; do echo "-g ${b}=^${c}"; done | xargs)
r_primers=$(grep '^reverse' "$PRIMER_FILE" | while read a b c d; do echo "-G ${b}=^${c}"; done | xargs)
fr_primers=$(grep -e '^forward' -e '^reverse' "$PRIMER_FILE" | while read a b c d; do echo "-g ${b}=^${c}"; done | xargs)

if [[ "$MODE" == "PE" ]]; then
  cutadapt_opts="$f_primers $r_primers --revcomp -j 1"
  echo "$SAMPLE_LIST" | rush -j "$THREADS" -v r1="$R1_SUFFIX",r2="$R2_SUFFIX",opt="${cutadapt_opts}" \
    --continue --eta --succ-cmd-file cutadapt.rush.finished \
    'cutadapt {opt} -o 02_cutadapt/{1}{r1} -p 02_cutadapt/{1}{r2} 01_fastp/{1}{r1} 01_fastp/{1}{r2} &> 02_cutadapt/{1}.cutadapt.log'
else
  cutadapt_opts="$fr_primers --revcomp -j 1"
  echo "$SAMPLE_LIST" | rush -j "$THREADS" -v r1="$R1_SUFFIX",opt="${cutadapt_opts}" \
    --continue --eta --succ-cmd-file cutadapt.rush.finished \
    'cutadapt {opt} -o 02_cutadapt/{1}{r1} 01_fastp/{1}{r1} &> 02_cutadapt/{1}.cutadapt.log'
fi

CUTADAPT_COUNT=$(wc -l < cutadapt.rush.finished)
if [[ $SAMPLE_COUNT -ne $CUTADAPT_COUNT ]]; then
  ls 02_cutadapt/*$R1_SUFFIX| sed "s|$INPUT_DIR/||;s|$R1_SUFFIX||" | grep -w -v -f - <(echo $SAMPLE_LIST | tr ' ' '\n') | sort -u > cutadapt.rush.failed.list
  err "cutadapt failed"
  exit 1
fi

if [[ $MODE == "PE" ]]; then
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

#─────────────── STEP 4: DADA2 ────────────────
log "STEP 4: dada2.R"
start_t=$(date +%s)
mkdir -p 03_dada2
dd_cmd="dada2.R -i 02_cutadapt --output_dir 03_dada2 --mode $MODE --reads1_suffix $R1_SUFFIX --threads $THREADS --platform $PLATFORM"
[[ "$MODE" == "PE" ]] && dd_cmd+=" --reads2_suffix $R2_SUFFIX"

# classifier
if [[ "$CLASSIFIER" == true ]]; then
  CLASSIFIER_REF="$(readlink -f "$(dirname -- "$(realpath "${BASH_SOURCE[0]}")")/../data/gtdb_both_ssu_reps_r226.assignTaxonomy.fna")"
  if [[ ! -f "$CLASSIFIER_REF" ]]; then
    err "Classifier reference not found: $CLASSIFIER_REF"
    exit 1
  fi
  dd_cmd+=" --classifier $CLASSIFIER_REF"
fi

if ! eval "$dd_cmd" 2>&1 | tee dd2.log; then
  err "dada2.R failed"
  exit 1
fi
log "$(elapsed $start_t)"

if [[ ! -f 03_dada2/seqtab.nochim.rds || ! -f 03_dada2/track.summary.tsv ]]; then
  err "dada2 failed"
  exit 1
fi

log "--------------------- dada2 finished. $(elapsed $start_t)"
echo ""

#─────────────── Cleanup ───────────────────────
log "STEP cleanup: 00_fq 01_fastp 02_cutadapt"
for d in 00_fq 01_fastp 02_cutadapt; do
  [[ -d "$d" && "$d" != "/" ]] && rm -rf -- "$d"
done

log "dd2_pipeline finished."
touch dd2_finished.note
exit 0
