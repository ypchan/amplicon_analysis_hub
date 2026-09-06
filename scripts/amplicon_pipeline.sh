#!/usr/bin/env bash

# End-to-end 16S/ITS amplicon workflow for amplicon_analysis_hub.
# Input FASTQ files are read-only: exclusions are represented by a sample list,
# and cleanup removes only generated intermediates inside --output-dir.

set -Eeuo pipefail
shopt -s nullglob
ORIGINAL_ARGS=("$@")

VERSION="2.0.0"
# setup.sh installs this command as a symlink. Resolve that link before looking
# for bundled primers, classifiers, and helper scripts.
SCRIPT_PATH="$(realpath -- "${BASH_SOURCE[0]}")"
SCRIPT_DIR="$(cd -- "$(dirname -- "$SCRIPT_PATH")" && pwd -P)"
ROOT_DIR="$(cd -- "$SCRIPT_DIR/.." && pwd -P)"

THREADS=4
MODE="auto"
MARKER="16s"
PLATFORM="illumina"
R1_SUFFIX=""
R2_SUFFIX="_2.fastq.gz"
OUTPUT_DIR="amplicon_analysis_results"
PRIMER_FILE="auto"
PRIMER_MODE="trim"
DISCARD_UNTRIMMED=false
SCREEN="auto"
SCREEN_ACTION="report"
BLAST_DB="$ROOT_DIR/data/arc_bac_16s_blastDB/arch_bac_16s_ref_90"
CLASSIFIER="auto"
MIN_BOOT=50
FASTP="auto"
FASTP_MIN_LENGTH="auto"
FASTP_QUALIFIED_PHRED=15
FASTP_UNQUALIFIED_PERCENT=40
CUTADAPT_ERROR_RATE=0.10
CUTADAPT_OVERLAP=10
POOL="independent"
CLEANUP="intermediate"
FORCE=false
PRINT_PROFILE=false
declare -a DADA_OVERRIDES=()

usage_text() {
  cat <<'EOF'
amplicon_analysis_hub: end-to-end 16S/ITS ASV workflow

Usage:
  amplicon_pipeline.sh -i DIR [options]

Required:
  -i, --input-dir DIR          Directory of demultiplexed FASTQ files (never modified)

Profile selection:
  -M, --marker NAME           16s|its|other (default: 16s)
  -p, --platform NAME         illumina|mgi|element|aviti|iontorrent|454|
                              pacbio_ccs|nanopore (default: illumina)
  -m, --mode auto|pe|se       Read layout (default: auto)
  -1, --r1-suffix STR         R1/SE suffix (auto: _1.fastq.gz for PE,
                              .fastq.gz for SE)
  -2, --r2-suffix STR         R2 suffix (default: _2.fastq.gz)
  -t, --threads INT           Total concurrent CPU budget (default: 4)
  -o, --output-dir DIR        Run directory; may contain the input directory,
                              but cannot equal it or be located inside it
                              (default: amplicon_analysis_results)

Primer handling:
      --primer-file FILE      TSV: direction, name, sequence, region.
                              auto => data/16s_primer.tsv or data/its_primer.tsv
      --primer-mode trim|none Cutadapt primer/read-through removal (default: trim)
      --skip-cutadapt         Skip Cutadapt (alias for --primer-mode none)
      --discard-untrimmed     Keep only reads/pairs with a configured primer match
                              (default: off; avoids ecological primer-selection bias)
      --cutadapt-error NUM    Maximum primer error rate (default: 0.10)
      --cutadapt-overlap INT  Minimum primer match (default: 10)

Marker screening:
      --screen auto|yes|no    16S BLAST screen (default: auto: on for 16S when DB exists)
      --screen-action ACTION  report|exclude (default: report). exclude omits samples
                              from generated outputs but never deletes source FASTQs
      --blast-db PREFIX       16S nucleotide BLAST DB prefix

Quality control:
      --fastp auto|yes|no     auto => yes for short reads, no for PacBio/Nanopore
      --fastp-min-length INT  auto => 100 for 16S, 50 for ITS/other
      --fastp-qualified INT   Qualified-base Phred cutoff (default: 15)
      --fastp-unqualified NUM Maximum low-quality base percent (default: 40)

DADA2 profile and overrides:
      --pool MODE             independent|pseudo|true (default: independent)
      --trunc-len-f INT       Fixed R1/SE truncation (default: 0; ITS stays untruncated)
      --trunc-len-r INT       Fixed R2 truncation (default: 0)
      --trim-left INT         Leading bases to remove (Ion Torrent default: 15)
      --max-ee-f NUM          Expected-error limit for R1/SE (short-read default: 2)
      --max-ee-r NUM          Expected-error limit for R2 (short-read default: 2)
      --trunc-q INT           DADA2 quality truncation (short-read default: 2)
      --min-q INT             Per-base minimum (PacBio CCS default: 3)
      --min-len INT           Profile minimum length
      --max-len INT           Profile maximum length; 0 disables
      --learn-nbases NUM      Bases used to learn errors (default: 100000000)
      --min-overlap INT       PE merge overlap (default: 12)
      --max-mismatch INT      PE overlap mismatches (default: 0)
      --chimera METHOD        consensus|pooled|per-sample|none (default: consensus)
      --keep-filtered         Keep DADA2-filtered FASTQ files
      --print-profile         Print the resolved DADA2 profile and exit

Taxonomy and run control:
  -c, --classifier VALUE      auto|none|FASTA (default: auto). auto uses bundled
                              GTDB for 16S; ITS needs an explicit UNITE training FASTA
      --min-boot INT          Taxonomy bootstrap cutoff (default: 50)
      --cleanup MODE          intermediate|none (default: intermediate)
      --force                 Rerun steps even when their completion marker exists
  -h, --help                  Show this help
  -V, --version               Show version

Resolved profile defaults:
  Illumina/MGI/Element/AVITI  short PE/SE; maxEE 2/2, truncQ 2; minLen 100
                              for 16S or 50 for ITS; truncLen 0
  Ion Torrent                 SE; above + trimLeft 15 and homopolymer alignment
  Roche 454                   SE; homopolymer alignment; tune maxLen to chemistry
  PacBio CCS 16S             SE; maxEE 3, minQ 3, 1000..1800 bp, PacBioErrfun
  PacBio CCS ITS             SE; maxEE 5, minQ 3, 100..3000 bp, PacBioErrfun
  Nanopore                    Experimental SE; 16S 1000..1800 / ITS 100..3000 bp;
                              noqualErrfun. Validate with a mock community.

Examples:
  amplicon_pipeline.sh -i raw -M 16s -p illumina -m pe \
    -1 _R1.fastq.gz -2 _R2.fastq.gz -t 16
  amplicon_pipeline.sh -i raw -M its -p mgi -m pe --pool pseudo \
    -c unite_trainset.fa.gz
  amplicon_pipeline.sh -i ccs -M 16s -p pacbio_ccs -m se -1 .fastq.gz
EOF
}

usage() {
  local status="${1:-0}" line line_number=0
  local color_title='' color_section='' color_option='' color_reset=''
  if [[ -t 1 && -z "${NO_COLOR:-}" && "${TERM:-}" != "dumb" ]]; then
    color_title=$'\033[1;36m'
    color_section=$'\033[1;33m'
    color_option=$'\033[1;32m'
    color_reset=$'\033[0m'
  fi
  if [[ -z "$color_reset" ]]; then
    usage_text
  else
    while IFS= read -r line; do
      line_number=$((line_number + 1))
      if ((line_number == 1)); then
        printf '%s%s%s\n' "$color_title" "$line" "$color_reset"
      elif [[ "$line" =~ ^[^[:space:]].*:$ ]]; then
        printf '%s%s%s\n' "$color_section" "$line" "$color_reset"
      elif [[ "$line" =~ ^[[:space:]]+- ]]; then
        printf '%s%s%s%s\n' "$color_option" "${line:0:30}" "$color_reset" "${line:30}"
      else
        printf '%s\n' "$line"
      fi
    done < <(usage_text)
  fi
  exit "$status"
}

timestamp() { date '+%F %T'; }
log() { printf '[%s] [INFO] %s\n' "$(timestamp)" "$*"; }
warn() { printf '[%s] [WARN] %s\n' "$(timestamp)" "$*" >&2; }
die() { printf '[%s] [ERROR] %s\n' "$(timestamp)" "$*" >&2; exit 2; }
require_command() { command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"; }
positive_int() { [[ "$1" =~ ^[1-9][0-9]*$ ]] || die "$2 must be a positive integer"; }
blast_db_exists() {
  [[ -f "$1.nhr" || -f "$1.ndb" || -f "$1.00.nhr" ]]
}

parsed="$(getopt -o i:o:M:p:m:1:2:t:c:hV -l input-dir:,input_dir:,output-dir:,output_dir:,marker:,platform:,mode:,r1-suffix:,r1_suffix:,r2-suffix:,r2_suffix:,threads:,primer-file:,primer_file:,primer-mode:,skip-cutadapt,discard-untrimmed,cutadapt-error:,cutadapt-overlap:,screen:,screen-action:,blast-db:,fastp:,fastp-min-length:,fastp-qualified:,fastp-unqualified:,pool:,trunc-len-f:,trunc-len-r:,trim-left:,max-ee-f:,max-ee-r:,trunc-q:,min-q:,min-len:,max-len:,learn-nbases:,min-overlap:,max-mismatch:,chimera:,keep-filtered,classifier:,min-boot:,cleanup:,print-profile,force,help,version -- "$@")" || usage 2
eval "set -- $parsed"
while true; do
  case "$1" in
    -i|--input-dir|--input_dir) INPUT_DIR="$2"; shift 2 ;;
    -o|--output-dir|--output_dir) OUTPUT_DIR="$2"; shift 2 ;;
    -M|--marker) MARKER="$2"; shift 2 ;;
    -p|--platform) PLATFORM="$2"; shift 2 ;;
    -m|--mode) MODE="$2"; shift 2 ;;
    -1|--r1-suffix|--r1_suffix) R1_SUFFIX="$2"; shift 2 ;;
    -2|--r2-suffix|--r2_suffix) R2_SUFFIX="$2"; shift 2 ;;
    -t|--threads) THREADS="$2"; shift 2 ;;
    --primer-file|--primer_file) PRIMER_FILE="$2"; shift 2 ;;
    --primer-mode) PRIMER_MODE="$2"; shift 2 ;;
    --skip-cutadapt) PRIMER_MODE="none"; shift ;;
    --discard-untrimmed) DISCARD_UNTRIMMED=true; shift ;;
    --cutadapt-error) CUTADAPT_ERROR_RATE="$2"; shift 2 ;;
    --cutadapt-overlap) CUTADAPT_OVERLAP="$2"; shift 2 ;;
    --screen) SCREEN="$2"; shift 2 ;;
    --screen-action) SCREEN_ACTION="$2"; shift 2 ;;
    --blast-db) BLAST_DB="$2"; shift 2 ;;
    --fastp) FASTP="$2"; shift 2 ;;
    --fastp-min-length) FASTP_MIN_LENGTH="$2"; shift 2 ;;
    --fastp-qualified) FASTP_QUALIFIED_PHRED="$2"; shift 2 ;;
    --fastp-unqualified) FASTP_UNQUALIFIED_PERCENT="$2"; shift 2 ;;
    --pool) POOL="$2"; DADA_OVERRIDES+=(--pool "$2"); shift 2 ;;
    --trunc-len-f) DADA_OVERRIDES+=(--trunc_len_f "$2"); shift 2 ;;
    --trunc-len-r) DADA_OVERRIDES+=(--trunc_len_r "$2"); shift 2 ;;
    --trim-left) DADA_OVERRIDES+=(--trim_left "$2"); shift 2 ;;
    --max-ee-f) DADA_OVERRIDES+=(--max_ee_f "$2"); shift 2 ;;
    --max-ee-r) DADA_OVERRIDES+=(--max_ee_r "$2"); shift 2 ;;
    --trunc-q) DADA_OVERRIDES+=(--trunc_q "$2"); shift 2 ;;
    --min-q) DADA_OVERRIDES+=(--min_q "$2"); shift 2 ;;
    --min-len) DADA_OVERRIDES+=(--min_len "$2"); shift 2 ;;
    --max-len) DADA_OVERRIDES+=(--max_len "$2"); shift 2 ;;
    --learn-nbases) DADA_OVERRIDES+=(--learn_nbases "$2"); shift 2 ;;
    --min-overlap) DADA_OVERRIDES+=(--min_overlap "$2"); shift 2 ;;
    --max-mismatch) DADA_OVERRIDES+=(--max_mismatch "$2"); shift 2 ;;
    --chimera) DADA_OVERRIDES+=(--chimera "$2"); shift 2 ;;
    --keep-filtered) DADA_OVERRIDES+=(--keep_filtered); shift ;;
    -c|--classifier) CLASSIFIER="$2"; shift 2 ;;
    --min-boot) MIN_BOOT="$2"; shift 2 ;;
    --cleanup) CLEANUP="$2"; shift 2 ;;
    --print-profile) PRINT_PROFILE=true; shift ;;
    --force) FORCE=true; shift ;;
    -h|--help) usage 0 ;;
    -V|--version) printf 'amplicon_pipeline.sh %s\n' "$VERSION"; exit 0 ;;
    --) shift; break ;;
    *) die "Internal option parser error: $1" ;;
  esac
done
[[ $# -eq 0 ]] || die "Unexpected positional arguments: $*"
positive_int "$THREADS" "--threads"
[[ "$MIN_BOOT" =~ ^[0-9]+$ ]] && ((MIN_BOOT <= 100)) || die "--min-boot must be an integer in [0,100]"
positive_int "$CUTADAPT_OVERLAP" "--cutadapt-overlap"
[[ "$FASTP_QUALIFIED_PHRED" =~ ^[0-9]+$ ]] && ((FASTP_QUALIFIED_PHRED <= 93)) || die "--fastp-qualified must be an integer in [0,93]"
[[ "$FASTP_UNQUALIFIED_PERCENT" =~ ^[0-9]+$ ]] && ((FASTP_UNQUALIFIED_PERCENT <= 100)) || die "--fastp-unqualified must be an integer in [0,100]"
awk -v value="$CUTADAPT_ERROR_RATE" 'BEGIN{exit !(value>0 && value<=1)}' || die "--cutadapt-error must be in (0,1]"

MARKER="${MARKER,,}"
MODE="${MODE,,}"
PLATFORM="${PLATFORM,,}"
PLATFORM="${PLATFORM//-/_}"
case "$PLATFORM" in
  bgi|bgiseq|mgiseq|dnbseq) PLATFORM="mgi" ;;
  roche454|roche_454) PLATFORM="454" ;;
  ion_torrent) PLATFORM="iontorrent" ;;
  pacbio|ccs|hifi) PLATFORM="pacbio_ccs" ;;
  ont|oxford_nanopore) PLATFORM="nanopore" ;;
esac
[[ "$MARKER" =~ ^(16s|its|other)$ ]] || die "--marker must be 16s, its, or other"
[[ "$MODE" =~ ^(auto|pe|se)$ ]] || die "--mode must be auto, pe, or se"
[[ "$PLATFORM" =~ ^(illumina|mgi|element|aviti|iontorrent|454|pacbio_ccs|nanopore)$ ]] || die "Unsupported --platform: $PLATFORM"
[[ "$PRIMER_MODE" =~ ^(trim|none)$ ]] || die "--primer-mode must be trim or none"
[[ "$SCREEN" =~ ^(auto|yes|no)$ ]] || die "--screen must be auto, yes, or no"
[[ "$SCREEN_ACTION" =~ ^(report|exclude)$ ]] || die "--screen-action must be report or exclude"
[[ "$FASTP" =~ ^(auto|yes|no)$ ]] || die "--fastp must be auto, yes, or no"
[[ "$CLEANUP" =~ ^(intermediate|none)$ ]] || die "--cleanup must be intermediate or none"
[[ "$POOL" =~ ^(independent|pseudo|true)$ ]] || die "--pool must be independent, pseudo, or true"
if [[ "$PRINT_PROFILE" == true ]]; then
  profile_mode="$MODE"
  [[ "$profile_mode" == "auto" ]] && profile_mode="pe"
  exec Rscript "$SCRIPT_DIR/dada2.R" -M "$MARKER" -P "$PLATFORM" -m "$profile_mode" \
    "${DADA_OVERRIDES[@]}" --print_profile
fi

[[ -n "${INPUT_DIR:-}" ]] || die "--input-dir is required (see --help)"
[[ -d "$INPUT_DIR" ]] || die "Input directory not found: $INPUT_DIR"

INPUT_DIR="$(cd -- "$INPUT_DIR" && pwd -P)"
require_command realpath
OUTPUT_DIR="$(realpath -m -- "$OUTPUT_DIR")"
[[ "$OUTPUT_DIR" != "/" ]] || die "--output-dir cannot be the filesystem root"
case "$OUTPUT_DIR/" in "$INPUT_DIR/"*) die "--output-dir cannot equal or be inside --input-dir" ;; esac
mkdir -p -- "$OUTPUT_DIR"
OUTPUT_DIR="$(cd -- "$OUTPUT_DIR" && pwd -P)"
if [[ "$FORCE" == false && -f "$OUTPUT_DIR/amplicon_analysis_hub.finished" ]]; then
  log "Completed run found in $OUTPUT_DIR; use --force to rerun"
  exit 0
fi
STATE_DIR="$OUTPUT_DIR/.state"
mkdir -p -- "$STATE_DIR"

discover_with_suffix() {
  local suffix="$1" file
  while IFS= read -r -d '' file; do
    [[ "${file##*/}" == *"$suffix" ]] && printf '%s\0' "$file"
  done < <(find "$INPUT_DIR" -maxdepth 1 -type f -print0)
}

if [[ -z "$R1_SUFFIX" ]]; then
  if [[ "$MODE" == "pe" ]]; then
    R1_SUFFIX="_1.fastq.gz"
  elif [[ "$MODE" == "se" ]]; then
    R1_SUFFIX=".fastq.gz"
  else
    candidates=("$INPUT_DIR"/*_1.fastq.gz)
    if ((${#candidates[@]})); then R1_SUFFIX="_1.fastq.gz"; else R1_SUFFIX=".fastq.gz"; fi
  fi
fi

mapfile -d '' -t R1_FILES < <(discover_with_suffix "$R1_SUFFIX" | sort -z)
((${#R1_FILES[@]})) || die "No FASTQ files match suffix '$R1_SUFFIX' in $INPUT_DIR"
declare -a SAMPLES=()
for file in "${R1_FILES[@]}"; do
  base="${file##*/}"
  sample="${base%"$R1_SUFFIX"}"
  [[ -n "$sample" ]] || die "Suffix '$R1_SUFFIX' consumes an entire filename: $base"
  SAMPLES+=("$sample")
done

if [[ "$MODE" == "auto" ]]; then
  paired=true
  for sample in "${SAMPLES[@]}"; do
    [[ -f "$INPUT_DIR/$sample$R2_SUFFIX" ]] || { paired=false; break; }
  done
  if [[ "$paired" == true && "$R1_SUFFIX" != "$R2_SUFFIX" ]]; then MODE="pe"; else MODE="se"; fi
fi
if [[ "$MODE" == "pe" ]]; then
  for sample in "${SAMPLES[@]}"; do
    [[ -f "$INPUT_DIR/$sample$R2_SUFFIX" ]] || die "Missing R2: $INPUT_DIR/$sample$R2_SUFFIX"
  done
fi
if [[ "$MODE" == "pe" && "$PLATFORM" =~ ^(iontorrent|454|pacbio_ccs|nanopore)$ ]]; then
  die "$PLATFORM is supported only in SE mode"
fi

if [[ "$PRIMER_FILE" == "auto" ]]; then
  case "$MARKER" in
    16s) PRIMER_FILE="$ROOT_DIR/data/16s_primer.tsv" ;;
    its) PRIMER_FILE="$ROOT_DIR/data/its_primer.tsv" ;;
    other) PRIMER_MODE="none"; PRIMER_FILE="none" ;;
  esac
fi
[[ "$PRIMER_MODE" == "none" || -f "$PRIMER_FILE" ]] || die "Primer file not found: $PRIMER_FILE"

if [[ "$FASTP_MIN_LENGTH" == "auto" ]]; then
  if [[ "$MARKER" == "16s" ]]; then FASTP_MIN_LENGTH=100; else FASTP_MIN_LENGTH=50; fi
fi
positive_int "$FASTP_MIN_LENGTH" "--fastp-min-length"
if [[ "$FASTP" == "auto" ]]; then
  if [[ "$PLATFORM" =~ ^(pacbio_ccs|nanopore)$ ]]; then FASTP="no"; else FASTP="yes"; fi
fi
if [[ "$SCREEN" == "auto" ]]; then
  if [[ "$MARKER" == "16s" ]] && blast_db_exists "$BLAST_DB"; then
    SCREEN="yes"
  else
    [[ "$MARKER" != "16s" ]] || warn "Optional 16S BLAST DB is absent; auto screening is disabled"
    SCREEN="no"
  fi
fi
[[ "$SCREEN" == "yes" || "$SCREEN_ACTION" == "report" ]] || warn "--screen-action has no effect because screening is disabled"
[[ "$PRIMER_MODE" == "trim" || "$DISCARD_UNTRIMMED" == false ]] || warn "--discard-untrimmed has no effect because primer trimming is disabled"

case "$CLASSIFIER" in
  auto)
    if [[ "$MARKER" == "16s" && -f "$ROOT_DIR/data/gtdb_both_ssu_reps_r226.assignTaxonomy.fna.gz" ]]; then
      CLASSIFIER="$ROOT_DIR/data/gtdb_both_ssu_reps_r226.assignTaxonomy.fna.gz"
    elif [[ "$MARKER" == "16s" && -f "$ROOT_DIR/data/gtdb_both_ssu_reps_r226.assignTaxonomy.fna" ]]; then
      CLASSIFIER="$ROOT_DIR/data/gtdb_both_ssu_reps_r226.assignTaxonomy.fna"
    else
      CLASSIFIER="none"
    fi ;;
  none) ;;
  *) [[ -f "$CLASSIFIER" ]] || die "Classifier not found: $CLASSIFIER" ;;
esac

require_command seqkit
require_command Rscript
require_command cmp
require_command cksum
[[ -f "$SCRIPT_DIR/dada2.R" ]] || die "DADA2 engine not found: $SCRIPT_DIR/dada2.R"
[[ "$FASTP" == "no" ]] || require_command fastp
[[ "$PRIMER_MODE" == "none" ]] || require_command cutadapt
if [[ "$SCREEN" != "no" || "$PRIMER_MODE" != "none" ]]; then require_command python3; fi
[[ "$SCREEN" == "no" ]] || { require_command blastn; [[ -f "$SCRIPT_DIR/is_16s_amplicon.py" ]] || die "16S screen script not found"; }

parameter_tmp="$OUTPUT_DIR/.run_parameters.$$.tsv"
printf 'parameter\tvalue\n' > "$parameter_tmp"
sample_signature="$(printf '%s\0' "${SAMPLES[@]}" | cksum)"
override_string=""
if ((${#DADA_OVERRIDES[@]})); then
  printf -v override_string '%q ' "${DADA_OVERRIDES[@]}"
fi
for pair in "marker=$MARKER" "platform=$PLATFORM" "mode=$MODE" "threads=$THREADS" \
            "input_dir=$INPUT_DIR" "output_dir=$OUTPUT_DIR" "r1_suffix=$R1_SUFFIX" \
            "r2_suffix=$R2_SUFFIX" "primer_file=$PRIMER_FILE" "primer_mode=$PRIMER_MODE" \
            "discard_untrimmed=$DISCARD_UNTRIMMED" "cutadapt_error_rate=$CUTADAPT_ERROR_RATE" \
            "cutadapt_overlap=$CUTADAPT_OVERLAP" "screen=$SCREEN" \
            "screen_action=$SCREEN_ACTION" "blast_db=$BLAST_DB" "fastp=$FASTP" \
            "fastp_min_length=$FASTP_MIN_LENGTH" "fastp_qualified_phred=$FASTP_QUALIFIED_PHRED" \
            "fastp_unqualified_percent=$FASTP_UNQUALIFIED_PERCENT" "classifier=$CLASSIFIER" \
            "min_boot=$MIN_BOOT" "cleanup=$CLEANUP" "sample_count=${#SAMPLES[@]}" \
            "sample_name_cksum=$sample_signature" "dada2_overrides=$override_string"; do
  printf '%s\t%s\n' "${pair%%=*}" "${pair#*=}" >> "$parameter_tmp"
done
if [[ "$FORCE" == false && -f "$OUTPUT_DIR/run_parameters.tsv" ]] && \
   ! cmp -s "$parameter_tmp" "$OUTPUT_DIR/run_parameters.tsv"; then
  rm -f -- "$parameter_tmp"
  die "Parameters differ from this partial output directory; use a new --output-dir or --force"
fi
mv -f -- "$parameter_tmp" "$OUTPUT_DIR/run_parameters.tsv"
printf '%q ' "$0" "${ORIGINAL_ARGS[@]}" > "$OUTPUT_DIR/command.txt"
printf '\n' >> "$OUTPUT_DIR/command.txt"

log "amplicon_analysis_hub $VERSION: ${#SAMPLES[@]} samples; marker=$MARKER platform=$PLATFORM mode=$MODE"
log "Input is read-only: $INPUT_DIR"

step_done() { [[ "$FORCE" == false && -f "$STATE_DIR/$1.done" ]]; }
finish_step() { printf '%s\n' "$(timestamp)" > "$STATE_DIR/$1.done"; }

if [[ "$SCREEN" == "yes" ]]; then
  [[ "$MARKER" == "16s" ]] || die "BLAST screening is currently available only for --marker 16s"
  blast_db_exists "$BLAST_DB" || die "BLAST database files not found for prefix: $BLAST_DB"
  if ! step_done screen; then
    log "Step 0/4: reporting 16S reference hits"
    for sample in "${SAMPLES[@]}"; do printf '%s\n' "$INPUT_DIR/$sample$R1_SUFFIX"; done \
      | python3 "$SCRIPT_DIR/is_16s_amplicon.py" - --db "$BLAST_DB" \
      --threads 1 --concurrent "$THREADS" --nreads 1000 --format tsv \
      --output "$OUTPUT_DIR/16s_screen.tsv" --out-format tsv >/dev/null
    finish_step screen
  fi
  if [[ "$SCREEN_ACTION" == "exclude" ]]; then
    declare -A excluded=()
    while IFS=$'\t' read -r sample_id _ _ _ _ verdict; do
      [[ "$verdict" == "NO" ]] && excluded["$sample_id"]=1
    done < <(tail -n +2 "$OUTPUT_DIR/16s_screen.tsv")
    kept=()
    for sample in "${SAMPLES[@]}"; do
      [[ -n "${excluded["$sample$R1_SUFFIX"]:-}" ]] || kept+=("$sample")
    done
    SAMPLES=("${kept[@]}")
    ((${#SAMPLES[@]})) || die "All samples failed the 16S screen"
    warn "Screen exclusion retained ${#SAMPLES[@]} samples; source FASTQs were not changed"
  fi
fi

SOURCE_DIR="$INPUT_DIR"
if [[ "$FASTP" == "yes" ]]; then
  FASTP_DIR="$OUTPUT_DIR/01_fastp"
  mkdir -p -- "$FASTP_DIR"
  if ! step_done fastp; then
    log "Step 1/4: fastp QC"
    run_fastp() {
      local sample="$1"
      local -a cmd=(fastp -i "$INPUT_DIR/$sample$R1_SUFFIX" -o "$FASTP_DIR/$sample$R1_SUFFIX"
        --thread 1 --length_required "$FASTP_MIN_LENGTH" --n_base_limit 0
        --qualified_quality_phred "$FASTP_QUALIFIED_PHRED"
        --unqualified_percent_limit "$FASTP_UNQUALIFIED_PERCENT"
        --disable_adapter_trimming --html /dev/null
        --json "$FASTP_DIR/$sample.fastp.json")
      if [[ "$MODE" == "pe" ]]; then
        cmd+=(-I "$INPUT_DIR/$sample$R2_SUFFIX" -O "$FASTP_DIR/$sample$R2_SUFFIX")
      fi
      "${cmd[@]}" > "$FASTP_DIR/$sample.fastp.log" 2>&1
    }
    rc=0; running=0
    for sample in "${SAMPLES[@]}"; do
      run_fastp "$sample" &
      ((running+=1))
      if ((running >= THREADS)); then wait -n || rc=1; ((running-=1)); fi
    done
    while ((running > 0)); do wait -n || rc=1; ((running-=1)); done
    ((rc == 0)) || die "fastp failed; inspect $FASTP_DIR/*.fastp.log"
    finish_step fastp
  fi
  SOURCE_DIR="$FASTP_DIR"
fi

if [[ "$PRIMER_MODE" == "trim" ]]; then
  CUTADAPT_DIR="$OUTPUT_DIR/02_cutadapt"
  mkdir -p -- "$CUTADAPT_DIR"
  declare -a R1_ADAPTERS=() R2_ADAPTERS=() SE_ADAPTERS=()
  reverse_complement() { printf '%s' "$1" | tr 'ACGTRYKMSWBDHVNacgtrykmswbdhvn' 'TGCAYRMKSWVHDBNtgcayrmkswvhdbn' | rev; }
  while IFS=$'\t ' read -r direction name sequence _; do
    [[ -n "$direction" && "${direction:0:1}" != "#" && -n "$sequence" ]] || continue
    rc_sequence="$(reverse_complement "${sequence^^}")"
    safe_name="${name//[^A-Za-z0-9_.-]/_}"
    case "${direction,,}" in
      forward)
        R1_ADAPTERS+=(-g "${safe_name}=^${sequence^^}")
        R2_ADAPTERS+=(-A "${safe_name}_readthrough=$rc_sequence")
        SE_ADAPTERS+=(-g "${safe_name}=^${sequence^^}" -a "${safe_name}_readthrough=$rc_sequence") ;;
      reverse)
        R2_ADAPTERS+=(-G "${safe_name}=^${sequence^^}")
        R1_ADAPTERS+=(-a "${safe_name}_readthrough=$rc_sequence")
        SE_ADAPTERS+=(-g "${safe_name}=^${sequence^^}" -a "${safe_name}_readthrough=$rc_sequence") ;;
    esac
  done < "$PRIMER_FILE"
  ((${#SE_ADAPTERS[@]})) || die "No valid primers found in $PRIMER_FILE"
  if ! step_done cutadapt; then
    log "Step 2/4: cutadapt primer and read-through removal"
    run_cutadapt() {
      local sample="$1"
      local -a common=(--cores 1 --times 2 --revcomp -e "$CUTADAPT_ERROR_RATE"
        -O "$CUTADAPT_OVERLAP" --minimum-length "$FASTP_MIN_LENGTH")
      [[ "$DISCARD_UNTRIMMED" == false ]] || common+=(--discard-untrimmed)
      if [[ "$MODE" == "pe" ]]; then
        cutadapt "${common[@]}" "${R1_ADAPTERS[@]}" "${R2_ADAPTERS[@]}" \
          -o "$CUTADAPT_DIR/$sample$R1_SUFFIX" -p "$CUTADAPT_DIR/$sample$R2_SUFFIX" \
          "$SOURCE_DIR/$sample$R1_SUFFIX" "$SOURCE_DIR/$sample$R2_SUFFIX" \
          > "$CUTADAPT_DIR/$sample.cutadapt.log" 2>&1
      else
        cutadapt "${common[@]}" "${SE_ADAPTERS[@]}" \
          -o "$CUTADAPT_DIR/$sample$R1_SUFFIX" "$SOURCE_DIR/$sample$R1_SUFFIX" \
          > "$CUTADAPT_DIR/$sample.cutadapt.log" 2>&1
      fi
    }
    rc=0; running=0
    for sample in "${SAMPLES[@]}"; do
      run_cutadapt "$sample" &
      ((running+=1))
      if ((running >= THREADS)); then wait -n || rc=1; ((running-=1)); fi
    done
    while ((running > 0)); do wait -n || rc=1; ((running-=1)); done
    ((rc == 0)) || die "cutadapt failed; inspect $CUTADAPT_DIR/*.cutadapt.log"
    if ! python3 "$SCRIPT_DIR/summarize_cutadapt.py" --dir "$CUTADAPT_DIR" --mode "${MODE^^}" \
      --threads "$THREADS" --output-dir "$OUTPUT_DIR"; then
      warn "Cutadapt completed, but one or more text reports were incomplete; inspect cutadapt_details.tsv"
    fi
    finish_step cutadapt
  fi
  SOURCE_DIR="$CUTADAPT_DIR"
fi

log "Step 3/4: FASTQ statistics"
stats_inputs=()
for sample in "${SAMPLES[@]}"; do
  stats_inputs+=("$SOURCE_DIR/$sample$R1_SUFFIX")
  [[ "$MODE" == "se" ]] || stats_inputs+=("$SOURCE_DIR/$sample$R2_SUFFIX")
done
seqkit stats --all --tabular --threads "$THREADS" "${stats_inputs[@]}" > "$OUTPUT_DIR/seqkit.stats.tsv"

if ! step_done dada2; then
  log "Step 4/4: DADA2 ASV inference"
  dada_cmd=(Rscript "$SCRIPT_DIR/dada2.R" --input_dir "$SOURCE_DIR" --output_dir "$OUTPUT_DIR/03_dada2"
    --mode "$MODE" --marker "$MARKER" --platform "$PLATFORM"
    --reads1_suffix "$R1_SUFFIX" --threads "$THREADS" --min_boot "$MIN_BOOT"
    "${DADA_OVERRIDES[@]}")
  [[ "$MODE" == "se" ]] || dada_cmd+=(--reads2_suffix "$R2_SUFFIX")
  [[ "$CLASSIFIER" == "none" ]] || dada_cmd+=(--classifier "$CLASSIFIER")
  "${dada_cmd[@]}" 2>&1 | tee "$OUTPUT_DIR/dada2.log"
  [[ -s "$OUTPUT_DIR/03_dada2/seqtab.nochim.rds" ]] || die "DADA2 did not create seqtab.nochim.rds"
  finish_step dada2
fi

if [[ "$CLEANUP" == "intermediate" ]]; then
  log "Cleanup: removing generated fastp/cutadapt FASTQs (logs and source input are retained)"
  [[ "$FASTP" == "no" ]] || find "$OUTPUT_DIR/01_fastp" -maxdepth 1 -type f \( -name '*.fastq' -o -name '*.fastq.gz' -o -name '*.fq' -o -name '*.fq.gz' \) -delete
  [[ "$PRIMER_MODE" == "none" ]] || find "$OUTPUT_DIR/02_cutadapt" -maxdepth 1 -type f \( -name '*.fastq' -o -name '*.fastq.gz' -o -name '*.fq' -o -name '*.fq.gz' \) -delete
fi
printf '%s\n' "$(timestamp)" > "$OUTPUT_DIR/amplicon_analysis_hub.finished"
log "Finished. Main outputs: $OUTPUT_DIR/03_dada2"
