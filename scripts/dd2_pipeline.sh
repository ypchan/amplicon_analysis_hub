#!/usr/bin/env bash

# ───────────────────────────────────────────────
# 🧪 DADA2 Amplicon Pipeline (Full Auto Version)
# Author: yanpengch@qq.com
# Description: From raw FASTQ to ASV table using fastp, cutadapt and DADA2
# ───────────────────────────────────────────────

# ─────────────── Default Parameters ─────────────
threads=4
mode="PE"
platform="illumina"
primer_file="/home/chenyanpeng/database/16s_primer.tsv"
partition="cn"
mem_gb=500
walltime="10-00:00:00"
slurm=false
classifier=false


# ──────────────── Usage Function ────────────────
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
  -1, --r1_suffix STR         R1 FASTQ suffix (e.g. _R1.fastq.gz)

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

# ─────────────── Parse Arguments ────────────────
# NOTE: short options -1/-2 are supported by GNU getopt
ARGS=$(getopt -o i:1:2:t:m:p:h -l input_dir:,r1_suffix:,r2_suffix:,threads:,mode:,platform:,primer_file:,slurm,partition:,mem:,request_time:,classifier,help -n "dada2_pipeline.sh" -- "$@") || { echo "Try --help for usage." >&2; exit 1; }
eval set -- "$ARGS"
while true; do
  case "$1" in
    -i|--input_dir) input_dir="$2";   shift 2;;
    -1|--r1_suffix) r1_suffix="$2";   shift 2;;
    -2|--r2_suffix) r2_suffix="$2";   shift 2;;
    -t|--threads)   threads="$2";     shift 2;;
    -m|--mode)      mode="$2";        shift 2;;
    -p|--platform)  platform="$2";    shift 2;;
    --primer_file)  primer_file="$2"; shift 2;;
    --slurm)        slurm=true;       shift;;
    --partition)    partition="$2";   shift 2;;
    --classifier)   classifier=true;  shift;;
    --mem)          mem_gb="$2";      shift 2;;
    --request_time) walltime="$2";    shift 2;;
    -h|--help)      usage;            exit 0;;
    --)             shift;            break;;
    *) echo "Internal error: $1" >&2; exit 1;;
  esac
done

# ─────────────── Validate Input ────────────────
[[ -z "${input_dir:-}" ]] && echo "❌ Missing --input_dir" && usage
[[ -z "${r1_suffix:-}" ]] && echo "❌ Missing --r1_suffix" && usage


mode=$(echo "$mode" | tr '[:upper:]' '[:lower:]')
platform=$(echo "$platform" | tr '[:upper:]' '[:lower:]')

if [[ "$mode" != "pe" && "$mode" != "se" ]]; then
  echo "❌ --mode must be SE or PE (case-insensitive)"; usage
fi
if [[ "$mode" == "pe" && -z "${r2_suffix:-}" ]]; then
  echo "❌  --r2_suffix is required in PE mode"; usage
fi

# Correct precedence: fail only when file is missing
[[ -f "$primer_file" ]] || { echo "❌ not found $primer_file"; usage; }

# ─────────────── Logging Functions ─────────────
log() {
  echo "$(date '+[%F %T]') $*"
}
elapsed() {
  local s=$1; local e=$(date +%s)
  printf "Elapsed time: %02d:%02d:%02d\n" $(( (e-s)/3600 )) $(( ((e-s)%3600)/60 )) $(( (e-s)%60 ))
}

echo '''
          🧬 DADA2 Amplicon Pipeline
╭──────────────────────────────────────────────╮
│  Raw Reads  → fastp  →  cutadapt  →  dada2   │
╰──────────────────────────────────────────────╯
'''

# ─────────────── Step 0: fastp ─────────────────
log "🧼 Step 1: stat fq statisics using seqkit stats"
start_t=$(date +%s)
fqfiles=$(find "$input_dir" -type f \( -name "*$r1_suffix" -o -name "*$r2_suffix" \))

[[ -z "$fqfiles" ]] && { echo "No matching files found."; exit 1; }
seqkit stats -j "$threads" $fqfiles | sed "s|$input_dir\/||;s|$r1_suffix||;s|$r2_suffix||" > seqkit.stat.tsv
if [ $? -eq 0 ]; then
  echo "--------------------- seqkit finished. $(elapsed $start_t)"
else
  echo "seqkit error"
  exit 1
fi

echo ""
log "🧼Step 2: QC using fastp"
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

# Check failed
sample_count=$(printf '%s\n' $sample_list | sed '/^$/d' | wc -l)
fastp_finished_count=$(wc -l < fastp.rush.finished)

if (( $sample_count != $fastp_finished_count )) ; then
	echo "Sample count: $sample_count, fastp finished $fastp_finished_count"
    find 01_fastp -maxdepth 2 -name "*$r1_suffix" -exec basename {} \; | sed "s/$r1_suffix//" | grep -w -v -f - <(echo $sample_list | tr ' ' '\n') | sort -u > fastp.rush.failed.list
    echo "  ❌ fastp failed";
    exit 1
fi

# Summarize logs
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
  sample=$sample
  awk -v sample="$sample" "$awk_cmd" "$f" >> fastp.filter.tsv
done
echo "fastp resummary -> fastp.filter.tsv"
echo "--------------------- fastp finished. $(elapsed $start_t)"


# ─────────────── Step 3: cutadapt ──────────────
echo ""
log "✂️ Step 3: cutadapt primer trimming"
start_t=$(date +%s)

mkdir -p 02_cutadapt

f_primers=$(grep '^forward' "$primer_file" | while read a b c d; do echo "-g ${b}=^${c}";done | xargs)
r_primers=$(grep '^reverse' "$primer_file" | while read a b c d; do echo "-G ${b}=^${c}";done | xargs)
fr_primers=$(grep -e '^forward' -e '^reverse' "$primer_file" | while read a b c d; do echo "-g ${b}=^${c}";done | xargs)

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

# Check failed
cutadapt_finished_count=$(wc -l < cutadapt.rush.finished)

if (( $sample_count != $cutadapt_finished_count )) ; then
    find 02_cutadapt -maxdepth 2 -name "*$r1_suffix" -exec basename {} \; | sed "s/$r1_suffix//" | grep -w -v -f - <(echo $sample_list | tr ' ' '\n') | sort -u > cutadapt.rush.failed.list
    echo "  ❌ cutadapt failed";
    exit 1
fi

#echo "    summarizing 16S rRNA gene primer use"
if [[ $mode == "pe" ]];then
    summarize_cutadapt.py -d 02_cutadapt/ -m PE -t $threads
else
    summarize_cutadapt.py -d 02_cutadapt/ -m SE -t $threads
fi
if [ $? -ne 0 ]; then
  echo "    ❌ summarize_cutadapt.py"
  exit 1
fi
echo "--------------------- cutadapt finished. $(elapsed $start_t)"



# ─────────────── Step 4: DADA2 ────────────────
echo ""
log "🧬 Step 4: dada2.R"
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
    log "❌ dada2.R failed"
    exit 1
  fi
  log "$(elapsed $start_t)"
else
  log "sbatch dada2.slurm.sh"
  sbatch dada2.slurm.sh
fi

log "🧬 step check, should pe -> se?"
if [[ !-f 03_dada2/track.summary.tsv ]];then
  echo "    dada2.R error"
  exit 1
fi
amplicon_reads_lost_check.sh -i 03_dada2/track.summary.tsv

if [ $? -eq 0 ]; then
  log "amplicon_reads_lost_check.sh finished"
else
  echo "amplicon_reads_lost_check.sh error"
  exit 1
fi

if [[ -f 03_dada2/suggestion.pe2se.note ]];then
	echo "    PE -> SE"
	dada2.R -i 02_cutadapt --output_dir 03_dada2 --mode SE --reads1_suffix $r1_suffix --threads $threads --platform $platform
fi

[[ -f seqtab.nochim.rds ]] || echo "Error: ❌ dada2 failed"; exit 1
[[ -f track.summary.tsv ]] || echo "Error: ❌ dada2 failed"; exit 1
echo "--------------------- dada2 finished. $(elapsed $start_t)"

log "🧬 cleanup 00_fq 01_fastp 02_cutadapt"

rm -rf 00_fq 01_fastp 02_cutadapt
log "dd2_pipeline finished."
touch dd2_finished.note
exit 0