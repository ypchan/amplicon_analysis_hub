#!/usr/bin/env bash

# ───────────────────────────────────────────────
# 🧪 DADA2 Amplicon Pipeline (Full Auto Version)
# Author: yanpengch@qq.com
# Description: From raw FASTQ to ASV table using fastp, cutadapt and DADA2
# ───────────────────────────────────────────────

# ──────────────── Usage Function ────────────────
usage() {
  cat <<EOF
dada2_pipeline.sh: Process amplicon sequencing data with DADA2.

Steps:
    1. fastp filtering & QC
    2. Primer detection and trimming (cutadapt)
    3. DADA2 denoising & ASV generation

Usage: dada2_pipeline.sh [options]

Required options:
    -i, --input_dir           Input directory with raw FASTQ files
    -1, --r1_suffix           R1 FASTQ suffix (e.g. _R1.fq.gz)
Optional options:
    -2, --r2_suffix           R2 suffix (required if mode=PE)
    -t, --threads             Number of threads (default: 4)
    -m, --mode                SE or PE (default: PE)
    -p, --platform            Platform: illumina, 454, iontorrent (default: illumina)
    -s, --slurm              Submit DADA2 job via SLURM (default: false)
    -h, --help                Show this message
EOF
  exit 1
}

# ─────────────── Default Parameters ─────────────
threads=4
mode="PE"
platform="illumina"

# ─────────────── Parse Arguments ────────────────
ARGS=$(getopt -o i:1:2:t:m:p:hs --long input_dir:,r1_suffix:,r2_suffix:,threads:,mode:,platform:,help,slurm -n 'dada2_pipeline.sh' -- "$@")
[[ $? -ne 0 ]] && usage
eval set -- "$ARGS"

while true; do
  case "$1" in
    -i|--input_dir)   input_dir="$2"; shift 2 ;;
    -1|--r1_suffix)   r1_suffix="$2"; shift 2 ;;
    -2|--r2_suffix)   r2_suffix="$2"; shift 2 ;;
    -t|--threads)     threads="$2"; shift 2 ;;
    -m|--mode)        mode="$2"; shift 2 ;;
    -p|--platform)    platform="$2"; shift 2 ;;
    -s|--slurm)       slurm=true; shift ;;
    -h|--help)        usage ;;
    --) shift; break ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

# ─────────────── Validate Input ────────────────
mode=$(echo "$mode" | tr '[:upper:]' '[:lower:]')
[[ -z "$input_dir" ]] && echo "❌ Missing --input_dir" && usage
[[ -z "$r1_suffix" ]] && echo "❌ Missing --r1_suffix" && usage
[[ "$mode" == "pe" && -z "$r2_suffix" ]] && echo " ❌ --r2_suffix is required in PE mode" && usage

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
│   Raw Reads   →   QC   →   Trim   →   ASV    │
│      │            │         │          │     │
│     🔬           🧼        ✂️        📊    │
╰──────────────────────────────────────────────╯
'''

# ─────────────── Step 1: fastp ─────────────────
log "🧼 Step 1: stat fq statisics using seqkit stats"
start_t=$(date +%s)
if [[ "$mode" == "se" ]]; then
    seqkit stats -j "$threads" $(ls "$input_dir"/*"$r1_suffix") | sed "s/$input_dir\///;s/$r1_suffix//" > seqkit.stat.tsv
else
    seqkit stats -j "$threads" $(ls "$input_dir"/*"$r1_suffix"; ls "$input_dir"/*"$r2_suffix" ) | sed "s/$input_dir\///;s/$r1_suffix//;s/$r2_suffix//" > seqkit.stat.tsv
fi

mkdir -p 01_fastp

log "🧼 fastp to filter bad reads"
sample_list=$(find "$input_dir" -name "*$r1_suffix" -exec basename {} \; | sed "s/$r1_suffix//")

if [[ "$mode" == "pe" ]]; then
  echo "$sample_list" | rush -j "$threads" -v r1="$r1_suffix",r2="$r2_suffix",input_dir="$input_dir" -c --eta --succ-cmd-file fastp.rush.finished 'fastp -i {input_dir}/{1}{r1} -I {input_dir}/{1}{r2} -o 01_fastp/{1}{r1} -O 01_fastp/{1}{r2} --thread 1 --length_required 100 --n_base_limit 0 --cut_tail --qualified_quality_phred 20 --unqualified_percent_limit 20 --html /dev/null --json /dev/null &> 01_fastp/{1}.fastp.log'
else
  echo "$sample_list" | rush -j "$threads" -v r1="$r1_suffix",input_dir="$input_dir" -c --eta --succ-cmd-file fastp.rush.finished 'fastp -i {input_dir}/{1}{r1} -o 01_fastp/{1}{r1} --thread 1 --length_required 100 --n_base_limit 0 --cut_tail --qualified_quality_phred 20 --unqualified_percent_limit 20 --html /dev/null --json /dev/null &> 01_fastp/{1}.fastp.log'
fi

# Summarize logs
log "🧼 summary fastp results"
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

log "$(elapsed $start_t)"

# ─────────────── Step 2: cutadapt ──────────────
log "✂️  Step 2: cutadapt primer trimming"
start_t=$(date +%s)
mkdir -p 02_cutadapt
primer_file="/home/data/t170527/database/16s_primer.tsv"

f_primers=$(grep '^forward' "$primer_file" | while read a b c d; do echo "-g ${b}=^${c}";done | xargs)
r_primerss=$(grep '^reverse' "$primer_file" | while read a b c d; do echo "-G ${b}=^${c}";done | xargs)
fr_primers=$(grep -e '^forward' -e '^reverse' "$primer_file" | while read a b c d; do echo "-g ${b}=^${c}";done | xargs)

if [[ "$mode" == "pe" ]]; then
  cutadapt_opts="$f_primers $r_primerss --revcomp -j 1"
  echo "$sample_list" | rush -j "$threads" -v r1="$r1_suffix",r2="$r2_suffix",opt="${cutadapt_opts}" -c --eta --succ-cmd-file cutadapt.rush.finished 'cutadapt {opt} -o 02_cutadapt/{1}{r1} -p 02_cutadapt/{1}{r2} 01_fastp/{1}{r1} 01_fastp/{1}{r2} &> 02_cutadapt/{1}.cutadapt.log'
else
  cutadapt_opts="$fr_primers --revcomp -j 1"
  echo "$sample_list" | rush -j "$threads" -v r1="$r1_suffix",opt="${cutadapt_opts}" -c --eta --succ-cmd-file cutadapt.rush.finished 'cutadapt {opt} -o 02_cutadapt/{1}{r1} 01_fastp/{1}{r1} &> 02_cutadapt/{1}.cutadapt.log'
fi
log "$(elapsed $start_t)"

log "summarize 16S rRNA gene primer use"
if [[ $mode == "pe" ]];then
    summarize_cutadapt.py -d 02_cutadapt/ -m PE -t 24
else
    summarize_cutadapt.py -d 02_cutadapt/ -m SE -t 24
fi

# ─────────────── Step 3: DADA2 ────────────────
log "🧬 Step 3: Run DADA2 (R script)"
mkdir -p 03_dada2
dd_cmd="dada2.R -i 02_cutadapt --output_dir 03_dada2 --mode $mode --reads1_suffix $r1_suffix --threads $threads --platform $platform"
[[ "$mode" == "pe" ]] && dd_cmd+=" --reads2_suffix $r2_suffix"

cat > dada2.slurm.sh <<EOF
#!/bin/bash
#SBATCH --job-name=dada2
#SBATCH --partition=cn
#SBATCH --output=%x.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=$threads
#SBATCH --mem=500G
#SBATCH --time=8:00:00

exec 2>&1
source /home/software/miniconda3/etc/profile.d/conda.sh
conda activate dada2
$dd_cmd
EOF

if [[ "$slurm" != true ]]; then
  log " run DADA2: $dd_cmd"
  rm -f dada2.slurm.sh
  if ! eval "$dd_cmd" 2>&1 | tee dd2.log; then
    log "❌ dada2.R failed"
    exit 1
  fi
  exit 0
else
  log "sbatch dada2.slurm.sh"
  if [[ ! -x dada2.slurm.sh ]]; then
    echo "❌ dada2.slurm.sh is not executable, please check permissions"
    exit 1
  fi
  sbatch dada2.slurm.sh
fi

exit 0