#!/usr/bin/env bash

set -Eeuo pipefail
root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
work="$(mktemp -d "${TMPDIR:-/tmp}/amplicon-hub-tests.XXXXXX")"
trap 'rm -rf -- "$work"' EXIT

printf '@read.a.1\nACGT\n+\nIIII\n@read.a.2\nTGCA\n+\nIIII\n' > "$work/mixed.fastq"
bash "$root/scripts/split_fq12.sh" -i "$work/mixed.fastq" -1 "$work/r1.fastq" -2 "$work/r2.fastq" >/dev/null
[[ "$(wc -l < "$work/r1.fastq")" -eq 4 && "$(wc -l < "$work/r2.fastq")" -eq 4 ]]

mkdir "$work/rename"
cp "$work/r1.fastq" "$work/rename/sample_1.fastq"
python3 "$root/scripts/unify_fq_suffix.py" -i "$work/rename" -1 _1.fastq -f _R1.fastq.gz >/dev/null
gzip -t "$work/rename/sample_R1.fastq.gz"

printf '\tinput\tfiltered\tdenoisedF\tdenoisedR\tmerged\tnonchim\nS1\t100\t90\t85\t80\t40\t35\nS2\t100\t90\t85\t80\t80\t75\n' > "$work/track.tsv"
bash "$root/scripts/amplicon_reads_lost_check.sh" -i "$work/track.tsv" -o "$work/retention.tsv" >/dev/null
grep -q $'review_pe_as_se\tYES' "$work/retention.summary.tsv"
[[ -f "$work/suggestion.pe2se.note" && ! -e "$work/suggestion.is_pe.note" ]]
bash "$root/scripts/amplicon_reads_lost_check.sh" -i "$work/track.tsv" -o "$work/retention.tsv" \
  --retained-threshold 30 >/dev/null
[[ -f "$work/suggestion.is_pe.note" && ! -e "$work/suggestion.pe2se.note" ]]

mkdir "$work/bin"
printf '#!/usr/bin/env bash\nawk '\''/^>/{print substr($0,2) "\\tbacteria_mock\\t99\\t100\\t100"}'\''\n' > "$work/bin/blastn"
chmod 755 "$work/bin/blastn"
: > "$work/mockdb.nhr"
PATH="$work/bin:$PATH" python3 "$root/scripts/is_16s_amplicon.py" "$work/r1.fastq" \
  --db "$work/mockdb" --format tsv --output "$work/screen.tsv" --out-format tsv >/dev/null
grep -q $'r1.fastq\t1\t0\t1\t100.000\tYES' "$work/screen.tsv"

mkdir "$work/cutadapt"
printf "Total read pairs processed:              100\n=== First read: Adapter 'FWD' ===\nSequence: ACGT; Type: regular 5'; Length: 4; Trimmed: 80 times\n=== Second read: Adapter 'REV' ===\nSequence: TGCA; Type: regular 5'; Length: 4; Trimmed: 75 times\n" > "$work/cutadapt/sample.cutadapt.log"
python3 "$root/scripts/summarize_cutadapt.py" -d "$work/cutadapt" -m PE -o "$work" >/dev/null
grep -q $'sample\t100\tFWD' "$work/cutadapt_details.tsv"

mkdir "$work/dispatch_in"
cp "$work/r1.fastq" "$work/dispatch_in/SRR1_1.fastq.gz"
cp "$work/r2.fastq" "$work/dispatch_in/SRR1_2.fastq.gz"
printf 'run\tproject\tlayout\tplatform\nSRR1\tPRJ1\tPAIRED\tDNBSEQ-G400\n' > "$work/metadata.tsv"
python3 "$root/scripts/fastq_dispatcher.py" -m "$work/metadata.tsv" -f "$work/dispatch_in" \
  -o "$work/dispatched" --header --run-col 1 --bioproject-col 2 \
  --layout-col 3 --platform-col 4 --action copy -t 2 >/dev/null
[[ -f "$work/dispatched/PRJ1_pe_bgi/00_fq/SRR1_1.fastq.gz" ]]

Rscript -e 'x <- matrix(c(10,0,3,7), nrow=2, byrow=TRUE, dimnames=list(c("S1","S2"),c("ACGT","TGCA"))); saveRDS(x, commandArgs(TRUE)[1])' "$work/seqtab.rds"
printf '\tKingdom\tPhylum\tGenus\nACGT\tBacteria\tFirmicutes\tBacillus\nTGCA\tBacteria\tProteobacteria\t\n' > "$work/taxonomy.tsv"
Rscript "$root/scripts/count_abundance.R" -s "$work/seqtab.rds" -t "$work/taxonomy.tsv" \
  -o "$work/abundance" -r Genus >/dev/null
grep -q 'Unclassified_Proteobacteria' "$work/abundance/abundance_Genus_counts.tsv"

mkdir -p "$work/path_guard/input" "$work/path_guard/output/input"
touch "$work/path_guard/input/sample.fastq.gz" "$work/path_guard/output/input/sample.fastq.gz"
if bash "$root/scripts/amplicon_analysis" -i "$work/path_guard/input" \
  -o "$work/path_guard/input/results" -m se -1 .fastq.gz --primer-mode none \
  --fastp no --screen no -c none >/dev/null 2>&1; then
  echo "nested output path was not rejected" >&2
  exit 1
fi

# The project run directory may safely contain a read-only input subdirectory.
# Invoke through a symlink too, matching setup.sh, to verify that bundled data
# paths are resolved relative to the real script rather than the install prefix.
mkdir -p "$work/path_guard/bin"
ln -s "$root/scripts/amplicon_analysis" "$work/path_guard/bin/amplicon_analysis"
printf '#!/usr/bin/env bash\nexit 0\n' > "$work/path_guard/bin/seqkit"
printf '#!/usr/bin/env bash\nwhile (($#)); do\n  case "$1" in\n    -o|-O) printf stub > "$2"; shift 2 ;;\n    --json) printf "{}\\n" > "$2"; shift 2 ;;\n    *) shift ;;\n  esac\ndone\n' \
  > "$work/path_guard/bin/fastp"
printf '#!/usr/bin/env bash\nout=""\nwhile (($#)); do\n  case "$1" in\n    --output_dir|-o) out="$2"; shift 2 ;;\n    *) shift ;;\n  esac\ndone\nif [[ -n "$out" ]]; then mkdir -p -- "$out"; printf stub > "$out/seqtab.nochim.rds"; fi\n' \
  > "$work/path_guard/bin/Rscript"
chmod 755 "$work/path_guard/bin/seqkit" "$work/path_guard/bin/fastp" \
  "$work/path_guard/bin/Rscript"
PATH="$work/path_guard/bin:$PATH" "$work/path_guard/bin/amplicon_analysis" \
  -i "$work/path_guard/output/input" -o "$work/path_guard/output" -M its \
  -m se -1 .fastq.gz --skip-cutadapt --fastp yes --screen no -c none \
  --cleanup none >/dev/null
[[ -f "$work/path_guard/output/input/sample.fastq.gz" ]]
[[ -s "$work/path_guard/output/01_fastp/sample.fastq.gz" ]]
[[ -s "$work/path_guard/output/amplicon_analysis_hub.finished" ]]
grep -Fq $'primer_file\t'"$root/data/its_primer.tsv" "$work/path_guard/output/run_parameters.tsv"
[[ ! -e "$work/path_guard/output/02_cutadapt" ]]

printf 'Helper integration tests passed.\n'
