#!/usr/bin/env bash

# Fast, data-free syntax and CLI smoke tests for amplicon_analysis_hub.

set -Eeuo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
root_dir="$(cd -- "$script_dir/.." && pwd -P)"

for script in "$root_dir/setup.sh" "$script_dir/amplicon_analysis" "$script_dir"/*.sh; do
  bash -n "$script"
done
python3 -m compileall -q "$script_dir"
for script in "$script_dir"/*.R; do
  Rscript -e 'parse(file=commandArgs(TRUE)[1])' "$script" >/dev/null
done

bash "$script_dir/amplicon_analysis" --help >/dev/null
python3 "$script_dir/is_16s_amplicon.py" --help >/dev/null
python3 "$script_dir/fastq_dispatcher.py" --help >/dev/null
python3 "$script_dir/get_ena_fq_url_by_sra.py" --help >/dev/null
python3 "$script_dir/summarize_cutadapt.py" --help >/dev/null
python3 "$script_dir/unify_fq_suffix.py" --help >/dev/null
python3 "$script_dir/ontology_train_cv.py" --help >/dev/null
python3 "$script_dir/ontology_infer.py" --help >/dev/null
Rscript "$script_dir/dada2.R" --help >/dev/null
Rscript "$script_dir/dada2.R" -M its -P pacbio_ccs -m se --print_profile >/dev/null
Rscript "$script_dir/asv_annotator.R" --help >/dev/null
Rscript "$script_dir/count_abundance.R" --help >/dev/null
Rscript "$script_dir/infer_16s_regions.R" --help >/dev/null

printf 'All syntax and CLI smoke tests passed.\n'
