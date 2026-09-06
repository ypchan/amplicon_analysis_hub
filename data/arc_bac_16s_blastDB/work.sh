#!/usr/bin/env bash

# Build a non-redundant bacterial/archaeal 16S BLAST database from RefSeq.

set -Eeuo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
threads="${THREADS:-4}"
memory_mb="${MEMORY_MB:-0}"

for command in wget gzip awk cd-hit-est makeblastdb blastn; do
  command -v "$command" >/dev/null 2>&1 || { printf 'Missing command: %s\n' "$command" >&2; exit 127; }
done
[[ "$threads" =~ ^[1-9][0-9]*$ ]] || { echo "THREADS must be a positive integer" >&2; exit 2; }

arch_url="https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Archaea/archaea.16SrRNA.fna.gz"
bac_url="https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz"
wget --continue --directory-prefix "$script_dir" "$arch_url" "$bac_url"

combined="$script_dir/arch_bac_16s_ref.fna"
gzip -dc "$script_dir/archaea.16SrRNA.fna.gz" \
  | awk '/^>/{sub(/^>/, ">archaea_")} {print}' > "$combined"
gzip -dc "$script_dir/bacteria.16SrRNA.fna.gz" \
  | awk '/^>/{sub(/^>/, ">bacteria_")} {print}' >> "$combined"

clustered="$script_dir/arch_bac_16s_ref_90.fna"
cd-hit-est -i "$combined" -o "$clustered" -c 0.90 -n 8 -aS 0.8 \
  -T "$threads" -M "$memory_mb"
makeblastdb -in "$clustered" -dbtype nucl -parse_seqids \
  -out "$script_dir/arch_bac_16s_ref_90"

{
  printf 'built_at\t%s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')"
  printf 'blast_version\t%s\n' "$(blastn -version | head -n 1)"
  printf 'identity_cluster\t0.90\n'
  printf 'coverage_shorter\t0.8\n'
  printf 'archaea_source\t%s\n' "$arch_url"
  printf 'bacteria_source\t%s\n' "$bac_url"
} > "$script_dir/database_manifest.tsv"
printf 'Built database: %s\n' "$script_dir/arch_bac_16s_ref_90"
