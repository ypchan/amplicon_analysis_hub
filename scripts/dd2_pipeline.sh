#!/usr/bin/env bash

# Backward-compatible entry point. New workflows should use
# amplicon_analysis, which supports both 16S and ITS profiles.

set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
printf 'WARNING: dd2_pipeline.sh is deprecated; use amplicon_analysis.\n' >&2
exec bash "$script_dir/amplicon_analysis" "$@"
