#!/usr/bin/env bash

# Install amplicon_analysis_hub command symlinks and validate dependencies.
# Source files are never rewritten with machine-specific absolute paths.

set -Eeuo pipefail

VERSION="2.0.0"
ROOT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
PREFIX="${HOME}/.local/bin"
CHECK_ONLY=false
BUILD_16S_DB=false

usage() {
  cat <<'EOF'
Usage: bash setup.sh [options]

Options:
  --prefix DIR       Command symlink directory (default: ~/.local/bin)
  --check-only       Validate dependencies without installing symlinks
  --build-16s-db     Download RefSeq 16S loci and build the optional BLAST DB
  -h, --help         Show this help
  -V, --version      Show version

Core dependencies:
  bash >=4.3, Python >=3.9, R, DADA2, getopt (R), fastp, Cutadapt >=4.1,
  seqkit, and standard POSIX/GNU utilities (including cmp, cksum, and realpath).

Optional dependencies:
  BLAST+ and cd-hit-est for the 16S content screen/database build; vsearch and
  tidyverse packages for region/abundance helpers; scikit-learn,
  sentence-transformers, pandas, pyarrow and joblib for ontology helpers.
EOF
}

while (($#)); do
  case "$1" in
    --prefix) [[ $# -ge 2 ]] || { echo "--prefix requires a value" >&2; exit 2; }; PREFIX="$2"; shift 2 ;;
    --check-only) CHECK_ONLY=true; shift ;;
    --build-16s-db) BUILD_16S_DB=true; shift ;;
    -h|--help) usage; exit 0 ;;
    -V|--version) printf 'setup.sh %s\n' "$VERSION"; exit 0 ;;
    *) printf 'Unknown option: %s\n' "$1" >&2; usage >&2; exit 2 ;;
  esac
done

log() { printf '[INFO] %s\n' "$*"; }
ok() { printf '[ OK ] %s\n' "$*"; }
warn() { printf '[WARN] %s\n' "$*" >&2; }
failures=0
require_command() {
  if command -v "$1" >/dev/null 2>&1; then ok "$1: $(command -v "$1")"; else warn "missing command: $1"; ((failures+=1)); fi
}

for command in bash getopt python3 Rscript fastp cutadapt seqkit awk sed find sort gzip rev tr cmp cksum realpath; do
  require_command "$command"
done

if ((BASH_VERSINFO[0] < 4 || (BASH_VERSINFO[0] == 4 && BASH_VERSINFO[1] < 3))); then
  warn "bash >=4.3 is required (found $BASH_VERSION)"
  ((failures+=1))
else
  ok "bash $BASH_VERSION"
fi

if command -v python3 >/dev/null 2>&1; then
  if python3 -c 'import sys; raise SystemExit(sys.version_info < (3, 9))'; then
    ok "$(python3 --version)"
  else
    warn "Python >=3.9 is required"
    ((failures+=1))
  fi
fi

if command -v Rscript >/dev/null 2>&1; then
  for package in dada2 getopt; do
    if Rscript -e "quit(status=!requireNamespace('$package', quietly=TRUE))"; then
      ok "R package: $package"
    else
      warn "missing R package: $package"
      ((failures+=1))
    fi
  done
fi

if command -v cutadapt >/dev/null 2>&1; then
  cutadapt_version="$(cutadapt --version | head -n 1)"
  cutadapt_major="${cutadapt_version%%.*}"
  if [[ "$cutadapt_major" =~ ^[0-9]+$ ]] && ((cutadapt_major >= 4)); then
    ok "cutadapt $cutadapt_version"
  else
    warn "Cutadapt >=4.1 is recommended; found $cutadapt_version"
    ((failures+=1))
  fi
fi

if [[ "$BUILD_16S_DB" == true ]]; then
  for command in blastn makeblastdb cd-hit-est wget; do require_command "$command"; done
  if ((failures == 0)); then
    log "Building the optional 16S BLAST database"
    bash "$ROOT_DIR/data/arc_bac_16s_blastDB/work.sh"
  fi
elif [[ ! -f "$ROOT_DIR/data/arc_bac_16s_blastDB/arch_bac_16s_ref_90.nhr" ]]; then
  warn "optional 16S BLAST DB is absent; run: bash setup.sh --build-16s-db"
fi

if ((failures > 0)); then
  warn "$failures required check(s) failed; see README.md for installation guidance"
  exit 1
fi

if [[ "$CHECK_ONLY" == false ]]; then
  mkdir -p -- "$PREFIX"
  for script in "$ROOT_DIR"/scripts/*; do
    [[ -f "$script" ]] || continue
    case "$script" in */amplicon_analysis|*.py|*.R|*.sh) ;; *) continue ;; esac
    chmod 755 -- "$script"
    ln -sfn -- "$script" "$PREFIX/${script##*/}"
  done
  ok "commands linked into $PREFIX"
  case ":$PATH:" in
    *":$PREFIX:"*) ;;
    *) warn "$PREFIX is not in PATH; add: export PATH=\"$PREFIX:\$PATH\"" ;;
  esac
fi

log "amplicon_analysis_hub setup complete"
bash "$ROOT_DIR/scripts/amplicon_analysis" --version
