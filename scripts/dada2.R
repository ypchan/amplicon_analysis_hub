#！/usr/bin/env bash

set -euo pipefail

# ---- Logging function ----
log() {
    local info="$1"
    echo -e "\033[1;32m[$(date +'%Y-%m-%d %H:%M:%S')]\033[0m $info"
}
# ---- Check if required commands exist ----
check_dependencies() {
    local tools=("prefetch" "ParaFly" "awk" "cut" "sort" "tail" "mkdir" "wc" "mv")
    local missing=0

    for cmd in "${tools[@]}"; do
        if ! command -v "$cmd" &>/dev/null; then
            echo "❌ Required command not found: $cmd"
            missing=1
        fi
    done

    if [[ $missing -eq 1 ]]; then
        echo "🛑 Please install the missing tools before running this script."
        exit 1
    fi
}

check_dependencies

# ---- Argument check ----
if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <bioproject2sra.tsv>"
    exit 1
fi

bioproject2sra_tsv="$1"

if [[ ! -f $bioproject2sra_tsv ]]; then
    echo "❌ Input file '$bioproject2sra_tsv' not found."
    exit 1
fi

# ---- Extract unique BioProject IDs ----
project_lst=$(tail -n +2 "$bioproject2sra_tsv" | cut -f1 | sort -u)

# ---- Create output folders for each project ----
log "📁 Creating '00_fq' directories for each project..."
for project in $project_lst; do
    mkdir -p "${project}/00_fq"
done

# ---- Generate sra.list file for each project ----
log "📝 Generating 'sra.list' for each project..."
while read -r project; do
    awk -F '\t' -v pid="$project" '$1 == pid {print $2}' "$bioproject2sra_tsv" > "${project}/sra.list"
done <<< "$project_lst"

# ---- Download SRA using prefetch and ParaFly ----
log "⬇️  Starting SRA downloads using ParaFly..."
while read -r project; do
    log "🚀 Processing project: $project"

    cd "$project"

    # Generate prefetch command list
    awk '{print "prefetch " $1 " -O 00_fq"}' sra.list > prefetch.sh
    log "📄 Created prefetch.sh with $(wc -l < prefetch.sh) jobs"

    # Run ParaFly with automatic retries for failed downloads
    retry=0
    while [[ ! -f prefetch.sh.failed || -s prefetch.sh.failed ]]; do
        ((retry++))
        log "🔁 Running ParaFly (attempt $retry)..."
        ParaFly -c prefetch.sh -CPU 5 -failed_cmds prefetch.sh.failed
        mv -f prefetch.sh.failed prefetch.sh || true
    done

    log "✅ Finished project: $project"
    cd - > /dev/null
done <<< "$project_lst"

log "🎉 All SRA downloads completed successfully."
