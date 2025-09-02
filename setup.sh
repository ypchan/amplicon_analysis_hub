# gh repo clone ypchan/5M16S #not public yet
# if failed, download this repo manually and place it in the INSTALL_DIR
# cd 5M16S || exit 1

# colorful output log
if [[ -t 1 && -z "${NO_COLOR:-}" && "${TERM:-}" != "dumb" ]]; then
  C0=$'\033[0m'; Cg=$'\033[1;32m'; Cy=$'\033[1;33m'; Cr=$'\033[1;31m'; Cb=$'\033[1;34m'
else
  C0=""; Cg=""; Cy=""; Cr=""; Cb=""
fi
ts(){ date '+[%F %T]'; }
log(){  printf '%s [%sINFO%s] %s\n'  "$(ts)" "$Cb" "$C0" "$*"; }
ok(){   printf '%s [ %sOK%s ] %s\n'   "$(ts)" "$Cg" "$C0" "$*"; }
err(){  printf '%s [%sERR%s ] %s\n'  "$(ts)" "$Cr" "$C0" "$*" >&2;exit 1; }


INSTALL_HOME=$(realpath .)
log "Installing to $INSTALL_HOME"
# ------ biosoftwares dependency
REQUIREMENTS=("python3" "Rscript" "blastn" "makeblastdb" \
            "cd-hit" "fastp" "seqkit" "cutadapt" "rush" "dos2unix")

for cmd in "${REQUIREMENTS[@]}"; do
    command -v "$cmd" &>/dev/null
    if [[ $? -ne 0 ]]; then
        err "$cmd is required, but it's not installed"
    else
        ok "$cmd ok"
    fi
done

if [[ ! $(cutadapt --version) == "5.1" ]];then
    err "cutadapt version 5.1 is required, but found $(cutadapt --version), please install cutadapt v5.1"
else
    ok "cutadapt ok"  
fi

## ---- R packages dependency
Rscript -e 'quit(status = !requireNamespace("dada2", quietly = TRUE))' 
if [[ $? -eq 0 ]]; then 
    ok "dada2 ok"
else
    err "dada2 not installed, try installing it"
fi

Rscript -e 'quit(status = !requireNamespace("getopt", quietly = TRUE))' 
if [[ $? -eq 0 ]]; then 
    ok "getopt ok"
else
    err "getopt not installed, try installing it"
fi 

## ---- python packages dependency scripts
python -c "import sys, importlib.util as u; sys.exit(0 if u.find_spec('pandas') else 1)" 
if [[ $? -eq 0 ]]; then 
    ok "pandas ok"
else
    err "pandas not installed, try installing it"
fi 

# -- make scripts executable
dos2unix scripts/*  
chmod 755 scripts/*.py scripts/*.R scripts/*.sh

if [ $? -ne 0 ]; then
    log "Please execute the following commands:"
    mkdir -p "$HOME/bin"
    #echo "export PATH=$HOME/bin:$PATH >>~/.bashrc"
    #echo "source ~/.bashrc"
fi
ls scripts | while read a;do 
    rm -f "$HOME/bin/$a" && ln -s "$INSTALL_HOME"/scripts/"$a" "$HOME/bin/$a"
    if [[ $? -ne 0 ]]; then 
        err "Failed to link $a to $HOME/bin, please check"
    else
        ok "$a linked to $HOME/bin"
    fi
done

ls scripts | while read a;do 
    command -v "$a" &>/dev/null
    if [[ $? -ne 0 ]]; then 
        err "$a not in PATH, please check"
    else
        ok "$a ok"
    fi
done

#  -- prepare databases
log "Constructing 16S rRNA gene blast db"
DB_NOTE="blastn.$(blastn -version | head -n 1 |awk '{print $2}')"
if [ -f "data/arc_bac_16s_blastDB/$DB_NOTE" ];then
    ok "blast db already prepared, skip"
else
    cd data/arc_bac_16s_blastDB
    bash work.sh
fi
cd "$INSTALL_HOME"

# -- set default data paths in dd2_pipeline.sh
PRIMER_FILE="$INSTALL_HOME/data/16s_primer.tsv"
if [ ! -f "$PRIMER_FILE" ];then
    err "$PRIMER_FILE not found"
fi
sed -i "s|^PRIMER_FILE=.*$|PRIMER_FILE=\"$PRIMER_FILE\"|" scripts/dd2_pipeline.sh

BLASTDB_16S="$INSTALL_HOME/data/arc_bac_16s_blastDB/arch_bac_16s_ref_90.nhr"
if [[ ! -f "$BLASTDB_16S" ]];then
    err "$BLASTDB_16S not found, please check"
fi

BLASTDB_16S=${BLASTDB_16S%.ndb}
sed -i "s|^BLASTDB_16S=.*$|BLASTDB_16S=\"$BLASTDB_16S\"|" scripts/dd2_pipeline.sh
log "Installation completed."
echo ""
log "Take a quick look"
dd2_pipeline.sh --help