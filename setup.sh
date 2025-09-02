# gh repo clone ypchan/5M16S #not public yet
# if failed, download this repo manually and place it in the INSTALL_DIR
# cd 5M16S || exit 1
INSTALL_HOME=$(realpath .)
# ------ biosoftwares dependency
REQUIREMENTS=("python3" "Rscript" "blastn" "makeblastdb" \
            "cd-hit" "fastp" "seqkit" "cutadapt" "rush")

for cmd in "${REQUIREMENTS[@]}"; do
    command -v "$cmd" &>/dev/null
    if [[ $? -ne 0 ]]; then
        echo "$cmd is required but it's not installed"
        exit 1
    else
        echo "$cmd ok"
    fi
done

if [[ ! $(cutadapt --version) == "5.1" ]];then
    echo "cutadapt version 5.1 is required but found $(cutadapt --version), please install cutadapt v5.1"
    exit 1
else
    echo "cutadapt ok"  
fi

## ---- R packages dependency
Rscript -e 'quit(status = !requireNamespace("dada2", quietly = TRUE))' 
if [[ $? -eq 0 ]]; then 
    echo "dada2 ok"
else
    echo "dada2 not installed, try installing it"
    exit 1
fi

Rscript -e 'quit(status = !requireNamespace("getopt", quietly = TRUE))' 
if [[ $? -eq 0 ]]; then 
    echo "getopt ok"
else
    echo "getopt not installed, try installing it"
fi 

## ---- python packages dependency scripts
python -c "import sys, importlib.util as u; sys.exit(0 if u.find_spec('pandas') else 1)" 
if [[ $? -eq 0 ]]; then 
    echo "pandas ok"
else
    echo "pandas not installed, try installing it"
fi 

# -- make scripts executable
chmod 755 scripts/*.py scripts/*.R scripts/*.sh

if [ $? -ne 0 ]; then
    echo "Please execute the following commands:"
    mkdir -p "$HOME/bin"
    echo "export PATH=$HOME/bin:$PATH >>~/.bashrc"
    echo "source ~/.bashrc"
fi
ls scripts | while read a;do 
    rm -f "$HOME/bin/$a" && ln -s "$INSTALL_HOME"/scripts/"$a" "$HOME/bin/$a"
    if [[ $? -ne 0 ]]; then 
        echo "Failed to link $a to $HOME/bin, please check"
        exit 1
    else
        echo "$a linked to $HOME/bin"
    fi
done

ls scripts | while read a;do 
    command -v "$a"
    if [[ $? -ne 0 ]]; then 
        echo "$a not in PATH, please check"
        exit 1
    else
        echo "$a ok"
    fi
done

#  -- prepare databases
DB_NOTE="blastn.$(blastn -version | head -n 1 |awk '{print $2}')"
if [ -f "data/arc_bac_16s_blastDB/$DB_NOTE" ];then
    echo "blast db already prepared, skip"
else
    cd data/arc_bac_16s_blastDB
    bash work.sh
fi
cd "$INSTALL_HOME"

# -- set default data paths in dd2_pipeline.sh
PRIMER_FILE="$INSTALL_HOME/data/16s_primer.tsv"
if [ ! -f "$PRIMER_FILE" ];then
    echo "$PRIMER_FILE not found"
    exit 1
fi
sed -i "s|^PRIMER_FILE=.*$|PRIMER_FILE=\"$PRIMER_FILE\"|" scripts/dd2_pipeline.sh

BLASTDB_16S="$INSTALL_HOME/data/arc_bac_16s_blastDB/arch_bac_16s_ref_90.ndb"
if [[ ! -f "$BLASTDB_16S" ]];then
    echo "$BLASTDB_16S not found, please check"
    exit 1
fi

BLASTDB_16S=${BLASTDB_16S%.ndb}
sed -i "s|^BLASTDB_16S=.*$|BLASTDB_16S=\"$BLASTDB_16S\"|" scripts/dd2_pipeline.sh
echo "Installation completed."
dd2_pipeline.sh --help