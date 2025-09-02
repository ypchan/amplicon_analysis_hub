INSTALL_DIR=""

cd "$INSTALL_DIR" || exit 1

gh repo clone ypchan/5M16S
# if failed, download this repo manually and place it in the INSTALL_DIR
cd 5M16S || exit 1

# ------ biosoftwares dependency
REQUIREMENTS=("python3" "Rscript" "blastn" "makeblastdb" \
            "cd-hit" "fastp" "seqtk" "cutadapt" "rush")

for cmd in "${REQUIREMENTS[@]}"; do
    command -v "$cmd" >/dev/null 2>&1 || { echo >&2 "$cmd is required but it's not installed"; exit 1; }
done    
if [[ $(cutadapt --version) == "5.1" ]];then
    echo "cutadapt version 5.1 is required but found $(cutadapt --version), please install cutadapt v5.1"
    exit 1
fi

## ---- R packages dependency
Rscript -e 'quit(status = !requireNamespace("dada2", quietly = TRUE))' \
  && echo "dada2 installed" || echo "dada2 NOT installed"
[ $? -ne 0 ] && exit 1
Rscript -e 'quit(status = !requireNamespace("optparse", quietly = TRUE))' \
  && echo "optparse installed" || echo "optparse NOT installed"
[ $? -ne 0 ] && exit 1
Rscript -e 'quit(status = !requireNamespace("getopt", quietly = TRUE))' \
  && echo "getopt installed" || echo "getopt NOT installed"
[ $? -ne 0 ] && exit 1 

## ---- python packages dependency scripts
python -c "import sys, importlib.util as u; sys.exit(0 if u.find_spec('pandas') else 1)" \
  && echo "pandas installed" || echo "pandas NOT installed"
[ $? -ne 0 ] && exit 1

# -- make scripts executable
chmod 755 scripts/*.py scripts/*.R scripts/*.sh

if [ $? -ne 0 ]; then
    echo "Please execute the following commands:"
    mkdir -p "$HOME/bin"
    echo "export PATH=$HOME/bin:$PATH >>~/.bashrc"
    echo "source ~/.bashrc"
fi
ls scripts | while read a;do rm -f "$HOME/bin/$a";ln -s "$INSTALL_DIR"/5M16S/scripts/"$a" "$HOME/bin/$a";done

ls scripts | while read a;do command -v "$a" >/dev/null 2>&1 || { echo >&2 "$a is not executable, pls debug"; exit 1; }; done

#  -- prepare databases
cd data/arc_bac_16s_blastDB
bash work.sh

# -- set default data paths in dd2_pipeline.sh
PRIMER_FILE=$(realpath data/16s_primer.tsv)
if [ ! -f "$PRIMER_FILE" ];then
    echo "$PRIMER_FILE not found, please check"
    exit 1
fi
sed -i "s|^PRIMER_FILE=.*$|PRIMER_FILE=\"$PRIMER_FILE\"|" scripts/dd2_pipeline.sh

BLASTDB_16S=$(realpath data/arc_bac_16s_blastDB/arch_bac_16s_ref_90.ndb)
if [ ! -f "$BLASTDB_16S".nin ];then
    echo "$BLASTDB_16S not found, please check"
    exit 1
fi
BLASTDB_16S=${BLASTDB_16S%.ndb}
sed -i "s|^BLASTDB_16S=.*$|BLASTDB_16S=\"$BLASTDB_16S\"|" scripts/dd2_pipeline.sh
echo "Installation completed."
dd2_pipeline.sh --help



