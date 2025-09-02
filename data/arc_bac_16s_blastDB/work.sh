wget -c https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Archaea/archaea.16SrRNA.fna.gz
wget -c https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz
gunzip *
sed 's/^>/>archaea_/' archaea.16SrRNA.fna > arch_16s_ref.fna
sed 's/^>/>bacteria_/' bacteria.16SrRNA.fna > bac_16s_ref.fna
cat arch_16s_ref.fna bac_16s_ref.fna > arch_bac_16s_ref.fna
# Add NR version
cd-hit -i arch_bac_16s_ref.fna -o arch_bac_16s_ref_90.fna -c 0.90 -T 4 -M 16000
# Build blast db
makeblastdb -in arch_bac_16s_ref_90.fna -dbtype nucl -input_type fasta -out arch_bac_16s_ref_90