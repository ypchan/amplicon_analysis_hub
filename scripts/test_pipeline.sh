# test is_16s_amplicon.py
ls 00_fq/*_1.fastq.gz | is_16s_amplicon.py - --threads 4 --concurrent 3 --nreads 100 --out-format tsv --output is_16s.tsv

# seqkit
seqkit stats -j 20 00_fq/*.gz > seqkit.stats.tsv

# fastp
mkdir -p 01_fastp
ls 00_fq/ | sed 's/_1.fastq.gz//;s/_2.fastq.gz//' | sort -u | rush -j 20 --eta 'fastp -i 00_fq/{1}_1.fastq.gz -I 00_fq/{1}_2.fastq.gz -o 01_fastp/{1}_1.fastq.gz -O 01_fastp/{1}_2.fastq.gz --thread 1 --length_required 100 --n_base_limit 0 --cut_tail --qualified_quality_phred 20 --unqualified_percent_limit 20 --html /dev/null --json /dev/null &> 01_fastp/{1}.fastp.log'

# cutadapt
mkdir -p 02_cutadapt
PRIMER_FILE=
ls 01_fastp/ | sed 's/_1.fastq.gz//;s/_2.fastq.gz//' | sort -u | rush -j 20 --eta 'cutadapt -g ^CCTACGGGNGGCWGCAG...AGAGTTTGATCMTGGCTCAG -G ^GACTACHVGGGTATCTAATCC...TACGGYTACCTTGTTACGACT -o 02_cutadapt/{1}_1.fastq.gz -p 02_cutadapt/{1}_2.fastq.gz 01_fastp/{1}_1.fastq.gz 01_fastp/{1}_2.fastq.gz --discard-untrimmed --minimum-length 100 --max-n 0 --quality-cutoff 20 --trim-n --cores=0 &> 02_cutadapt/{1}.cutadapt.log'