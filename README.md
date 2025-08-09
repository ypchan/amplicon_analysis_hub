# ProkaAtlas ![GitHub](https://img.shields.io/badge/GitHub-Atlas-blue?logo=github)

A comprehensive pipeline for processing 16S rRNA gene amplicon data, including DADA2 tables, scripts, and results.




- [Software Dependencies](#software-dependencies)
- [Reference Data](#reference-data)
- [Workflow](#command-line-workflow)
    - [1. Select SRA Records](#1-select-sra-records-by-bioproject-accession)
    - [2. Exclude Non-16S Records & Download](#2-exclude-non-16s-sra-records-and-download)
    - [3. Convert SRA to FASTQ](#3-convert-sra-to-fastq)
    - [4. Downstream Processing](#4-downstream-processing-with-in-house-pipeline)


## Challenges
> 1: How to download all available 16S amplicon sequencing data as comprehensively as possible?

> 2: How to tankle different sequencing platforms? Roche 454 | Illumina | Ion Torrent | Pacbio?

> 3: To or not cut primers?

> 4: Different regions?

> 5: Ecological inchs?

> 6: threshold for a quality data

## Software Dependencies

| Software  | Version  | Purpose                                      |
|-----------|----------|----------------------------------------------|
| ![sratools](https://img.shields.io/badge/sratools-3.1.1-blue?logo=github) | 3.1.1    | Download SRA data and split into FASTQ files |
| ![fastp](https://img.shields.io/badge/fastp-0.24.0-green)           | 0.24.0   | Quality control of FASTQ files               |
| ![cutadapt](https://img.shields.io/badge/cutadapt-5.1-orange)       | 5.1      | Primer detection and trimming                |
| ![python3](https://img.shields.io/badge/python-3.9%2B-yellow?logo=python) | ≥3.9     | Run Python scripts                           |
| ![R](https://img.shields.io/badge/R-4.4.3-blue?logo=r)              | 4.4.3    | Run R packages                               |
| ![seqkit](https://img.shields.io/badge/seqkit-2.10.0-lightgrey)     | 2.10.0   | FASTQ file statistics                        |
| ![rush](https://img.shields.io/badge/rush-0.6.1-lightgrey)          | 0.6.1    | Parallel task execution                      |
| ![csvtk](https://img.shields.io/badge/csvtk-0.33.0-lightgrey)       | 0.33.0   | Table filtering                              |

### R Package Dependencies

| Package | Version | Purpose                           |
|---------|---------|-----------------------------------|
| ![dada2](https://img.shields.io/badge/dada2-1.34.0-blue) | 1.34.0  | Amplicon sequence analysis        |
| getopt  |         | Command-line argument parsing     |

---

## Dataset

### 1. 16S rRNA gene reference database
> To build a 16S rRNA gene reference database for verifying FASTQ data

```bash
# download archaeal 16S rRNA gene sequence
wget -c https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Archaea/archaea.16SrRNA.fna.gz

# bacterial 
wget -c https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz
gunzip *.gz

# formatting header line
cat archaea.16SrRNA.fna | awk '{print $1}' | sed 's/>/>archaea__/' > arch_bac_nr_16s_ref.fna

# makeblastdb
makeblastdb -in arch_bac_nr_16s_ref.fna -input_type fasta -db_type nucl -out arch_bac_nr_16s_ref
```

## Workflow
### 00. Obtaining all 16S rRNA amplicon data

```bash
# all accessions
wget -c https://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata/SRA_Accessions.tab

# all metadata, one submission, one folder
wget -c https://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata/NCBI_SRA_Metadata_Full_20250619.tar.gz

tar -zxvf NCBI_SRA_Metadata_Full_20250619.tar.gz

# parse each submission and get metadata for each accession
head -n 1 SRA_Accessions.tab > SRA_Accessions.tab.live.run.public

# the SRA_Accessions.tab includes all submissions such as run, expriment, bioproject, study, ...
## if the bioproject is marked as '-', meaning the bioproject is not public
sed '1d' SRA_Accessions.tab | awk -F '\t' '$3 ~ "live" && $7 ~ "RUN" && $19 !~ "-" && $9 ~ "public" {print}' >> SRA_Accessions.tab.live.run.public

# one submission might includs multiple runs or experiments, and not every submission includes all related files. sometimes some of the xml files are missing.
cat SRA_Accessions.tab.live.run.public | awk '{print $2}' | sed '1d' | sort -u > SRA_Accessions.tab.live.run.public.submission_acc.uniq

# xml files
# submission_accession.experiment.xml  ~.run.xml  ~.sample.xml  ~.study.xml	~.submission.xml
# search exparment.xml

python3 extract_experiment_info.py -o submission.accession.unique.experiment.tsv -t 24 NCBI_SRA_Metadata_Full_20250619

Rscript merge_table.R SRA_Accessions.tab.live.run.public submission.accession.unique.experiment.tsv SRA_Accessions.tab.live.run.public.add_experiment 24

cat SRA_Accessions.tab.live.run.public.add_experiment | csvtk -t filter2 -j 8 -f '$lib_strategy == "AMPLICON" && ($lib_source == "METAGENOMIC" || $lib_source == "GENOMIC" || $lib_source == "OTHER")' > SRA_Accessions.tab.live.run.public.add_experiment.amplicon

## 16S rRNA gene amplicon sequencing data, missing, mimatch possible
head -n 1 SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics > SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics.16s
cat non_16s_keywords.list | grep -f - SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics > SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics.non16s
cat non_16s_keywords.list | grep -v -f - SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics > SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics.16s
grep -i '16s' SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics.non16s >> SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics.16s
```
text file: ```SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics.16s```, *mismatch or missing possisble*

---

***For large project, better to splitting it into batchs***

**Detail steps examplified by batch2**
```text
02_batch.bioproject.list #bioproject accession list, one per line
```

### Step 1: Select SRA Records by BioProject Accession
Like fishing, bioproject accession like baits, and file ```SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics.16s``` is the lake.

```bash
# using csvtk
csvtk grep -t -f BioProject --pattern-file 02_batch.bioproject.list SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics.16s > 02_batch. SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics.16s

# or use grep
head -n 1 SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics > 02_batch.bioplicon.metagenomics # header line

cat 02_batch.bioproject.list | grep -w -f - SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics.16s >> SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics.16s
```
***if need, please add more specific filtering string, exclusively match non-16s records***

### Step2: Download SRA using prefetch
*Adjust the number of threads according to the network download speed. Too many threads will waste computing resources, while too few will underutilize the available bandwidth.*

```bash
# If the size of the amplicon data exceeds 1GB, exercise caution—it is likely not genuine amplicon data.
# $2 is SRA accession
cat 02_batch.SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics.16s | awk -F '\t' '{print $2}' | rush -j 24 --continue --eta --succ-cmd-file rush_prefetch.finished 'prefetch {} -O sra &> /dev/null'
```
***Error-prone areas:***
*Is downloading successful?*
```bash
wc -l rush_prefetch.finished # finished jobs
cat 02_batch.SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics.16s | sed '1d' | wcl # how many jobs
# if == Yes, prefetch successfully
# if != No, some failed

#  if !=, checked failed ones
awk '{print $2}' rush_prefetch.finished | grep -w -v -f - selectd.SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics.16s 

# if failed due to internet broken,
# remove the broken ones, if exist,
awk '{print $2}' rush_prefetch.finished | grep -w -v -f - selectd.SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics.16s | awk -F '\t' '{print $2}' | sed '1d' | rush -j 4 'rm -rf sra/{1}'

#  continue the failed ones
cat 02_batch.SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics.16s | awk -F '\t' '{print $2}' | rush -j 24 --continue --eta --succ-cmd-file rush_prefetch.finished 'prefetch {} -O sra &> /dev/null'
```

### Step 3: Split SRA to fastq

```bash
mkdir fq
cat 02_batch.SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics.16s | awk -F '\t' '{print $2}' | sed '1d' | rush -j 48 --continue --eta --succ-cmd-file rush_fastq_dump.finished 'fastq-dump --threads 1 --split-3 --outdir fq sra/{1}/{1}.sra &>/dev/null && rm -rf sra/{1}'

# if faild
cat 02_batch.SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics.16s | awk -F '\t' '{print $2}' | sed '1d' | rush -j 48 --continue --eta --succ-cmd-file rush_fastq_dump.finished 'fastq-dump sra/{1}/{1}.sra --threads 1 --split-3 --outdir fq'

# for saving disk, split, gzip and remove sra simultaneously
cat 02_batch.SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics.16s | awk -F '\t' '{print $2}' | sed '1d' | rush -j 48 --continue --eta --succ-cmd-file rush_fastq_dump.finished 'fastq-dump  --threads 1 --split-3 --outdir fq --gzip sra/{1}/{1}.sra && rm -rf sra/{1}'

# remove empty sra
rmdir sra
```

### Arrange fq by bioproject

```bash
# create bioproject directory name after biopject accession
cat selectd.SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics.16s | sed '1d' | awk -F '\t' '{print $19}' | sort -u | xargs -I {} mkdir -p {}/00_fq

cat selectd.SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics.16s | sed '1d' | awk -F '\t' '{print $2,$19}' | rush -j 4 --eta --verbose 'mv fq/{}(_[12])?.fastq(.gz)? {2}/00_fq/'
```

Remove SRA files and their parent directories (use with caution):

```bash
cat sra2bioproject.list | awk -F '\t' '{print $1}' | rush -j 48 'rm -rf 02_batch/sra/{}'
```

> **Note:**  
> Some projects contain both paired-end and single-end reads. Split these into separate BioProjects (e.g., `bioproject`, `bioproject_2`). If sequencing platforms differ, split again.

---

### dada2 by bioproject
**Paired-end reads:**  
```bash
dd2_pipeline.sh --input_dir 00_fq \
    --r1_suffix _1.fastq --r2_suffix _2.fastq \
    --threads 48 \
    --mode PE \
    --platform illumina
```
**Single-end reads:**  
```bash
# illumina
dd2_pipeline.sh --input_dir 00_fq \
    --r1_suffix _1.fastq \
    --threads 48 \
    --mode SE \
    --platform illumina

# Roche 454
d2_pipeline.sh --input_dir 00_fq \
    --r1_suffix _1.fastq \
    --threads 48 \
    --mode SE \
    --platform 454

# Ion torrent
d2_pipeline.sh --input_dir 00_fq \
    --r1_suffix _1.fastq \
    --threads 48 \
    --mode SE \
    --platform iontorrent
```
---
## Troubleshooting

> **Large-Scale Amplicon Data Analysis: Problems and Practical Fixes**

### 1. Is the SRA data derived from 16S amplicon sequencing?
```bash
is_16s_amplicon.sh -i in.fq -t 16 
```
```text
$ is_16s_amplicon.sh --input ERR6876596.fastq.gz --threads 12
sample_id            bac_hits   arch_hits  total_hits   total_percent   is_16S
ERR6876596.fastq.gz  1000       0          1000         100.0           YES
```



### 2. Were the paired-end reads properly merged?
When dada2 finished, to check the track.summary.tsv, to check how many reads were ```filtered```, how many reads are ```merged successfuly or failed```?, how many reads are ```chimeras``` and were removed. if the rest reads are less than half of the input reads, **change dada2 se mode using forward reads**

```bash
# scripts/amplicon_reads_lost_check.sh
amplicon_reads_lost_check.sh -i track.summary.tsv -o reads_lost_ratio.details.tsv
```

```text
$ cat reads_lost_ratio.summary.tsv
Sample Count                : 502
  nonchim reads left ≥ 50%  : 125
  merged reads ≥ 50%        : 130
  reads retained < 50%      : 372

⚠️  Suggestion: More than 25% of samples have low merged and nonchim rates. Switch to SE analysis may improve results.
```




