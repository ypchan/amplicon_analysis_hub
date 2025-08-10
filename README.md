# 5M16S ![GitHub](https://img.shields.io/badge/GitHub-5M16S-blue?logo=github) [![Active Development](https://img.shields.io/badge/status-active--development-orange?style=flat-square&logo=github)](https://github.com/yourusername/5M16S)

***An end-to-end pipeline for 16S rRNA gene amplicon analysis** — from raw SRA to DADA2 ASV tables, with reproducible scripts, references, and troubleshooting.*

>🚧 **Project under active development — features and docs may change.**


![5M16S](imgs/5M16S_mainpage.png)
---

## 📑 Table of Contents
- [Overview](#overview)
- [Key Challenges](#key-challenges)
- [Software Requirements](#software-requirements)
  - [R Packages](#r-packages)
- [Reference Data](#reference-data)
- [Workflow](#workflow)
  - [00. Obtain Global 16S SRA Candidates](#00-obtain-global-16s-sra-candidates)
  - [Batching Strategy](#batching-strategy)
  - [Step 1: Select SRA Records by BioProject](#step-1-select-sra-records-by-bioproject)
  - [Step 2: Download with prefetch](#step-2-download-with-prefetch)
  - [Step 3: Convert SRA → FASTQ](#step-3-convert-sra--fastq)
  - [Arrange FASTQ by BioProject](#arrange-fastq-by-bioproject)
  - [Run DADA2 per BioProject](#run-dada2-per-bioproject)
- [Common Pitfalls](#common-pitfalls)
- [Troubleshooting](#troubleshooting)
- [Notes](#notes)

---

## Overview
ProkaAtlas is a scalable, cross-platform pipeline for 16S rRNA gene amplicon processing. It standardizes **retrieval → QC → primer trimming → DADA2** across **Illumina, Roche 454, and Ion Torrent** datasets, producing comparable ASV tables and summary reports.

---

## 🚧 Key Challenges
1. **Comprehensive Retrieval** – Capturing all available 16S amplicon datasets.
2. **Platform Diversity** – Handling Roche 454 | Illumina | Ion Torrent | PacBio.
3. **Primer Strategy** – Whether, when, and how to trim primers.
4. **Target Regions** – Managing different 16S variable regions.
5. **Ecology/Metadata** – Heterogeneous ecological metadata.
6. **Quality Thresholds** – Consistent, defensible QC criteria.

---

## 🛠 Software Requirements

| Software  | Version  | Purpose |
|-----------|----------|---------|
| ![sratools](https://img.shields.io/badge/sratools-3.1.1-blue?logo=github) | 3.1.1 | Retrieve and split SRA data into FASTQ |
| ![fastp](https://img.shields.io/badge/fastp-0.24.0-green) | 0.24.0 | FASTQ quality control and filtering |
| ![cutadapt](https://img.shields.io/badge/cutadapt-5.1-orange) | 5.1 | Primer detection and trimming |
| ![python3](https://img.shields.io/badge/python-3.9%2B-yellow?logo=python) | ≥3.9 | Run Python utilities |
| ![R](https://img.shields.io/badge/R-4.4.3-blue?logo=r) | 4.4.3 | Run R-based analysis |
| ![seqkit](https://img.shields.io/badge/seqkit-2.10.0-lightgrey) | 2.10.0 | FASTQ statistics |
| ![rush](https://img.shields.io/badge/rush-0.6.1-lightgrey) | 0.6.1 | Parallel execution |
| ![csvtk](https://img.shields.io/badge/csvtk-0.33.0-lightgrey) | 0.33.0 | Table filtering and manipulation |

### R Packages
| Package | Version | Purpose |
|---------|---------|---------|
| ![dada2](https://img.shields.io/badge/dada2-1.34.0-blue) | 1.34.0 | Amplicon sequence variant inference |
| getopt  | — | Command-line argument parsing |

---

## Reference Data
Build a 16S rRNA reference database for verifying FASTQ content (optional but recommended).

```bash
# Download archaeal & bacterial 16S rRNA reference sequences
wget -c https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Archaea/archaea.16SrRNA.fna.gz
wget -c https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz
gunzip *.gz

# Example: normalize headers and build BLAST DB
cat archaea.16SrRNA.fna | awk '{print $1}' | sed 's/>/>archaea__/' > arch_bac_nr_16s_ref.fna
makeblastdb -in arch_bac_nr_16s_ref.fna -input_type fasta -db_type nucl -out arch_bac_nr_16s_ref
```

---

## Workflow

### 00. Obtain Global 16S SRA Candidates
Fetch metadata snapshots and derive a working set of public, live RUN accessions likely to be amplicons.

```bash
# All accessions (master table)
wget -c https://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata/SRA_Accessions.tab

# Full metadata (one submission per folder)
wget -c https://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata/NCBI_SRA_Metadata_Full_20250619.tar.gz

tar -zxvf NCBI_SRA_Metadata_Full_20250619.tar.gz

# Keep public, live RUN entries with a BioProject
head -n 1 SRA_Accessions.tab > SRA_Accessions.tab.live.run.public
sed '1d' SRA_Accessions.tab | awk -F '\t' '$3 ~ "live" && $7 ~ "RUN" && $19 !~ "-" && $9 ~ "public" {print}' >> SRA_Accessions.tab.live.run.public

# One submission may contain multiple runs/experiments; not all XMLs may exist
cat SRA_Accessions.tab.live.run.public | awk '{print $2}' | sed '1d' | sort -u > SRA_Accessions.tab.live.run.public.submission_acc.uniq

# Extract experiment info from XMLs
python3 extract_experiment_info.py -o submission.accession.unique.experiment.tsv -t 24 NCBI_SRA_Metadata_Full_20250619

# Merge tables in R
Rscript merge_table.R SRA_Accessions.tab.live.run.public submission.accession.unique.experiment.tsv SRA_Accessions.tab.live.run.public.add_experiment 24

# Keep likely amplicon libraries (strategy/source)
cat SRA_Accessions.tab.live.run.public.add_experiment | \
  csvtk -t filter2 -j 8 -f '$lib_strategy == "AMPLICON" && ($lib_source == "METAGENOMIC" || $lib_source == "GENOMIC" || $lib_source == "OTHER")' \
  > SRA_Accessions.tab.live.run.public.add_experiment.amplicon

# Heuristic filter to keep 16S (allow mismatches/missing)
head -n 1 SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics > SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics.16s
cat non_16s_keywords.list | grep -f - SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics > SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics.non16s
cat non_16s_keywords.list | grep -v -f - SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics >> SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics.16s
# recover possible 16S
grep -i '16s' SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics.non16s >> SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics.16s
```

> **File:** `SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics.16s` *(mismatches/missing possible)*

### Batching Strategy
> For large projects, split into batches to improve throughput and resilience.

```
02_batch.bioproject.list   # BioProject accessions, one per line
```

### Step 1: Select SRA Records by BioProject
Treat the 16S candidate table as the lake; BioProject accessions are your baits.

```bash
# Using csvtk
csvtk grep -t -f BioProject --pattern-file 02_batch.bioproject.list \
  SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics.16s \
  > 02_batch.SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics.16s

# Or use grep
head -n 1 SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics \
  > 02_batch.bioplicon.metagenomics  # header line

cat 02_batch.bioproject.list | grep -w -f - \
  SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics.16s \
  >> 02_batch.SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics.16s
```

> **Tip:** add more specific filters to exclude non‑16S records.

### Step 2: Download with prefetch
Adjust concurrency based on your network. Too many threads can waste CPU; too few underuse bandwidth.

```bash
# If an amplicon run is >1 GB, double‑check its identity — often not true amplicon data.
# $2 = SRA run accession
cat 02_batch.SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics.16s \
  | awk -F '\t' '{print $2}' \
  | rush -j 24 --continue --eta --succ-cmd-file rush_prefetch.finished \
      'prefetch {1} -O sra &> /dev/null'
```

**Error‑prone checks:**
```bash
# Did all downloads finish?
wc -l rush_prefetch.finished                      # finished jobs
cat 02_batch.SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics.16s | sed '1d' | wc -l  # total jobs

# List failed runs
awk '{print $2}' rush_prefetch.finished \
  | grep -w -v -f - 02_batch.SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics.16s

# Cleanup broken directories
awk -F '\t' '{print $2}' 02_batch.SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics.16s | sed '1d' \
  | grep -w -v -f <(awk '{print $2}' rush_prefetch.finished) \
  | rush -j 4 'rm -rf sra/{1}'

# Retry failed
cat 02_batch.SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics.16s \
  | awk -F '\t' '{print $2}' \
  | rush -j 24 --continue --eta --succ-cmd-file rush_prefetch.finished 'prefetch {1} -O sra &> /dev/null'
```

### Step 3: Convert SRA → FASTQ

```bash
mkdir -p fq
# Split and delete SRA to save space
awk -F '\t' 'NR>1{print $2}' 02_batch.SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics.16s \
  | rush -j 48 --continue --eta --succ-cmd-file rush_fastq_dump.finished \
    'fastq-dump --threads 1 --split-3 --outdir fq sra/{1}/{1}.sra &>/dev/null && rm -rf sra/{1}'

# If any failed, try again without redirecting logs
awk -F '\t' 'NR>1{print $2}' 02_batch.SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics.16s \
  | rush -j 48 --continue --eta --succ-cmd-file rush_fastq_dump.finished \
    'fastq-dump sra/{1}/{1}.sra --threads 1 --split-3 --outdir fq'

# Or gzip on the fly, then remove SRA
awk -F '\t' 'NR>1{print $2}' 02_batch.SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics.16s \
  | rush -j 48 --continue --eta --succ-cmd-file rush_fastq_dump.finished \
    'fastq-dump --threads 1 --split-3 --outdir fq --gzip sra/{1}/{1}.sra && rm -rf sra/{1}'

# Remove empty sra dir if any
rmdir sra || true
```

### Step 4: Arrange FASTQ by BioProject
*Arrange fq files by lib_layout=PAIRD|SINGLE, platfprm=illumina|454|ion torrent, bioproject/00_fq*

```bash
fq_sorter.py --metadata 02_batch.SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics.16s --fq-dir fq --threads 4
```
***Pitfalls***
```bash
# if some fq files are not moved in fq. Check them 
ls fq | sed -E 's/(_[12])\.fastq(.gz)?//' | grep -w -f - 02_batch.SRA_Accessions.tab.live.run.public.add_experiment.amplicon.metagenomics.16s

# 1 sometimes, sra will be split into 3 files
# if the left fq files is the third fq for PE reads, remove it, like
target_project/accession1_1.fastq.gz 
target_project/accession1_1.fastq.gz 
fq/accession1.fastq.gz  # remove it

# 2 sometimes, some fq files are not moved properly, please check the metadata and move thme manually
```
### Step 5: run dd2_pipeline.sh 
dd2_pipeline.sh: seqkit -> fastp -> cutadapt -> dada2 pe| se -> check mereged reation -> if need, dada2 se -> rm 01_fastp 02_cutadapt 03_dada2/dada2_filtered

```bash
# PE
cd PAIRED/Illumina/
ls | while read project;do cd project && \
    dd2_pipeline.sh --input_dir 00_fq --r1_suffix _1.fastq --r1_suffix _2.fastq --threads 60 --mode PE --platform illumina;done

# SE
ls | while read project;do cd project && \
    dd2_pipeline.sh --input_dir 00_fq --r1_suffix _1.fastq --threads 60 --mode SE --platform illumina;done

# final check
find . -name 'track.summary.tsv' -type f 
find . -name '00_fq' -type d | xargs -I {} rm -rf {}
```

---

## Common Pitfalls
- **Shell globbing**: Always quote patterns; escape parentheses in `find` (`\(`, `\)`).
- **Mixed read types**: Do not mix PE and SE in a single BioProject directory.
- **Oversized runs**: Amplicon runs >1 GB likely misannotated — verify before processing.
- **SRA locks**: Remove stale `.sra.lock` files if `prefetch` aborts (ensure no process is running).

---

## Troubleshooting

### 1) Is this truly 16S amplicon data?
```bash
is_16s_amplicon.sh -i in.fq -t 16
```
```
$ is_16s_amplicon.sh --input ERR6876596.fastq.gz --threads 12
sample_id            bac_hits   arch_hits  total_hits   total_percent   is_16S
ERR6876596.fastq.gz  1000       0          1000         100.0           YES
```

### 2) Were paired-end reads merged well?
Check `track.summary.tsv` after DADA2: counts for `filtered`, `merged`, `nonchim`. If retained reads < 50% of input, try **SE analysis using forward reads only**.

```bash
# scripts/amplicon_reads_lost_check.sh
amplicon_reads_lost_check.sh -i track.summary.tsv -o reads_lost_ratio.details.tsv
```

```
$ cat reads_lost_ratio.summary.tsv
Sample Count                : 502
  nonchim reads left ≥ 50%  : 125
  merged reads ≥ 50%        : 130
  reads retained < 50%      : 372

⚠️  Suggestion: More than 25% of samples have low merged and nonchim rates. Switch to SE analysis may improve results.
```

---

## Notes
- **Batching** improves stability on large cohorts; reruns can target failed batches only.
- **Resource tuning**: Align `-j/--threads` with available CPU/IO bandwidth.
- **Reproducibility**: Pin software versions; export conda envs / R session info.

---

📌 *Built for scalability, reproducibility, and cross‑platform flexibility.*