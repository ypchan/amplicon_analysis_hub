<p align="center">
  <img src="imgs/amplicon_analysis_hub_mainpage.png" alt="amplicon_analysis_hub: 16S, ITS, short-read, and long-read amplicon analysis" width="860">
</p>

<h1 align="center">amplicon_analysis_hub</h1>

<p align="center">
  <a href="https://github.com/ypchan/amplicon_analysis_hub"><img alt="GitHub" src="https://img.shields.io/badge/GitHub-amplicon__analysis__hub-2f80ed?logo=github&logoColor=white"></a>
  <img alt="16S" src="https://img.shields.io/badge/marker-16S-2a9d8f">
  <img alt="ITS" src="https://img.shields.io/badge/marker-ITS-e9c46a">
  <img alt="Python" src="https://img.shields.io/badge/Python-%E2%89%A53.9-3776AB?logo=python&logoColor=white">
  <img alt="R" src="https://img.shields.io/badge/R-DADA2-276DC3?logo=r&logoColor=white">
</p>

`amplicon_analysis_hub` is an integrated toolkit for processing demultiplexed FASTQ files into non-chimeric ASV tables, taxonomic annotations, and abundance tables. The main workflow covers bacterial and archaeal 16S rRNA as well as fungal ITS amplicons. It supports Illumina, MGI/DNBSEQ, Element, AVITI, Ion Torrent, Roche 454, PacBio CCS, and an explicitly experimental Nanopore mode.

> [!IMPORTANT]
> The workflow selects defaults from the combination of `marker × platform × layout`, but no universal default can replace inspection of run-specific quality profiles, the actual primers, expected amplicon lengths, negative controls, and mock communities. Every run records its resolved settings in `run_parameters.tsv` and `03_dada2/effective_parameters.tsv`.

> [!CAUTION]
> The main workflow never modifies source FASTQ files. The 16S screen defaults to `report`, which neither deletes nor excludes samples. Samples are omitted from the current analysis only when `--screen-action exclude` is explicitly selected, and the original files are still retained. In contrast, `fastq_dispatcher.py --action move` and `unify_fq_suffix.py` move source files; use `--dry-run` first.

## Contents

- [Features and design principles](#features-and-design-principles)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Input conventions](#input-conventions)
- [Ecological and platform defaults](#ecological-and-platform-defaults)
- [Complete main-workflow parameters](#complete-main-workflow-parameters)
- [Workflow stages and outputs](#workflow-stages-and-outputs)
- [Running the DADA2 engine directly](#running-the-dada2-engine-directly)
- [Complete helper-script parameters](#complete-helper-script-parameters)
- [Ecological considerations for 16S and ITS](#ecological-considerations-for-16s-and-its)
- [Performance and resource tuning](#performance-and-resource-tuning)
- [Troubleshooting](#troubleshooting)
- [Testing and reproducibility](#testing-and-reproducibility)
- [Method references](#method-references)

## Features and design principles

The primary workflow is:

```text
demultiplexed FASTQ
  ├─ optional 16S content screen: report or non-destructive exclusion
  ├─ short reads: basic fastp QC; long reads skip fastp by default
  ├─ cutadapt: 5' primers and 3' read-through into the opposite primer
  ├─ seqkit: statistics for the FASTQs that actually enter DADA2
  └─ DADA2: filter → learn errors → denoise → merge PE → remove chimeras → taxonomy
```

Core design decisions:

1. Input data are read-only; every generated file is written under `--output-dir`.
2. ITS reads are not truncated to a fixed length by default, preserving genuine biological length variation.
3. PE merging depends on amplicon length and retained read length; the workflow does not apply blind fixed-length truncation.
4. DADA2 performs vectorized dereplication and denoising. Independent pooling is the default because its memory use scales more predictably with sample count.
5. A `.state/*.done` marker is written only after its stage succeeds. A completed workflow writes `amplicon_analysis_hub.finished`.
6. Historical entry points such as `dd2_pipeline.sh` and misspelled helper names remain as compatibility wrappers, but new analyses should use the canonical names.

## Installation

### Clone and validate

```bash
git clone https://github.com/ypchan/amplicon_analysis_hub.git
cd amplicon_analysis_hub

# Validate dependencies without creating command links
bash setup.sh --check-only

# Link commands into ~/.local/bin, the default prefix
bash setup.sh
export PATH="$HOME/.local/bin:$PATH"
amplicon_analysis --version
```

To use a custom command directory:

```bash
bash setup.sh --prefix "$HOME/bin"
```

### `setup.sh` parameters

| Parameter | Default | Behavior |
|---|---:|---|
| `--prefix DIR` | `~/.local/bin` | Creates or updates symbolic links for `amplicon_analysis` and the `.sh`, `.py`, and `.R` helpers in `scripts/`. |
| `--check-only` | off | Checks dependencies without creating directories or links. |
| `--build-16s-db` | off | Downloads bacterial and archaeal 16S loci from NCBI RefSeq, clusters them at 90% with `cd-hit-est`, and builds a BLAST database. Network access is required. |
| `-h, --help` | — | Shows help. |
| `-V, --version` | — | Shows the version. |

### Dependency tiers

The core workflow requires:

- Bash 4.3 or newer, Python 3.9 or newer, and R;
- the R packages `dada2` and `getopt`;
- `fastp`, Cutadapt 4.1 or newer, and `seqkit`;
- standard GNU/POSIX tools: `awk`, `sed`, `find`, `sort`, `gzip`, `getopt`, `cmp`, `cksum`, and `realpath`.

Optional features require:

- 16S screening: NCBI BLAST+;
- 16S database construction: NCBI BLAST+, `cd-hit-est`, and `wget`;
- 16S region inference: `vsearch`;
- ontology modeling: `numpy`, `pandas`, `pyarrow`, `joblib`, `scikit-learn`, and `sentence-transformers`, with optional CUDA;
- SRA download conversion: NCBI SRA Toolkit, installed separately when needed.

Install DADA2 by following the current [Bioconductor DADA2 page](https://bioconductor.org/packages/dada2/) rather than relying on a version number frozen in this README.

### Optional 16S BLAST database

```bash
THREADS=8 MEMORY_MB=32000 bash setup.sh --build-16s-db
```

`THREADS` defaults to 4. `MEMORY_MB` defaults to 0, with the meaning defined by `cd-hit-est -M`. Build provenance is written to `data/arc_bac_16s_blastDB/database_manifest.tsv`. This database is used only for content screening and does not replace an ASV taxonomy database.

## Quick start

### Illumina or MGI paired-end 16S

```bash
amplicon_analysis \
  --input-dir raw_fastq \
  --output-dir run_16s \
  --marker 16s \
  --platform illumina \
  --mode pe \
  --r1-suffix _R1.fastq.gz \
  --r2-suffix _R2.fastq.gz \
  --threads 16
```

When the bundled GTDB training FASTA is present, `--classifier auto` automatically performs 16S taxonomic assignment.

### Illumina or MGI paired-end ITS

```bash
amplicon_analysis \
  -i raw_fastq -o run_its \
  -M its -p mgi -m pe \
  -1 _R1.fastq.gz -2 _R2.fastq.gz \
  --pool pseudo \
  --classifier unite_trainset.fa.gz \
  -t 16
```

The ITS profile uses `truncLenF/R=0`, preventing systematic loss of genuine short ITS variants through fixed-length truncation.

### PacBio CCS full-length 16S

```bash
amplicon_analysis \
  -i ccs_fastq -o run_ccs \
  -M 16s -p pacbio_ccs -m se \
  -1 .fastq.gz -t 12
```

PacBio and Nanopore profiles default to `--fastp no`. Primer removal and DADA2 long-read filtering still run.

### Inspect a resolved profile without running the workflow

```bash
amplicon_analysis -M its -p pacbio_ccs -m se --print-profile
```

## Input conventions

1. FASTQs must already be demultiplexed: one SE file or one R1/R2 pair per sample.
2. Do not mix different markers, sequencing platforms, or incompatible library batches in the same input directory.
3. For PE data, a sample name is derived by removing the complete R1 suffix from the filename. Its mate must be named `sample name + R2 suffix`.
4. Suffixes are literal strings, not regular expressions. Filenames may contain spaces but must not contain newline characters.
5. FASTQs may be gzip-compressed; `.fastq.gz` is recommended for consistency.
6. PE inputs must remain synchronized. Any upstream program that independently removes R1 and R2 reads can corrupt pairing.

Example:

```text
raw_fastq/
├── lake_A_R1.fastq.gz
├── lake_A_R2.fastq.gz
├── soil_B_R1.fastq.gz
└── soil_B_R2.fastq.gz
```

The matching options are `-m pe -1 _R1.fastq.gz -2 _R2.fastq.gz`.

`--mode auto` works as follows: when `-1` is omitted, it first searches for `*_1.fastq.gz`. If every R1 has a matching `_2.fastq.gz`, it selects PE; otherwise it selects SE. Explicitly set the mode and suffixes for non-standard names to avoid ambiguity.

## Ecological and platform defaults

### Exact DADA2 profiles

`maxLen=0` disables the maximum-length limit. `maxEE=Inf` disables expected-error filtering. For SE data, only the F value in an F/R pair is used.

| Platform | Layout | Marker | maxEE F/R | truncQ | minQ | minLen | maxLen | trimLeft | Error function | DADA alignment |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|
| Illumina/MGI/Element/AVITI | PE/SE | 16S | 2 / 2 | 2 | 0 | 100 | 0 | 0 | `loessErrfun` | `BAND_SIZE=16` |
| Illumina/MGI/Element/AVITI | PE/SE | ITS/other | 2 / 2 | 2 | 0 | 50 | 0 | 0 | `loessErrfun` | `BAND_SIZE=16` |
| Ion Torrent | SE | 16S | 2 | 2 | 0 | 100 | 0 | 15 | `loessErrfun` | band 32, homopolymer penalty -1 |
| Ion Torrent | SE | ITS/other | 2 | 2 | 0 | 50 | 0 | 15 | `loessErrfun` | band 32, homopolymer penalty -1 |
| Roche 454 | SE | 16S | 2 | 2 | 0 | 100 | 0 | 0 | `loessErrfun` | band 32, homopolymer penalty -1 |
| Roche 454 | SE | ITS/other | 2 | 2 | 0 | 50 | 0 | 0 | `loessErrfun` | band 32, homopolymer penalty -1 |
| PacBio CCS | SE | 16S | 3 | 0 | 3 | 1000 | 1800 | 0 | `PacBioErrfun` | band 32, self-consistent |
| PacBio CCS | SE | ITS/other | 5 | 0 | 3 | 100 | 3000 | 0 | `PacBioErrfun` | band 32, self-consistent |
| Nanopore, experimental | SE | 16S | Inf | 0 | 0 | 1000 | 1800 | 0 | `noqualErrfun` | band 32, homopolymer -1, self-consistent |
| Nanopore, experimental | SE | ITS/other | Inf | 0 | 0 | 100 | 3000 | 0 | `noqualErrfun` | band 32, homopolymer -1, self-consistent |

Defaults shared by all profiles:

| Parameter | Default | Meaning |
|---|---:|---|
| `truncLenF`, `truncLenR` | 0 / 0 | No fixed-length truncation. `truncQ` and `maxEE` handle low-quality tails. |
| `learn_nbases` | 100,000,000 | Total bases used to learn each error model. |
| `seed` | 100 | Random seed for `learnErrors(randomize=TRUE)`. |
| `pool` | `independent` | Samples are denoised independently, with compute scaling approximately linearly by sample. |
| `minOverlap` | 12 | Absolute minimum PE overlap. Library design and truncation should preferably leave at least 20 bp plus a margin for biological length variation. |
| `maxMismatch` | 0 | No mismatches are allowed in the PE overlap. |
| `chimera` | `consensus` | Chimeras are detected per sample and resolved by consensus. |
| `rm.phix` | on for short reads; off for long reads | Removes PhiX-like reads. |

DADA2 recommends selecting a chemistry-specific `maxLen` for Roche 454, so this project does not invent one universal upper-length boundary. Inspect the observed length distribution and pass `--max-len` explicitly.

Nanopore is not a formally accuracy-guaranteed target of this workflow. `noqualErrfun` only prevents unreliable or incomparable quality scores from dominating the error model. A mock community, negative controls, and an independent method are required for validation. Prefer a dedicated, validated Nanopore workflow when species- or strain-level full-length 16S quantification is the objective.

### fastp defaults

| Parameter | 16S short reads | ITS/other short reads | PacBio/Nanopore |
|---|---:|---:|---:|
| Run fastp | yes | yes | no |
| `length_required` | 100 | 50 | — |
| `qualified_quality_phred` | 15 | 15 | — |
| `unqualified_percent_limit` | 40% | 40% | — |
| `n_base_limit` | 0 | 0 | — |
| Adapter trimming | disabled | disabled | — |

fastp performs basic read-level QC. Cutadapt handles adapters and primers, while DADA2 `maxEE` performs the final stringent filter. This division avoids contradictory or redundant trimming by three different programs.

### Cutadapt defaults

| Parameter | Default | Explanation |
|---|---:|---|
| `--primer-mode` | `trim` | Enables primer removal. The `other` marker has no default primers and therefore automatically changes to `none`. |
| Error rate | 0.10 | Maximum error proportion in a primer match. |
| Minimum overlap | 10 | Requires at least 10 nt of primer overlap, reducing short random matches in large primer collections. |
| `--times` | 2 | Allows up to two trimming rounds to handle a 5' primer and read-through into the opposite primer. |
| `--revcomp` | enabled | Searches the reverse-complement orientation. |
| Minimum length | same as fastp | 100 for 16S; 50 for ITS/other. |
| Discard untrimmed | disabled | Avoids ecological or study-source selection bias when the broad primer list does not cover every library. |
| Cores per sample | 1 | Total `--threads` controls sample-level concurrency. |

For PE data, R1 searches for an anchored forward 5' primer and reverse-primer reverse-complement read-through. R2 uses the symmetric configuration. SE reads search for forward and reverse primers plus their read-through sequences. If a sequencing center has already removed primers reliably, use `--skip-cutadapt` (equivalent to `--primer-mode none`).

The default `data/16s_primer.tsv` and `data/its_primer.tsv` files are broad candidate collections intended for heterogeneous public data. When the experimental primer pair is known, use a study-specific TSV containing only that pair to reduce non-specific trimming.

## Complete main-workflow parameters

Command: `amplicon_analysis`. Hyphenated long options are recommended. Selected underscore aliases remain available for compatibility with historical commands.

### Input, profile, and resources

| Parameter | Default | Detailed behavior |
|---|---:|---|
| `-i, --input-dir DIR` | required | Directory of demultiplexed FASTQs. The main workflow never writes to or deletes from it. |
| `-o, --output-dir DIR` | `amplicon_analysis_results` | Root for logs, intermediates, state, and results. It may contain the read-only input directory (for example, output `.` with input `00_fq`), but it cannot equal the input directory or be located inside it. |
| `-M, --marker` | `16s` | `16s`, `its`, or `other`; controls lengths, primers, screening, and the automatic classifier policy. |
| `-p, --platform` | `illumina` | Supported platforms are listed above. `bgi/bgiseq/mgiseq/dnbseq` normalize to `mgi`; `pacbio/ccs/hifi` normalize to `pacbio_ccs`; `ont` normalizes to `nanopore`. |
| `-m, --mode` | `auto` | `auto`, `pe`, or `se`. Long-read, Ion Torrent, and 454 modes accept SE only. |
| `-1, --r1-suffix` | automatic | Explicit PE defaults to `_1.fastq.gz`; explicit SE defaults to `.fastq.gz`. |
| `-2, --r2-suffix` | `_2.fastq.gz` | PE mate suffix. |
| `-t, --threads` | 4 | Total CPU budget: concurrent fastp/Cutadapt samples and the thread count passed to seqkit/DADA2. Must be at least 1. |
| `--print-profile` | off | Prints the final DADA2 profile and exits without input data or the DADA2 package. Auto mode is displayed as PE. R and `getopt` are still required. |

### Primers, screening, and fastp

| Parameter | Default | Detailed behavior |
|---|---:|---|
| `--primer-file FILE` | `auto` | 16S uses `data/16s_primer.tsv`; ITS uses `data/its_primer.tsv`; other uses none. |
| `--primer-mode` | `trim` | `trim` or `none`. `--skip-cutadapt` is a flag alias for `--primer-mode none`. |
| `--discard-untrimmed` | off | Cutadapt discards reads or pairs without a configured primer match. This can alter community composition and should be used only when primer identity and orientation are known. |
| `--cutadapt-error NUM` | 0.10 | Passed to Cutadapt as `-e`. |
| `--cutadapt-overlap INT` | 10 | Passed to Cutadapt as `-O`. |
| `--screen` | `auto` | Runs when marker=16s and the bundled database exists; otherwise skips. May also be set to `yes` or `no`. |
| `--screen-action` | `report` | `report` writes `16s_screen.tsv`; `exclude` omits `NO` samples from this run. Source FASTQs are never deleted. |
| `--blast-db PREFIX` | bundled prefix | Nucleotide 16S BLAST database prefix without `.nhr`. Single-volume and multi-volume databases are detected. |
| `--fastp` | `auto` | Resolves to `yes` for short-read platforms and `no` for PacBio/Nanopore. |
| `--fastp-min-length` | automatic | 100 for 16S and 50 for ITS/other. It is also used as Cutadapt's minimum length. |
| `--fastp-qualified INT` | 15 | Phred threshold that defines a qualified base. |
| `--fastp-unqualified NUM` | 40 | Maximum percentage of bases below the qualified-base threshold. |

### DADA2 overrides

These parameters override only the matching value in the selected profile. Every unspecified setting retains its profile default.

| Main-workflow parameter | Passed to `dada2.R` | Default or meaning |
|---|---|---|
| `--pool` | `--pool` | `independent`, `pseudo`, or `true`. Pseudo-pooling performs approximately two independent passes; true pooling has the highest memory cost. |
| `--trunc-len-f` | `--trunc_len_f` | Default 0. Fixed-length truncation for F/SE. |
| `--trunc-len-r` | `--trunc_len_r` | Default 0. Fixed-length truncation for R. |
| `--trim-left` | `--trim_left` | Usually 0; the Ion Torrent profile uses 15. |
| `--max-ee-f` | `--max_ee_f` | Maximum expected errors for F/SE. |
| `--max-ee-r` | `--max_ee_r` | Maximum expected errors for R. |
| `--trunc-q` | `--trunc_q` | Truncates at the first base with quality less than or equal to Q. |
| `--min-q` | `--min_q` | Rejects a read containing any base below this quality. |
| `--min-len` | `--min_len` | Minimum retained length after DADA2 filtering. |
| `--max-len` | `--max_len` | Maximum retained length; 0 disables the limit. |
| `--learn-nbases` | `--learn_nbases` | Default 1e8. It may be lowered for small pilots or raised for large heterogeneous batches. |
| `--min-overlap` | `--min_overlap` | Minimum PE overlap, default 12 bp. |
| `--max-mismatch` | `--max_mismatch` | Maximum mismatches in the overlap, default 0. |
| `--chimera` | `--chimera` | `consensus`, `pooled`, `per-sample`, or `none`. Use `none` only for diagnostics. |
| `--keep-filtered` | `--keep_filtered` | Retains `03_dada2/dada2_filtered/`. |

### Taxonomy, cleanup, and reruns

| Parameter | Default | Detailed behavior |
|---|---:|---|
| `-c, --classifier` | `auto` | Uses the bundled GTDB `.fna.gz` or `.fna` when available for 16S. ITS/other resolves to none. May be set explicitly to `none` or a training FASTA. |
| `--min-boot INT` | 50 | `assignTaxonomy` bootstrap cutoff in the range 0–100. Lower values increase deeper-rank assignments but also increase uncertain assignments. |
| `--cleanup` | `intermediate` | Deletes generated fastp/Cutadapt FASTQs while retaining their JSON reports, logs, summaries, source inputs, and DADA2 results. `none` retains everything. |
| `--force` | off | Reruns stages even when `.done` or completion markers exist and replaces matching generated outputs. It never touches source input files. |
| `-h, --help` | — | Prints complete help synchronized with the code. |
| `-V, --version` | — | Prints the version. |

## Workflow stages and outputs

Default directory structure:

```text
amplicon_analysis_results/
├── .state/                         # successful stage markers
├── run_parameters.tsv              # resolved main-workflow parameters
├── command.txt                     # shell-escaped original command
├── 16s_screen.tsv                  # only when the 16S screen runs
├── 01_fastp/
│   ├── SAMPLE.fastp.json
│   └── SAMPLE.fastp.log
├── 02_cutadapt/
│   └── SAMPLE.cutadapt.log
├── cutadapt_details.tsv
├── cutadapt_summary.tsv
├── seqkit.stats.tsv                # reads that actually enter DADA2
├── dada2.log
├── 03_dada2/
│   ├── effective_parameters.tsv
│   ├── error_model_r1.rds
│   ├── error_model_r2.rds          # PE only
│   ├── seqtab.nochim.rds           # primary sample × sequence count matrix
│   ├── seqtab.nochim.tsv
│   ├── ASVs.fasta
│   ├── track.summary.tsv
│   ├── taxonomy.tsv / taxonomy.rds # only when a classifier is used
│   ├── suggestion.*.note           # PE merge-retention assessment
│   └── sessionInfo.txt
└── amplicon_analysis_hub.finished
```

Stage details:

1. The screen samples at most 1,000 R1 reads per sample. Each BLAST job uses one thread, and up to `--threads` samples run concurrently. Output order follows input order.
2. fastp uses one thread per sample and obtains the total CPU budget through sample concurrency. A failure makes the workflow exit nonzero while preserving the corresponding log.
3. Cutadapt also uses sample-level concurrency and supports IUPAC-degenerate primers, orientation checks, and read-through trimming.
4. `seqkit stats --all --tabular` summarizes the files entering DADA2 rather than the original input files.
5. DADA2 performs vectorized filtering, error learning, dereplication, and denoising. PE mode then merges pairs before chimera removal and optional taxonomy assignment.
6. `track.summary.tsv` contains `input`, `filtered`, `denoisedF`, optional `denoisedR`, optional `merged`, `nonchim`, and `retained_pct`.
7. When at least 25% of PE samples have `merged/input < 50%`, the workflow writes `suggestion.pe2se.note`. It does not silently replace PE results with an SE analysis. Inspect primers, direction, quality, amplicon length, and overlap before deciding to rerun forward reads as SE.

## Running the DADA2 engine directly

`dada2.R` accepts primer-free FASTQs and does not run fastp, Cutadapt, or the content screen. Direct commands use underscore-form long options such as `--input_dir`; short options are unchanged.

### Complete parameters

| Parameter | Default | Description |
|---|---:|---|
| `-i, --input_dir` | required | Directory containing primer-free FASTQs. |
| `-o, --output_dir` | required | DADA2 output directory. |
| `-m, --mode` | `pe` | `pe` or `se`. |
| `-M, --marker` | `16s` | `16s`, `its`, or `other`. |
| `-P, --platform` | `illumina` | Supported platforms and aliases match the main workflow. |
| `-1, --reads1_suffix` | `_1.fastq.gz` | F/SE suffix. |
| `-2, --reads2_suffix` | `_2.fastq.gz` | R suffix. |
| `-t, --threads` | 4 | DADA2 multithreading. |
| `-c, --classifier` | none | DADA2 taxonomy training FASTA. |
| `-f, --trunc_len_f` | profile | Fixed F/SE truncation. |
| `-r, --trunc_len_r` | profile | Fixed R truncation. |
| `--trim_left` | profile | Fixed leading-base removal. |
| `--max_ee_f`, `--max_ee_r` | profile | Expected-error limits. |
| `--trunc_q`, `--min_q` | profile | Quality truncation and rejection settings. |
| `--min_len`, `--max_len` | profile | Length window; maximum 0 disables the upper limit. |
| `--learn_nbases` | 1e8 | Total bases used for error learning. |
| `--pool` | `independent` | `independent`, `pseudo`, or `true`. |
| `--seed` | 100 | Error-learning random seed. |
| `--min_overlap` | 12 | Minimum PE overlap. |
| `--max_mismatch` | 0 | Maximum overlap mismatches. |
| `--chimera` | `consensus` | `consensus`, `pooled`, `per-sample`, or `none`. |
| `--min_boot` | 50 | Taxonomy bootstrap cutoff. |
| `--no_try_rc` | off | Disables reverse-complement taxonomy attempts. |
| `--keep_filtered` | off | Retains DADA2-filtered FASTQs. |
| `--print_profile` | off | Writes the resolved profile as TSV and exits. |
| `-h, --help` | — | Shows help. |
| `-V, --version` | — | Shows the version. |

Example:

```bash
dada2.R -i primer_free -o dada2_out -m pe -M its -P illumina \
  -1 _R1.fastq.gz -2 _R2.fastq.gz --pool pseudo
```

## Complete helper-script parameters

### `is_16s_amplicon.py`: 16S content screening

This command reads the first N records directly from plain or gzip FASTQs instead of using a `seqkit head | fq2fa` subprocess chain. The percentage denominator is the number of reads actually sampled, not the requested maximum.

| Parameter | Default | Description |
|---|---:|---|
| `inputs` | required | One or more FASTQs. `-` reads paths from standard input and may be mixed with explicit paths. Duplicate paths are removed. |
| `-d, --db` | bundled DB | BLAST database prefix. |
| `-n, --nreads` | 100 | Maximum sampled reads per FASTQ. The main workflow overrides this to 1,000. |
| `-p, --identity` | 70% | Minimum nucleotide identity for a passing hit. |
| `-q, --query-coverage` | 70% | Minimum query coverage for a passing hit. |
| `--hit-threshold` | 50% | Percentage of sampled reads requiring a hit for a `YES` verdict. |
| `-t, --threads` | 1 | BLAST threads per FASTQ. Approximate total CPUs are threads × concurrent jobs. |
| `-c, --concurrent` | 1 | Number of FASTQs processed concurrently. |
| `--format` | `table` | Standard-output format: `table`, `tsv`, or `csv`. |
| `-o, --output` | none | Optional simultaneous output file. |
| `--out-format` | inferred from suffix | `tsv` or `csv`. |
| `-h, --help` | — | Shows help. |
| `--version` | — | Shows the version. |

Output columns are `sample_id, bac_hits, arch_hits, total_hits, total_percent, is_16S`. This is a conservative content screen, not a taxonomic analysis. Manually review low-diversity, archaeal-rich, short-fragment, or highly divergent environmental samples before excluding them.

### `fastq_dispatcher.py`: dispatch by project, layout, and platform

```bash
fastq_dispatcher.py -m metadata.tsv -f fastq -o PROCESSING --header \
  --run-col 1 --bioproject-col 22 --layout-col 16 --platform-col 19 \
  --action move --dry-run
```

| Parameter | Default | Description |
|---|---:|---|
| `-m, --metadata TSV` | required | Tab-delimited metadata. |
| `-f, --fq-dir DIR` | required | Flat FASTQ source directory. Recognizes `ACC[_1/_2].fastq[.gz]` and `.fq[.gz]`. |
| `-o, --out-dir` | `.` | Output root; `--out-root` is retained as an alias. |
| `-t, --threads` | min(16, CPU) | Concurrent file operations. |
| `--header` | off | Skips the first metadata row; `--header_1strow` is retained as an alias. |
| `--run-col` | 1 | One-based run-accession column. |
| `--bioproject-col` | 22 | One-based BioProject column. |
| `--layout-col` | 16 | Layout column; 0 disables it. An observed file pair takes precedence over metadata. |
| `--platform-col` | 19 | Platform column; 0 disables it. |
| `--action` | `move` | `move`, `copy`, or `symlink`. Move changes the source; copy doubles storage; symlink requires the source path to remain available. |
| `--dry-run` | off | Performs complete validation and reports actions without creating, moving, copying, or linking files. |
| `-h, --help` | — | Shows help. |

Output directories follow `BioProject_layout_platform/00_fq/`. For the three-file anomaly `ACC`, `ACC_1`, and `ACC_2`, the pair is placed in `00_fq`, the SE file is placed in `sra_3_fq_discarded`, and a note is recorded. Reports include `bioproject_missing.summary.tsv`, `platform_unknown.tsv`, and a collision report.

### `get_ena_fq_url_by_sra.py`: resolve ENA links

| Parameter | Default | Description |
|---|---:|---|
| `accession_file` | required | One SRR/ERR/DRR per line; accepts `-`. Blank lines and `#` comments are ignored, and accessions are deduplicated. |
| `--out-tsv` | `links.tsv` | Writes `run_accession, kind, url, md5, bytes`. |
| `--missing` | `ENA_missing.tsv` | Accessions for which no usable link was returned. |
| `--prefer` | `fastq` | `fastq`, `sra`, `submitted`, or `any`. Falls back to any available kind when the preferred kind is absent. |
| `--scheme` | `https` | `https` or `ftp`. |
| `--batch` | 100 | Accessions per API request, from 1 through 1,000. |
| `--threads` | 8 | Concurrent API requests. Requests use up to five exponential-backoff retries. |
| `-h, --help` | — | Shows help. |
| `--version` | — | Shows the version. |

This command generates links and does not download files. ENA is an external service; excessive concurrency may be rate-limited. Concurrent responses are written in the original accession order for reproducible output.

### `summarize_cutadapt.py`: parse Cutadapt logs

| Parameter | Default | Description |
|---|---:|---|
| `-d, --dir` | required | Directory containing `*.cutadapt.log`. |
| `-m, --mode` | required | `PE` or `SE`; case is normalized automatically. |
| `-o, --output-dir` | `.` | Output directory. |
| `-t, --threads` | 4 | Concurrent log parsers. |
| `-h, --help` | — | Shows help. |
| `--version` | — | Shows the version. |

The command writes `cutadapt_details.tsv` and `cutadapt_summary.tsv` in stable sample order. It returns nonzero when a report lacks its total-read field, preventing an incomplete report from silently becoming an apparently valid summary.

### `unify_fq_suffix.py`: normalize names and compression

| Parameter | Default | Description |
|---|---:|---|
| `-i, --input-dir` | required | Flat FASTQ directory. |
| `-1, --reads1-suffix-in` | required | Current R1/SE suffix. |
| `-2, --reads2-suffix-in` | none | Current R2 suffix. When present, every mate is validated before any file changes. |
| `-f, --reads1-suffix-out` | `_R1.fastq.gz` | New R1/SE suffix. |
| `-r, --reads2-suffix-out` | `_R2.fastq.gz` | New R2 suffix. |
| `-t, --threads` | 4 | Concurrent transformations. |
| `-s, --trim` | empty | Comma-separated literal tokens removed from sample names. |
| `--dry-run` | off | Runs the complete preflight and prints mappings without changing files. |
| `-h, --help` | — | Shows help. |
| `--version` | — | Shows the version. |

The command compresses when a target ends in `.gz` and the source does not, and decompresses when a source ends in `.gz` and the target does not. It uses temporary files and atomic renames. Destination conflicts, missing mates, and many-to-one names are rejected before changes begin. Successful transformations move or replace source files.

### `split_fq12.sh`: split mixed FASTQ records

| Parameter | Default | Description |
|---|---:|---|
| `-i, --input` | required | Plain or gzip-compressed mixed FASTQ. |
| `-1, --r1-out` | `reads_R1.fastq.gz` | R1 output. |
| `-2, --r2-out` | `reads_R2.fastq.gz` | R2 output. |
| `--delimiter` | `[.]` | awk split regular expression applied to the first token of the read ID. |
| `--field` | 3 | One-based field containing the mate label. |
| `--r1-value` | `1` | Value identifying R1. |
| `--r2-value` | `2` | Value identifying R2. |
| `--force` | off | Replaces existing outputs. |
| `-h, --help` | — | Shows help. |

The command validates four-line FASTQ structure, `@` and `+` markers, and equal sequence/quality length. Records with unknown labels are counted and skipped. Compression is selected independently from each output suffix.

### `amplicon_reads_lost_check.sh`: read retention

| Parameter | Default | Description |
|---|---:|---|
| `-i, --input` | required | DADA2 `track.summary.tsv`. Columns are located by name rather than position. |
| `-o, --output` | `reads_lost_ratio.tsv` beside input | Per-sample detail table. |
| `--retained-threshold` | 50% | Low merged/non-chimeric retention threshold. |
| `--sample-fraction` | 25% | Percentage of low-retention PE samples required to recommend SE review. |
| `-h, --help` | — | Shows help. |

Both PE and SE tables are accepted. SE has no merged value, so the output reports `NA` rather than inventing a merge rate. Reruns remove obsolete suggestion notes so contradictory recommendations cannot coexist.

### `asv_annotator.R`: rerun taxonomy and combine counts

| Parameter | Default | Description |
|---|---:|---|
| `-s, --seqtab_rds` | required | DADA2 sample × sequence RDS. |
| `-t, --train_fasta` | required | DADA2 training FASTA from GTDB, SILVA, RDP, UNITE, or another compatible source. |
| `-o, --output` | required | Comma-delimited for `.csv`; TSV otherwise. |
| `-p, --species_fasta` | none | Exact-match `addSpecies` reference, most commonly used for 16S. |
| `-M, --marker` | `16s` | `16s`, `its`, or `other`; used for applicability warnings. |
| `-b, --min_boot` | 50 | Bootstrap cutoff from 0 through 100. |
| `-n, --threads` | min(8, CPU) | Taxonomy threads. |
| `-d, --tax_delim` | `;` | Delimiter for collapsed taxonomy. |
| `--no_try_rc` | off | Disables reverse-complement taxonomy attempts. |
| `-h, --help` | — | Shows help. |
| `-V, --version` | — | Shows the version. |

Each output row contains an ASV ID, sequence, collapsed taxonomy, separate rank columns, and sample counts.

### `count_abundance.R`: taxonomic-rank abundance

This canonical command replaces the historical misspelling `dd2_count_aboundance.R`.

| Parameter | Default | Description |
|---|---:|---|
| `-s, --seqtab` | required | Sample × sequence RDS. |
| `-t, --taxonomy` | required | DADA2 taxonomy TSV or another TSV with a sequence-key column. |
| `-r, --rank` | `Genus` | One of the seven standard taxonomic ranks. |
| `-k, --key` | automatic | Sequence-key column. |
| `-o, --outdir` | `.` | Output directory. |
| `--taxon_col` | none | Column containing a delimited taxonomy string. |
| `--tax_sep` | `;` | Taxonomy-string delimiter. |
| `--prefix` | `abundance` | Output filename prefix. |
| `--unclassified` | `lineage` | `lineage`, `collapse`, or `drop`. Lineage generates labels such as `Unclassified_Proteobacteria`, preventing unrelated unknowns from being merged. |
| `--write_asv` | off | Also writes sample × ASV counts and an ASV map. |
| `-h, --help` | — | Shows help. |
| `-V, --version` | — | Shows the version. |

The implementation uses matrix `rowsum` rather than first expanding a potentially enormous sample-ASV long table. It writes count tables, relative abundance, a feature-oriented rank table, and an ASV-plus-taxonomy count table. Relative abundance is a compositional proportion, not absolute abundance.

### `infer_16s_regions.R`: infer 16S variable regions

This canonical command replaces the historical misspelling `infer_amplicon_segements.R`.

| Parameter | Default | Description |
|---|---:|---|
| `-f, --ref` | required | Full-length, coordinate-compatible 16S FASTA. |
| `-l, --rds_list` | none | Text file containing one seqtab RDS path per line. |
| `-r, --rds` | none | Additional comma-separated RDS paths. At least one of `--rds_list` and `--rds` is required. |
| `-o, --out` | `16s_regions.tsv` | Output table. |
| `-n, --sample_size` | 1,000 | ASVs sampled without replacement and weighted by abundance from each table. |
| `-i, --id_threshold` | 0.80 | vsearch identity fraction. |
| `-w, --window_bp` | 30 | Half-window around median start/end coordinates for the majority cluster. |
| `--overlap_min` | 10 | Minimum base-pair overlap required to call coverage of a variable region. |
| `-v, --vsearch_bin` | `vsearch` | Executable path or command. |
| `-t, --threads` | 1 | Concurrent RDS workers. |
| `-T, --vsearch_threads` | 1 | Threads per vsearch process. Approximate total CPUs are their product. |
| `-h, --help` | — | Shows help. |
| `-V, --version` | — | Shows the version. |

Output includes sampled and aligned ASV counts, alignment rate, majority coordinate-cluster rate, median coordinates, and a V1–V9 label. Coordinates follow the E. coli convention and are approximate for divergent taxa. This command cannot infer ITS1/ITS2 boundaries.

### Ontology models

`ontology_train_cv.py` parameters:

| Parameter | Default | Description |
|---|---:|---|
| `--train` | `train.parquet` | Input containing text and labels. |
| `--ontology` | `ontology_edges.tsv` | Child-parent edge TSV. |
| `--artifacts-dir` | `artifacts` | Model and metric directory. |
| `--text-col`, `--labels-col` | `text`, `labels` | Input columns. Labels may be a list or a comma-separated string. |
| `--embedder` | `intfloat/multilingual-e5-base` | SentenceTransformer model. |
| `--batch` | 256 | Embedding batch size. |
| `--workers` | 4 | Parallel one-vs-rest estimators. |
| `--cuda` | off | Uses CUDA only when PyTorch confirms it is available. |
| `--cv` | 5 | Shuffled K-fold cross-validation folds. |
| `--target-precision` | 0.90 | Precision target used to learn a threshold per label from out-of-fold predictions. |
| `--regularization-c` | 1.0 | Inverse logistic-regression regularization. |
| `--max-iter` | 2,000 | Maximum iterations. |
| `--seed` | 42 | CV and model seed. |
| `-h, --help` | — | Shows help. |
| `--version` | — | Shows the version. |

The model is `OneVsRestClassifier(LogisticRegression)`. Threshold selection and evaluation use out-of-fold probabilities, avoiding the leakage caused by learning thresholds on the training predictions and then reporting those same predictions as cross-validation performance.

`ontology_infer.py` parameters:

| Parameter | Default | Description |
|---|---:|---|
| `--input` | `big.parquet` | A Parquet file or a flat directory of Parquet files. |
| `--artifacts-dir` | `artifacts` | Training artifacts. |
| `--text-col` | `text` | Text column. |
| `--out-dir` | `pred_parts` | Output shard directory. |
| `--batch` | 256 | Embedding batch size. |
| `--shard` | 200,000 | Parquet streaming batch and output-shard size. |
| `--cuda` | off | Uses a GPU when available. |
| `--overwrite` | off | Replaces existing `pred_part_*.parquet` files. |
| `-h, --help` | — | Shows help. |
| `--version` | — | Shows the version. |

Inference uses PyArrow batch streaming rather than accumulating an entire input file in memory. Output includes direct labels, parent closure, source file and row identifiers, plus an inference manifest.

## Ecological considerations for 16S and ITS

### 16S

- ASV sequences from different variable regions cannot be merged directly by sequence. The same organism can yield incomparable fragments. Stratify analyses by region or move to a shared taxonomic rank with explicit database and resolution limitations.
- Primers differ in their affinity for bacteria, archaea, chloroplasts, and mitochondria. The BLAST screen asks only whether reads resemble bacterial or archaeal 16S; it does not automatically solve organelle contamination.
- rRNA operon copy number varies among taxa, so read proportions are not cell proportions.
- GTDB, SILVA, and RDP use different taxonomic systems and releases. Every project should record the database name, release, training FASTA checksum, and `minBoot`.
- Short V4 fragments are useful for broad community comparisons but often provide limited species resolution. Full-length 16S adds information but does not make every genus reliably identifiable to species or strain.

### ITS

- ITS1/ITS2 may range from roughly 200–600 bp and can vary even more widely. Fixed truncation creates length-selection bias, so the default is 0.
- Short amplicons can read through from one end into the opposite primer. Removing only the 5' primer leaves reverse-complement primer sequence that can affect denoising, merging, and chimera detection.
- ITS1 and ITS2 ASVs must not be merged directly by sequence. Cross-region studies should normally compare a shared, credible taxonomic rank using the same reference database and retain region as a covariate.
- Fungal rDNA copy number, nuclei per mycelium, DNA extraction efficiency, and primer mismatches all affect reads. Relative abundance is not biomass.
- Use a UNITE DADA2 training set matching ITS1/ITS2, orientation, and the targeted fungal scope. Interpret species labels together with SH/DOI information, bootstrap support, and ecological context.

### Controls and filtering

- Include extraction blanks, PCR negatives, and a mock community. A contamination model is especially important for low-biomass samples.
- Applying one rarefaction depth to every statistical question is not a universal solution. Alpha diversity, beta diversity, differential abundance, and occupancy require explicit strategies for depth, zeros, and compositionality.
- Choose ASV prevalence and abundance filters only after inspecting negative controls, sequencing depth, and study design. The workflow does not impose an arbitrary global deletion threshold.

## Performance and resource tuning

- `--threads` is a total concurrency budget, not N threads for every sample. fastp and Cutadapt each use one thread per sample and process up to N samples concurrently; DADA2 receives N internal threads.
- DADA2 `pool=independent` is fast and memory-stable. `pseudo` usually performs roughly twice the denoising work but improves sensitivity for low-abundance ASVs shared across samples. `true` can substantially increase memory and time with total reads and unique sequences.
- `learn_nbases=1e8` is a practical starting point for many batches. Very small datasets use the available bases. Do not learn one error model across mixed sequencing runs or chemistries; analyze those runs separately.
- Because fastp and Cutadapt operate concurrently by sample, excessive threads can reduce performance on shared mechanical disks. NVMe storage tolerates higher concurrency, while network filesystems generally require conservative settings.
- `seqtab.nochim.rds` uses xz compression, which writes more slowly but saves archival space. TSV is easier to inspect but can become large for high-dimensional ASV tables.
- `--cleanup intermediate` removes only generated fastp/Cutadapt FASTQs. Use `--cleanup none` for the first run when primer or quality diagnostics may be needed.

## Troubleshooting

### Low PE merge retention

Check the following before changing modes:

1. Confirm that R1/R2 files are truly paired and have the correct orientation.
2. Confirm that primer and read-through sequences were fully removed.
3. Check whether `truncLenF + truncLenR - amplicon_length` leaves sufficient overlap.
4. For ITS, determine whether a fraction of genuine amplicons is longer than the combined read span.
5. Inspect whether reverse-read quality is too poor for reliable merging.

Only when R1 alone addresses the biological question and PE merging cannot be made reliable should it be rerun independently as `-m se -1 <R1 suffix>`. The workflow never silently replaces a PE result with SE output.

### Extensive ITS loss during filtering

- Confirm that fixed `--trunc-len-f/r` values were not set.
- Check whether `--min-len` is too high.
- Check whether the broad primer table matches the actual library.
- Inspect `cutadapt_details.tsv`, length distributions, and `track.summary.tsv`.
- Distinguish PE merge loss from DADA2 filtering loss.

### SRA Lite quality scores

SRA Lite may simplify quality scores to a constant value across an entire read, disrupting standard quality-dependent error learning. Do not treat apparently uniform Q30 values as ordinary raw Illumina quality scores. When original qualities cannot be obtained, record the source explicitly, validate the error model, and do not learn one model jointly with conventional FASTQs.

### Rerunning a completed or partial analysis

Without `--force`, an existing `amplicon_analysis_hub.finished` marker causes an immediate exit. After confirming the output directory, add `--force` to rerun. A new `--output-dir` remains the clearest way to avoid mixing states.

For a partial run, the workflow compares the existing `run_parameters.tsv` with the newly resolved settings, including a sample-name checksum and all DADA2 overrides. When they differ, it refuses to reuse `.state` and requires either a new output directory or explicit `--force`. Input and output directories also cannot contain each other, preventing generated files or cleanup operations from entering the source-data tree.

When rerunning PE as SE or disabling taxonomy in the same DADA2 output directory, obsolete PE error models, suggestion notes, and taxonomy files are removed so stale optional outputs cannot be mistaken for current results.

### Missing taxonomy files

ITS `--classifier auto` intentionally resolves to none because the repository does not bind analyses to one UNITE release. Pass a compatible DADA2 training FASTA explicitly. 16S auto-classification runs only when the bundled GTDB file actually exists.

## Testing and reproducibility

```bash
# Data-free Shell/Python/R syntax checks and CLI smoke tests
bash scripts/test_pipeline.sh

# Small helper-script integration tests
bash tests/test_helpers.sh
```

A complete biological validation suite should additionally include:

- at least one mock community with known composition;
- separate short-read PE and SE fixtures for 16S and ITS;
- genuine PacBio CCS qualities and length distributions;
- a small regression dataset for every supported platform;
- comparisons of ASV count, mock false positives, read retention, chimera rate, and taxonomic accuracy.

Archive the original command, `run_parameters.tsv`, `effective_parameters.tsv`, `sessionInfo.txt`, database checksum, primer TSV, software-environment lockfile, and input FASTQ checksums.

## Method references

- [DADA2 official Illumina tutorial](https://benjjneb.github.io/dada2/tutorial.html): `maxN=0`, expected-error filtering, error learning, merging, chimera removal, and pooling.
- [DADA2 official ITS workflow](https://benjjneb.github.io/dada2/ITS_workflow): ITS length variation, avoidance of fixed truncation, paired primer/read-through removal, and UNITE taxonomy.
- [DADA2 FAQ](https://benjjneb.github.io/dada2/faq.html): Ion Torrent `trimLeft=15`, 454/Ion Torrent homopolymer and band settings, and PE overlap considerations.
- [DADA2 official big-data workflow](https://benjjneb.github.io/dada2/bigdata.html): linearly scaling independent sample inference.
- [DADA2 pseudo-pooling documentation](https://benjjneb.github.io/dada2/pseudo.html): sensitivity and computation tradeoffs among independent, pseudo, and true pooling.
- [DADA2 Bioconductor manual](https://bioconductor.org/packages/release/bioc/manuals/dada2/man/dada2.pdf): `PacBioErrfun`, `noqualErrfun`, `BAND_SIZE`, and function arguments.
- [Cutadapt official guide](https://cutadapt.readthedocs.io/en/stable/guide.html): anchored and linked adapters, IUPAC wildcards, paired-end behavior, and read-through trimming.

These defaults are interpretable starting points, not universally optimal settings across ecosystems, primers, platforms, and sequencing batches.
