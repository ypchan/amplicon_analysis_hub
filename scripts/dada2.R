#!/usr/bin/env Rscript

# ------------------------------------------------------------
# DADA2 Amplicon Processing Pipeline (PE or SE mode)
# Author: yanpengch@qq.com
# Date: 2025-08-13 (last update)
# Usage:
#   dada2.R -i <input_dir> -o <output_dir> -m pe|se [options]
# Description:
#   This script processes Illumina paired-end or single-end amplicon reads
#   using the DADA2 pipeline, and generates ASV tables and taxonomy assignment.
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(getopt)
  library(dada2)
  library(Biostrings)
})

# -------------------------------
# Argument specification
# -------------------------------
print_usage <- function() {
  cat("
DADA2 Amplicon Processing Pipeline
==================================

Usage:
  dada2.R -i <input_dir> -o <output_dir> -m pe|se [options]

Required:
  -i, --input_dir       Directory of input FASTQ files
  -o, --output_dir      Directory to write all outputs
  -m, --mode            Processing mode: pe or se

Options:
  -1, --reads1_suffix   Forward‑read suffix (default: _1.fastq.gz)
  -2, --reads2_suffix   Reverse‑read suffix (default: _2.fastq.gz; PE only)
  -t, --threads         Number of threads (default: 4)
  -P, --platform        Platform: illumina, 454, or iontorrent (default: illumina)
  -c, --classifier      Path to taxonomy classifier FASTA
  -h, --help            Show this help and exit

Examples:
  # Illumina SE 4 threads
  dada2.R -i 02_cutadapt -o 03_dada2 -m se -t 4 -1 _R1.fq.gz

  # iontorrent SE
  dada2.R -i 02_cutadapt -o 03_dada2 -m pe -1 _read1.fq.gz -2 _read2.fq.gz -P iontorrent

Output (in <output_dir>):
  dada2_filtered/       Filtered FASTQ files
  seqtab.nochim.rds     Non‑chimera ASV table (RDS)
  track.summary.tsv     Read counts at each step
  taxonomy.rds/.tsv     Taxonomy assignment (if -c given)
\n")
}

spec <- matrix(c(
  'input_dir',     'i', 1, "character",  'Path to directory with input FASTQ files, required',
  'output_dir',    'o', 1, "character",  'Path to output directory, required',
  'mode',          'm', 1, "character",  'Processing mode: "pe" or "se", required',
  'reads1_suffix', '1', 1, "character",  'Suffix for forward reads (default: _1.fastq.gz)',
  'reads2_suffix', '2', 1, "character",  'Suffix for reverse reads (only work for PE, default: _2.fastq.gz)',
  'threads',       't', 1, "integer",    'Number of CPU threads to use (default: 4)',
  'platform',      'P', 1, "character",  'Sequencing platform: illumina|454|iontorrent (default: illumina)',
  'classifier',    'c', 1, "character",  'Path to classifier FASTA for taxonomy',
  'help',          'h', 0, "logical",    'Show this help and exit'
), byrow = TRUE, ncol = 5)

opt <- getopt(spec, usage=FALSE)
# Print help and exit if requested or required arguments missing
if (!is.null(opt$help) || is.null(opt$input_dir) || is.null(opt$output_dir) || is.null(opt$mode) ) {
  print_usage()
  quit(status=1)
}

# Remove possible trailing slashes from input and output directories
input_dir <- sub("/+$", "", opt$input_dir)
output_dir <- sub("/+$", "", opt$output_dir)
mode <- tolower(opt$mode)
if (!mode %in% c("se", "pe")) stop("Unsupported mode: ", mode)

# Set default values for optional parameters
threads <- ifelse(is.null(opt$threads), 4, opt$threads)
reads1_suffix <- ifelse(is.null(opt$reads1_suffix), "_1.fastq.gz", opt$reads1_suffix)
reads2_suffix <- ifelse(is.null(opt$reads2_suffix), "_2.fastq.gz", opt$reads2_suffix)

platform <- tolower(ifelse(is.null(opt$platform), "illumina", opt$platform))
if (!platform %in% c("illumina", "454", "iontorrent")) {
  stop("Unsupported platform: ", platform)
}

# define function to calculate elapsed time
elapsed_time <- function(start_time) {
  end_time <- Sys.time()
  elapsed_secs <- as.numeric(difftime(end_time, start_time, units = "secs"))
  ours   <- elapsed_secs %/% 3600
  minutes <- (elapsed_secs %% 3600) %/% 60
  seconds <- round(elapsed_secs %% 60)
  return(sprintf("%02d:%02d:%02d", ours, minutes, seconds))
}

log_step <- function(..., sep = " ") {
  ts <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  msg <- paste(..., sep = sep)
  message(ts, " ", msg)
}

# Output start time and parameters
start_time0 <- Sys.time()
cat("\n           DADA2 Amplicon Analysis\n")
cat("----------------------------------------------\n")
cat("    Input directory :", input_dir, "\n")
cat("    Output directory:", output_dir, "\n")
cat("    Mode            :", toupper(mode), "\n")
cat("    Threads         :", threads, "\n")
cat("    Seq platform    :", platform, "\n")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) dir.create(output_dir, recursive=TRUE)

# Prepare input files
fastqFs <- list.files(input_dir, pattern = paste0(reads1_suffix,"$"), full.names = TRUE)
sample_names <- sub(paste0(reads1_suffix,"$"), "", basename(fastqFs))
names(fastqFs) <- sample_names
sample_count <- length(sample_names)
cat("    Sample number   :", sample_count, "\n")
cat("----------------------------------------------\n")

if (mode == "pe") {
  fastqRs <- file.path(input_dir, paste0(sample_names, reads2_suffix))
  names(fastqRs) <- sample_names
}
filtpath <- file.path(output_dir, "dada2_filtered")
dir.create(filtpath, showWarnings=FALSE)

filtFs <- file.path(filtpath, sub('.gz$', '', basename(fastqFs)))
names(filtFs) <- sample_names
if (mode == "pe") {
  filtRs <- file.path(filtpath, sub('.gz$', '', basename(fastqRs)))
  names(filtRs) <- sample_names
}
failed_sample_lst <- file.path(output_dir, "filterAndTrim_failed_samples.tsv")

# -------------------------------
# Single-End processing
# -------------------------------

if (mode == "se") {
  log_step("Step1: FilterAndTrim")
  start_time <- Sys.time()
  filter.params <- list(maxN = 0, maxEE = 1, truncQ = 11, rm.phix = TRUE, minLen = 100,
              compress = FALSE, multithread = threads, verbose = TRUE, n = 1e6)

  if (platform == "iontorrent") {
    filter_out <- do.call(filterAndTrim, c(
      list(fwd = fastqFs, filt = filtFs),
      filter.params, list(trimLeft=15)))
  } else {
    filter_out <- do.call(filterAndTrim, c(
      list(fwd = fastqFs, filt = filtFs),
      filter.params))
  }
  run_time <- elapsed_time(start_time)
  cat("    Elapsed time", run_time, "\n")

  # in case of some sample failed filtering, remove them from the list
  if (any(filter_out[,"reads.out"] == 0)) {
    failed_fqs <- rownames(filter_out)[filter_out[,"reads.out"] == 0]
    failed_samples <- sub(paste0(reads1_suffix,"$"), "", failed_fqs)
    if (length(failed_samples) >0 ) {
      write.table(failed_samples, file=failed_sample_lst, quote=FALSE, row.names=FALSE, col.names=FALSE)
      cat("    Warning: The following samples had no reads after filtering and will be removed:", paste(failed_samples, collapse=", "))
    } else {
      failed_samples <- character(0)
    }

    filtFs <- filtFs[!names(filtFs) %in% failed_samples]
    sample_names <- sample_names[!sample_names %in% failed_samples]
    filter_out <- filter_out[!rownames(filter_out) %in% failed_fqs, ]
  }

  log_step("Step2: learnErrors")
  errF <- learnErrors(filtFs, multithread=threads, randomize=TRUE, nbases = 1e+08)
  run_time <- elapsed_time(start_time)
  cat("    lElapsed time", run_time, "\n")

  # DADA
  log_step("Step3: derepFastq and dada denoising")
  start_time <- Sys.time()
  ddFs <- vector("list", length(sample_names))
  names(ddFs) <- sample_names
  if (platform == "illumina") {
    for (i in seq_along(sample_names)) {
      sam <- sample_names[i]
      cat("----------", sam, "\n")
      start_time_sub <- Sys.time()
      derep <- derepFastq(filtFs[[sam]])
      ddFs[[sam]] <- dada(derep, err=errF, multithread=threads)
      cat("----------", elapsed_time(start_time_sub), "progress", sprintf("%5d/%5d", i, sample_count),"\n")
    }
  } else {
    for (i in seq_along(sample_names)) {
      sam <- sample_names[i]
      cat("----------", sam, "\n")
      start_time_sub <- Sys.time()
      derep <- derepFastq(filtFs[[sam]])
      ddFs[[sam]] <- dada(derep, err=errF, multithread=threads, HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32)
      cat("----------", elapsed_time(start_time_sub), "progress", sprintf("%5d/%5d", i, sample_count),"\n")
      
    }
  }
  run_time <- elapsed_time(start_time)
  cat("    Elapsed time", run_time， "\n")

  # Chimera removal
  log_step("Step4: makeSequenceTable and removeBimeraDenovo")
  seqtab <- makeSequenceTable(ddFs)
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=threads)
  run_time <- elapsed_time(start_time)
  cat("    Elapsed time", run_time,"\n")

} else {
  # -------------------------------
  # Paired-End processing
  # -------------------------------
  log_step("Step1: FilterAndTrim")
  start_time <- Sys.time()

  if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

  filter_out <- filterAndTrim(fwd=fastqFs, filt=filtFs, rev=fastqRs, filt.rev = filtRs,
    maxEE=2, truncQ=11, maxN=0, rm.phix=TRUE,
    compress=FALSE, verbose=TRUE, multithread=threads, n=1e+06)
  run_time <- elapsed_time(start_time)
  cat("    Elapsed time",run_time,"\n")

  if (any(filter_out[,"reads.out"] == 0)) {
    failed_fqs <- rownames(filter_out)[filter_out[,"reads.out"] == 0]

    cat("                      Warning: The following fastq files had no reads after filtering and will be removed:", paste(failed_fqs, collapse=", "))
    # Remove failed samples from the lists
    suf_pat <- paste0("(", reads1_suffix, "|", opt$reads2_suffix, ")$")
    failed_samples <- unique(sub(suf_pat, "", failed_fqs))
    if (length(failed_samples) >0 ) {
      write.table(failed_samples, file=failed_sample_lst, quote=FALSE, row.names=FALSE, col.names=FALSE)
      cat("      Warning: The following samples had no reads after filtering and will be removed:", paste(failed_samples, collapse=", "))
    } else {
      failed_samples <- character(0)
    }
    filtFs <- filtFs[!names(filtFs) %in% failed_samples]
    filtRs <- filtRs[!names(filtRs) %in% failed_samples]
    failed_sample_fqs <- c(paste0(failed_samples, reads1_suffix),paste0(failed_samples, reads1_suffix))
    sample_names <- sample_names[!sample_names %in% failed_samples]
    filter_out <- filter_out[!rownames(filter_out) %in% failed_sample_fqs, ]
  }
  log_step("Step2: learnErrors")
  # Learn forward error rates
  errF <- learnErrors(filtFs, multithread=threads, randomize=TRUE, nbases = 1e+08)
  # Learn reverse error rates
  errR <- learnErrors(filtRs, multithread=threads, randomize=TRUE, nbases = 1e+08)
  # Sample inference and merger of paired-end reads
  run_time <- elapsed_time(start_time)
  cat("    Elapsed time",run_time,"\n")

  log_step("Step3: derepFastq, dada and mergePairs")
  # DADA for forward reads
  start_time <- Sys.time()
  mergers <- vector("list", length(sample_names))
  names(mergers) <- sample_names
  denoisedF_counts <- numeric(length(sample_names))
  denoisedR_counts <- numeric(length(sample_names))
  names(denoisedF_counts) <- sample_names
  names(denoisedR_counts) <- sample_names

  for (i in seq_along(sample_names)) {
    sam <- sample_names[i]
    cat("----------", sam, "\n")
    start_time_sub <- Sys.time()
    derepF <- derepFastq(filtFs[[sam]])
    ddF <- dada(derepF, err=errF, multithread=threads)
    derepR <- derepFastq(filtRs[[sam]])
    ddR <- dada(derepR, err=errR, multithread=threads)
    merger <- mergePairs(ddF, derepF, ddR, derepR)
    mergers[[sam]] <- merger
    denoisedF_counts[sam] <- sum(getUniques(ddF))
    denoisedR_counts[sam] <- sum(getUniques(ddR))
    cat("----------", elapsed_time(start_time_sub), "progress", sprintf("%5d/%5d", i, sample_count),"\n")
  }
  # allow garbage collection
  rm(derepF)
  rm(derepR)

  cat("    derepFastq, dada and mergePairs finished. Elapsed time", elapsed_time(start_time),"\n")

  start_time <- Sys.time()
  log_step("Step4: makeSequenceTable and removeBimeraDenovo")
  seqtab <- makeSequenceTable(mergers)
  # Chimera removal
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=threads, verbose=TRUE)
  run_time <- elapsed_time(start_time)
  cat("    Elapsed time", run_time, "\n")
}

# -------------------------------
# Save results
# -------------------------------
getN <- function(x) sum(getUniques(x))
if (mode == "se") {
  track <- data.frame(
    input    = filter_out[,"reads.in"],
    filtered = filter_out[,"reads.out"],
    denoised = sapply(ddFs, getN),
    nonchim  = rowSums(seqtab.nochim),
    row.names = sample_names
  )
} else {
  track <- data.frame(
    input     = filter_out[,"reads.in"],
    filtered  = filter_out[,"reads.out"],
    denoisedF = denoisedF_counts,
    denoisedR = denoisedR_counts,
    merged    = sapply(mergers, getN),
    nonchim   = rowSums(seqtab.nochim),
    row.names = sample_names
  )
}

saveRDS(seqtab.nochim, file=file.path(output_dir,"seqtab.nochim.rds"))
cat("    seqtab.nochim ->", file.path(output_dir,"seqtab.nochim.rds"),"\n")
write.table(track, file=file.path(output_dir,"track.summary.tsv"), sep="\t", quote=FALSE, col.names=NA)
cat("    track summary -> ", file.path(output_dir,"track.summary.tsv"),"\n")

if (!is.null(opt$classifier) && file.exists(opt$classifier)) {
  start_time <- Sys.time()
  log_step("Step5: Assign taxonomy using classifier")
  tax <- assignTaxonomy(seqtab.nochim, opt$classifier, multithread=threads, verbose=TRUE)
  # saveRDS(tax, file=file.path(output_dir,"taxonomy.rds"))
  write.table(tax, file=file.path(output_dir,"taxonomy.tsv"), sep="\t", quote=FALSE, col.names=NA)
  run_time <- elapsed_time(start_time)
  cat("    taxonomy result -> ", file.path(output_dir,"taxonomy.tsv"))
  cat("    Elapsed time",run_time," \n")
  
}

log_step('Final: cleanup', filtpath)
unlink(filtpath, recursive = TRUE, force = TRUE)
run_time <- elapsed_time(start_time0)
log_step(sprintf("DADA2 finished. Elapsed time %s", run_time))