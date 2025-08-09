#!/usr/bin/env Rscript

# ------------------------------------------------------------
# DADA2 Amplicon Processing Pipeline (PE or SE mode)
# Author: yanpengch@qq.com
# Date: 2025-08-09 (fixed)
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(getopt)
  library(dada2)
  library(Biostrings)
})

print_usage <- function() {
  cat("\nDADA2 Amplicon Processing Pipeline\n",
      "==================================\n\n",
      "Usage:\n",
      "  dada2.R -i <input_dir> -o <output_dir> -m pe|se [options]\n\n",
      "Required:\n",
      "  -i, --input_dir       Directory of input FASTQ files\n",
      "  -o, --output_dir      Directory to write all outputs\n",
      "  -m, --mode            Processing mode: pe or se\n\n",
      "Options:\n",
      "  -1, --reads1_suffix   Forward-read suffix (default: _R1.fastq.gz)\n",
      "  -2, --reads2_suffix   Reverse-read suffix (default: _R2.fastq.gz; PE only)\n",
      "  -t, --threads         Number of threads (default: 4)\n",
      "  -P, --platform        Platform: illumina, 454, or iontorrent (default: illumina)\n",
      "  -c, --classifier      Path to taxonomy classifier FASTA\n",
      "  -h, --help            Show this help and exit\n\n",
      "Examples:\n",
      "  # Single-end Illumina, 4 threads\n",
      "  dada2.R -i 02_cutadapt -o 03_dada2 -m se -t 4 -1 .fastq.gz\n\n",
      "  # Single-end IonTorrent with custom suffixes\n",
      "  dada2.R -i 02_cutadapt -o 03_dada2 -m se -1 _1.fastq.gz  -P iontorrent\n\n",
      "Output (in <output_dir>):\n",
      "  dada2_filtered/       Filtered FASTQ files\n",
      "  seqtab.nochim.rds     Non-chimera ASV table (RDS)\n",
      "  track.summary.tsv     Read counts at each step\n",
      "  taxonomy.tsv          Taxonomy assignment (if -c given)\n\n", sep = "")
}

spec <- matrix(c(
  'input_dir',     'i', 1, 'character',  'Path to directory with input FASTQ files, required',
  'output_dir',    'o', 1, 'character',  'Path to output directory, required',
  'mode',          'm', 1, 'character',  'Processing mode: "pe" or "se", required',
  'reads1_suffix', '1', 1, 'character',  'Suffix for forward reads (default: _1.fastq.gz)',
  'reads2_suffix', '2', 1, 'character',  'Suffix for reverse reads (PE only; default: _2.fastq.gz)',
  'threads',       't', 1, 'integer',    'Number of CPU threads to use (default: 4)',
  'platform',      'P', 1, 'character',  'Sequencing platform: illumina|454|iontorrent (default: illumina)',
  'classifier',    'c', 1, 'character',  'Path to classifier FASTA for taxonomy',
  'help',          'h', 0, 'logical',    'Show this help and exit'
), byrow = TRUE, ncol = 5)

opt <- getopt(spec)

if (is.null(opt) || isTRUE(opt$help)) {
  print_usage();
  quit(status = 0)
}

# Validate required args
if (is.null(opt$input_dir) || is.null(opt$output_dir) || is.null(opt$mode)) {
  cat("\n[ERROR] Missing required arguments.\n\n");
  print_usage();
  quit(status = 1)
}

# Normalize & defaults
input_dir  <- sub('/+$', '', opt$input_dir)
output_dir <- sub('/+$', '', opt$output_dir)
mode       <- tolower(opt$mode)
if (!mode %in% c('se','pe')) stop('Unsupported mode: ', mode)

reads1_suffix <- opt$reads1_suffix %||% '_1.fastq.gz'
reads2_suffix <- opt$reads2_suffix %||% '_2.fastq.gz'
threads       <- as.integer(opt$threads %||% 4)
platform      <- tolower(opt$platform %||% 'illumina')
if (!platform %in% c('illumina','454','iontorrent')) stop('Unsupported platform: ', platform)
classifier    <- opt$classifier %||% NULL

# Logging helpers
elapsed_time <- function(start_time) {
  end_time     <- Sys.time()
  elapsed_secs <- as.numeric(difftime(end_time, start_time, units = 'secs'))
  hours        <- elapsed_secs %/% 3600
  minutes      <- (elapsed_secs %% 3600) %/% 60
  seconds      <- round(elapsed_secs %% 60)
  sprintf('%02d:%02d:%02d', hours, minutes, seconds)
}
log_step <- function(...) {
  ts <- format(Sys.time(), '[%m-%d %H:%M:%S]')
  message(ts, ' ', paste(..., collapse = ' '))
}

start_time0 <- Sys.time()
log_step('Starting DADA2 pipeline')
cat('    Input directory: ',  input_dir,  "\n", sep='')
cat('    Output directory: ', output_dir, "\n", sep='')

if (!is.null(classifier)) cat('    Classifier: ', classifier, "\n", sep='')

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Input files
fastqFs <- list.files(input_dir, pattern = paste0(reads1_suffix, '$'), full.names = TRUE)
if (length(fastqFs) == 0) stop('No forward FASTQ files found with suffix ', reads1_suffix, ' in ', input_dir)
sample_names <- sub(paste0(reads1_suffix, '$'), '', basename(fastqFs))
names(fastqFs) <- sample_names

cat('    Sample count   : ', length(sample_names), "\n", sep='')
cat('    Processing mode: ',  mode,       "\n", sep='')
cat('    Threads: ',          threads,    "\n", sep='')
cat('    Platform: ',         platform,   "\n", sep='')

if (mode == 'pe') {
  fastqRs <- file.path(input_dir, paste0(sample_names, reads2_suffix))
  if (!all(file.exists(fastqRs))) {
    missing <- sample_names[!file.exists(fastqRs)]
    stop('Missing reverse FASTQs for samples: ', paste(missing, collapse = ', '))
  }
  names(fastqRs) <- sample_names
}

# Output/filtered paths
filtpath <- file.path(output_dir, 'dada2_filtered')
dir.create(filtpath, showWarnings = FALSE)

filtFs <- file.path(filtpath, basename(fastqFs))
names(filtFs) <- sample_names
if (mode == 'pe') {
  filtRs <- file.path(filtpath, basename(fastqRs))
  names(filtRs) <- sample_names
}

failed_sample_lst <- file.path(output_dir, 'filterAndTrim_failed_samples.tsv')

# --------------- Filtering ---------------
start_time <- Sys.time()
if (mode == 'se') {
  filter.params <- list(maxN = 0, maxEE = 1, truncQ = 11, rm.phix = TRUE, minLen = 100,
                        compress = TRUE, multithread = threads, verbose = TRUE, n = 1e6)
  if (platform == 'iontorrent') filter.params$trimLeft <- 15
  log_step('Step1: FilterAndTrim single-end reads...')
  filter_out <- do.call(filterAndTrim, c(list(fwd = fastqFs, filt = filtFs), filter.params))
} else {
  log_step('Step1: FilterAndTrim paired-end reads...')
  filter_out <- filterAndTrim(fwd = fastqFs, filt = filtFs,
                              rev = fastqRs, filt.rev = filtRs,
                              maxEE = 2, truncQ = 11, maxN = 0, rm.phix = TRUE,
                              compress = TRUE, verbose = TRUE,
                              multithread = threads, n = 1e6)
}
log_step('Step1: FilterAndTrim completed in', elapsed_time(start_time))

# Drop failed samples (reads.out == 0)
if (any(filter_out[, 'reads.out'] == 0)) {
  failed_fqs <- rownames(filter_out)[filter_out[, 'reads.out'] == 0]
  failed_samples <- unique(sub(paste0('(', reads1_suffix, '|', reads2_suffix, ')$'), '', basename(failed_fqs)))
  if (length(failed_samples) > 0) {
    write.table(failed_samples, file = failed_sample_lst, quote = FALSE, row.names = FALSE, col.names = FALSE)
    log_step('Warning: Removing samples with zero reads after filtering:', paste(failed_samples, collapse = ', '))
    keep <- !names(filtFs) %in% failed_samples
    filtFs <- filtFs[keep]
    sample_names <- sample_names[keep]
    if (mode == 'pe') filtRs <- filtRs[keep]
    filter_out <- filter_out[!basename(rownames(filter_out)) %in% paste0(failed_samples, reads1_suffix) &
                               (!grepl(reads2_suffix, rownames(filter_out)) | !basename(rownames(filter_out)) %in% paste0(failed_samples, reads2_suffix)), ]
  }
}

# --------------- Learn errors ---------------
start_time <- Sys.time()
log_step('Step2: learnErrors')
errF <- learnErrors(filtFs, multithread = threads, randomize = TRUE)
if (mode == 'pe') errR <- learnErrors(filtRs, multithread = threads, randomize = TRUE)
log_step('Step2: learnErrors completed in', elapsed_time(start_time))

# --------------- Denoise ---------------
log_step(if (mode == 'se') 'Step3: derepFastq and dada (SE)' else 'Step3: derepFastq, dada and merge pairs (PE)')
start_time <- Sys.time()

if (mode == 'se') {
  ddFs <- vector('list', length(sample_names)); names(ddFs) <- sample_names
  for (i in seq_along(sample_names)) {
    sam <- sample_names[i]
    st  <- Sys.time()
    derep <- derepFastq(filtFs[[sam]])
    if (platform == '454') {
      ddFs[[sam]] <- dada(derep, err = errF, multithread = threads, HOMOPOLYMER_GAP_PENALTY = -1, BAND_SIZE = 32)
    } else {
      ddFs[[sam]] <- dada(derep, err = errF, multithread = threads)
    }
    log_step('   processed', sam, 'elapsed', elapsed_time(st), sprintf('progress %5d/%5d', i, length(sample_names)))
  }
  seqtab <- makeSequenceTable(ddFs)
  seqtab.nochim <- removeBimeraDenovo(seqtab, method = 'consensus', multithread = threads)
} else {
  mergers <- vector('list', length(sample_names)); names(mergers) <- sample_names
  denoisedF_counts <- numeric(length(sample_names)); names(denoisedF_counts) <- sample_names
  denoisedR_counts <- numeric(length(sample_names)); names(denoisedR_counts) <- sample_names
  for (i in seq_along(sample_names)) {
    sam <- sample_names[i]
    st  <- Sys.time()
    derepF <- derepFastq(filtFs[[sam]])
    ddF <- dada(derepF, err = errF, multithread = threads)
    derepR <- derepFastq(filtRs[[sam]])
    ddR <- dada(derepR, err = errR, multithread = threads)
    merger <- mergePairs(ddF, derepF, ddR, derepR)
    mergers[[sam]] <- merger
    denoisedF_counts[sam] <- sum(getUniques(ddF))
    denoisedR_counts[sam] <- sum(getUniques(ddR))
    log_step('   processed', sam, 'elapsed', elapsed_time(st), sprintf('progress %5d/%5d', i, length(sample_names)))
  }
  rm(derepF); rm(derepR)
  seqtab <- makeSequenceTable(mergers)
  seqtab.nochim <- removeBimeraDenovo(seqtab, method = 'consensus', multithread = threads, verbose = TRUE)
}
log_step('Step3 completed in', elapsed_time(start_time))

# --------------- Save results ---------------
getN <- function(x) sum(getUniques(x))
track_file <- file.path(output_dir, 'track.summary.tsv')
seqtab_file <- file.path(output_dir, 'seqtab.nochim.rds')

if (mode == 'se') {
  track <- data.frame(
    input    = filter_out[, 'reads.in'],
    filtered = filter_out[, 'reads.out'],
    denoised = sapply(ddFs, getN),
    nonchim  = rowSums(seqtab.nochim),
    row.names = names(filtFs)
  )
} else {
  track <- data.frame(
    input     = filter_out[, 'reads.in'],
    filtered  = filter_out[, 'reads.out'],
    denoisedF = denoisedF_counts,
    denoisedR = denoisedR_counts,
    merged    = sapply(mergers, getN),
    nonchim   = rowSums(seqtab.nochim),
    row.names = names(filtFs)
  )
}

saveRDS(seqtab.nochim, file = seqtab_file)
log_step('Results: seqtab.nochim saved to', seqtab_file)
write.table(track, file = track_file, sep = '\t', quote = FALSE, col.names = NA)
log_step('Results: track summary saved to', track_file)

# --------------- Taxonomy ---------------
if (!is.null(classifier) && file.exists(classifier)) {
  start_time <- Sys.time()
  log_step('Step5: Assign taxonomy using classifier')
  tax <- assignTaxonomy(seqtab.nochim, classifier, multithread = threads, verbose = TRUE)
  tax_file <- file.path(output_dir, 'taxonomy.tsv')
  write.table(tax, file = tax_file, sep = '\t', quote = FALSE, col.names = NA)
  log_step('Step5: Assign taxonomy completed in', elapsed_time(start_time))
  log_step('Results: taxonomy saved to', tax_file)
} else if (!is.null(classifier)) {
  log_step('WARNING: Classifier file not found at', classifier, '- skipping taxonomy assignment')
}

log_step(sprintf('DADA2 pipeline completed in %s', elapsed_time(start_time0)))
