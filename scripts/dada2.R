#!/usr/bin/env Rscript

# ------------------------------------------------------------
# DADA2 Amplicon Processing Pipeline (PE or SE mode)
# Author: yanpengch@qq.com
# Refactor: function-based, unified logging, PE->SE suggestion
# Date: 2025-09-07 (last update)
# Usage:
#   dada2.R -i <input_dir> -o <output_dir> -m pe|se [options]
# Description:
#   Processes Illumina paired-end or single-end amplicon reads
#   using the DADA2 pipeline, generating ASV tables and taxonomy.
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(getopt)
  library(dada2)
  library(Biostrings)
})

# -------------------------------
# Usage / CLI
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
  -1, --reads1_suffix   Forward read suffix (default: _1.fastq.gz)
  -2, --reads2_suffix   Reverse read suffix (default: _2.fastq.gz; PE only)
  -t, --threads         Number of threads (default: 4)
  -P, --platform        Platform: illumina, 454, or iontorrent (default: illumina)
  -c, --classifier      Path to taxonomy classifier FASTA
  -f, --truncLengthf    Truncate forward reads at this length (default: 0, no truncation)
  -r, --truncLengthr    Truncate reverse reads at this length (PE only; default: 0, no truncation)
  -h, --help            Show this help and exit

Output (in <output_dir>):
  dada2_filtered        Folder of filtered FASTQ files
  seqtab.nochim.rds     Non-chimera ASV table (RDS)
  track.summary.tsv     Read counts at each step
  taxonomy.tsv          Taxonomy assignment (if -c given)

Examples:
  dada2.R -i 02_cutadapt -o 03_dada2 -m pe -t 24 -1 _1.fastq.gz -2 _2.fastq.gz -P illumina
  dada2.R -i 02_cutadapt -o 03_dada2 -m se -t 24 -1 .fastq.gz -P illumina
  # for some bioprojects, error, truncLengthf/r may need to be set, reference seqkit.stat.tsv, set -f 200 -r 160
  dada2.R -i 02_cutadapt -o 03_dada2 -m pe -t 24 -1 _1.fastq.gz -2 _2.fastq.gz -P illumina -f 200 -r 160
  dada2.R -i 02_cutadapt -o 03_dada2 -m se -t 24 -1 .fastq.gz -P illumina -f 200
\n")
}

spec <- matrix(c(
  'input_dir',     'i', 1, "character",  'Input FASTQ directory (required)',
  'output_dir',    'o', 1, "character",  'Output directory (required)',
  'mode',          'm', 1, "character",  'Processing mode: "pe" or "se" (required)',
  'reads1_suffix', '1', 1, "character",  'Forward read suffix (default: _1.fastq.gz)',
  'reads2_suffix', '2', 1, "character",  'Reverse read suffix (PE only, default: _2.fastq.gz)',
  'threads',       't', 1, "integer",    'CPU threads (default: 4)',
  'platform',      'P', 1, "character",  'Sequencing platform: illumina|454|iontorrent (default: illumina)',
  'classifier',    'c', 1, "character",  'Classifier FASTA for taxonomy',
  'truncLengthf',  'f', 1, "integer",    'Truncate reads after truncLen bases. Reads shorter than this are discarded',
  'truncLengthr',  'r', 1, "integer",    'Truncate reads after truncLen bases. Reads shorter than this are discarded (PE only)',
  'help',          'h', 0, "logical",    'Show help and exit'
), byrow = TRUE, ncol = 5)

opt <- getopt(spec, usage = FALSE)
if (!is.null(opt$help) || is.null(opt$input_dir) || is.null(opt$output_dir) || is.null(opt$mode)) {
  print_usage()
  quit(status = 1)
}

# -------------------------------
# Config / Defaults
# -------------------------------
threads        <- ifelse(is.null(opt$threads), 4, opt$threads)
reads1_suffix  <- ifelse(is.null(opt$reads1_suffix), "_1.fastq.gz", opt$reads1_suffix)
reads2_suffix  <- ifelse(is.null(opt$reads2_suffix), "_2.fastq.gz", opt$reads2_suffix)
truncLengthf   <- ifelse(is.null(opt$truncLengthf), 0, opt$truncLengthf)
truncLengthr   <- ifelse(is.null(opt$truncLengthr), 0, opt$truncLengthr)
platform       <- tolower(ifelse(is.null(opt$platform), "illumina", opt$platform))
if (!platform %in% c("illumina", "454", "iontorrent")) stop("Unsupported platform: ", platform)

input_dir   <- sub("/+$", "", opt$input_dir)
output_dir  <- sub("/+$", "", opt$output_dir)
mode        <- tolower(opt$mode)
if (!mode %in% c("se", "pe")) stop("Unsupported mode: ", mode)

# -------------------------------
# Utilities: time, logging, IO
# -------------------------------
ts_now <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

elapsed_time <- function(start_time) {
  end_time <- Sys.time()
  elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))
  h <- elapsed %/% 3600
  m <- (elapsed %% 3600) %/% 60
  s <- round(elapsed %% 60)
  sprintf("%02d:%02d:%02d", h, m, s)
}

log_msg <- function(level = "INFO", fmt, ...) {
  # level: INFO | WARN | ERROR | STEP
  ts <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  msg <- sprintf(fmt, ...)
  if (level %in% c("WARN", "ERROR")) {
    message(ts, " [", level, "] ", msg)
  } else {
    cat(ts, " [", level, "] ", msg, "\n", sep = "")
  }
}
log_info <- function(fmt, ...) log_msg("INFO", fmt, ...)
log_warn <- function(fmt, ...) log_msg("WARN", fmt, ...)
log_step <- function(fmt, ...) log_msg("STEP", fmt, ...)
log_error<- function(fmt, ...) log_msg("ERROR", fmt, ...)

ensure_dir <- function(d) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

write_and_echo <- function(lines, path) {
  # Write identical content to file, and echo to stdout
  ensure_dir(dirname(path))
  writeLines(lines, con = path)
  cat(paste(lines, collapse = "\n"), "\n")
  invisible(path)
}

# -------------------------------
# Core helpers
# -------------------------------
getN <- function(x) sum(getUniques(x))

build_track_se <- function(filter_out, ddFs, seqtab_nochim, sample_names) {
  data.frame(
    input     = filter_out[, "reads.in"],
    filtered  = filter_out[, "reads.out"],
    denoised  = sapply(ddFs, getN),
    nonchim   = rowSums(seqtab_nochim),
    row.names = sample_names
  )
}

build_track_pe <- function(filter_out, denoisedF_counts, denoisedR_counts, mergers, seqtab_nochim, sample_names) {
  data.frame(
    input     = filter_out[, "reads.in"],
    filtered  = filter_out[, "reads.out"],
    denoisedF = denoisedF_counts,
    denoisedR = denoisedR_counts,
    merged    = sapply(mergers, getN),
    nonchim   = rowSums(seqtab_nochim),
    row.names = sample_names
  )
}

# Suggest switching to SE if many reads fail to merge (PE only)
suggest_pe2se <- function(track,
                          out_dir = ".",
                          ratio_thr = 0.50,   # merged/input threshold
                          frac_thr  = 0.25,   # fraction threshold to trigger
                          write_note = TRUE) {
  # Only meaningful if 'merged' exists
  if (!("merged" %in% names(track))) return(FALSE)

  input  <- track[["input"]]
  merged <- track[["merged"]]
  valid  <- is.finite(input) & input > 0 & is.finite(merged)
  if (!any(valid)) return(FALSE)

  ratio   <- merged[valid] / input[valid]
  flags   <- ratio < ratio_thr
  frac    <- mean(flags, na.rm = TRUE)
  suggest <- isTRUE(frac >= frac_thr)

  if (write_note && suggest) {
    content <- c(
      "    Too many reads failed to merge; consider switching to single-end (se) mode.",
      sprintf("    Thresholds: merged/input < %.0f%% ; fraction >= %.0f%%", 100*ratio_thr, 100*frac_thr),
      sprintf("    Triggered samples: %d (%.1f%%)", sum(flags, na.rm = TRUE), 100*frac),
      sprintf("    Median merged/input: %.1f%%", 100*stats::median(ratio, na.rm = TRUE)),
      "    Action: re-run in se mode."
    )
    write_and_echo(content, file.path(out_dir, "suggestion.pe2se.note"))
  } else {
    content <- c(
      "    PE mode completed successfully; no need to switch to SE mode.",
      sprintf("    Thresholds: merged/input < %.0f%% ; fraction >= %.0f%%", 100*ratio_thr, 100*frac_thr),
      sprintf("    Triggered samples: %d (%.1f%%)", sum(flags, na.rm = TRUE), 100*frac),
      sprintf("    Median merged/input: %.1f%%", 100*stats::median(ratio, na.rm = TRUE)),
      "    Action: continue with PE results."
    )
    write_and_echo(content, file.path(out_dir, "suggestion.is_pe.note"))
  }
  return(suggest)
}

# -------------------------------
# Input discovery
# -------------------------------
discover_inputs <- function(input_dir, reads1_suffix, reads2_suffix, mode) {
  fastqFs <- list.files(input_dir, pattern = paste0(reads1_suffix, "$"), full.names = TRUE)
  sample_names <- sub(paste0(reads1_suffix, "$"), "", basename(fastqFs))
  names(fastqFs) <- sample_names
  if (mode == "pe") {
    fastqRs <- file.path(input_dir, paste0(sample_names, reads2_suffix))
    names(fastqRs) <- sample_names
  } else {
    fastqRs <- NULL
  }
  list(fastqFs = fastqFs, fastqRs = fastqRs, sample_names = sample_names)
}

# -------------------------------
# Filter & Trim
# -------------------------------
filter_and_trim_se <- function(fastqFs, filtFs, platform, threads,
                               failed_sample_lst) {
  params <- list(
    maxN = 0, maxEE = 1, truncQ = 11, rm.phix = TRUE, minLen = 100,
    compress = FALSE, multithread = threads, verbose = TRUE, n = 1e+08,
    truncLen = truncLengthf
  )
  if (platform == "iontorrent") {
    filter_out <- do.call(filterAndTrim, c(list(fwd = fastqFs, filt = filtFs), params, list(trimLeft = 15)))
  } else {
    filter_out <- do.call(filterAndTrim, c(list(fwd = fastqFs, filt = filtFs), params))
  }

  # Drop samples with zero reads after filtering
  if (any(filter_out[, "reads.out"] == 0)) {
    failed_fqs     <- rownames(filter_out)[filter_out[, "reads.out"] == 0]
    failed_samples <- sub("(.fastq|.fastq.gz)$", "", basename(failed_fqs))
    if (length(failed_samples) > 0) {
      write.table(failed_samples, file = failed_sample_lst, quote = FALSE, row.names = FALSE, col.names = FALSE)
      log_warn("Samples removed after filtering: %s", paste(failed_samples, collapse = ", "))
    }
    filtFs       <- filtFs[!names(filtFs) %in% failed_samples]
    filter_out   <- filter_out[!rownames(filter_out) %in% failed_fqs, , drop = FALSE]
  }
  list(filtFs = filtFs, filter_out = filter_out)
}

filter_and_trim_pe <- function(fastqFs, fastqRs, filtFs, filtRs, threads,
                               reads1_suffix, reads2_suffix, failed_sample_lst) {
  if (length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")
  filter_out <- filterAndTrim(
    fwd = fastqFs, filt = filtFs,
    rev = fastqRs, filt.rev = filtRs,
    maxEE = 2, truncQ = 11, maxN = 0, rm.phix = TRUE,
    compress = FALSE, verbose = TRUE, multithread = threads, n = 1e+08,
    truncLen=c(truncLengthf, truncLengthr)
  )

  # Drop samples with zero reads after filtering
  if (any(filter_out[, "reads.out"] == 0)) {
    failed_fqs <- rownames(filter_out)[filter_out[, "reads.out"] == 0]
    suf_pat    <- paste0("(", reads1_suffix, "|", reads2_suffix, ")$")
    failed_samples <- unique(sub(suf_pat, "", basename(failed_fqs)))
    if (length(failed_samples) > 0) {
      write.table(failed_samples, file = failed_sample_lst, quote = FALSE, row.names = FALSE, col.names = FALSE)
      log_warn("Samples removed after filtering: %s", paste(failed_samples, collapse = ", "))
    }
    filtFs     <- filtFs[!names(filtFs) %in% failed_samples]
    filtRs     <- filtRs[!names(filtRs) %in% failed_samples]
    keep_rows  <- !(basename(rownames(filter_out)) %in% c(paste0(failed_samples, reads1_suffix),
                                                          paste0(failed_samples, reads2_suffix)))
    filter_out <- filter_out[keep_rows, , drop = FALSE]
  }
  list(filtFs = filtFs, filtRs = filtRs, filter_out = filter_out)
}

# -------------------------------
# Error learning
# -------------------------------
learn_err_safe <- function(files, threads, nbases = 1e+08) {
  if (length(files) == 0) stop("No reads passed the filter.")
  learnErrors(files, multithread = threads, randomize = TRUE, nbases = nbases)
}

# -------------------------------
# Denoise SE
# -------------------------------
denoise_se <- function(filtFs, sample_names, errF, threads, platform) {
  ddFs <- vector("list", length(sample_names))
  names(ddFs) <- sample_names
  for (i in seq_along(sample_names)) {
    sam <- sample_names[i]
    t0  <- Sys.time()
    log_info("SE denoise: %s (%d/%d)", sam, i, length(sample_names))
    derep <- derepFastq(filtFs[[sam]], 1e+08)
    if (platform == "illumina") {
      ddFs[[sam]] <- dada(derep, err = errF, multithread = threads)
    } else {
      ddFs[[sam]] <- dada(derep, err = errF, multithread = threads,
                          HOMOPOLYMER_GAP_PENALTY = -1, BAND_SIZE = 32)
    }
    cat("    Elapsed:", elapsed_time(t0), "\n")
  }
  ddFs
}

# -------------------------------
# Denoise + Merge PE
# -------------------------------
denoise_merge_pe <- function(filtFs, filtRs, sample_names, errF, errR, threads) {
  mergers           <- vector("list", length(sample_names))
  names(mergers)    <- sample_names
  denoisedF_counts  <- numeric(length(sample_names))
  denoisedR_counts  <- numeric(length(sample_names))
  names(denoisedF_counts) <- sample_names
  names(denoisedR_counts) <- sample_names

  for (i in seq_along(sample_names)) {
    sam <- sample_names[i]
    t0  <- Sys.time()
    log_info("PE denoise/merge: %s (%d/%d)", sam, i, length(sample_names))
    derepF <- derepFastq(filtFs[[sam]], 1e+08)
    ddF    <- dada(derepF, err = errF, multithread = threads)
    derepR <- derepFastq(filtRs[[sam]], 1e+08)
    ddR    <- dada(derepR, err = errR, multithread = threads)
    mergers[[sam]] <- mergePairs(ddF, derepF, ddR, derepR, verbose=TRUE)
    denoisedF_counts[sam] <- sum(getUniques(ddF))
    denoisedR_counts[sam] <- sum(getUniques(ddR))
    cat("    Elapsed:", elapsed_time(t0), "\n")
  }
  list(mergers = mergers,
       denoisedF_counts = denoisedF_counts,
       denoisedR_counts = denoisedR_counts)
}

# -------------------------------
# Chimera removal
# -------------------------------
chimera_remove <- function(seqtab, threads, verbose = TRUE) {
  removeBimeraDenovo(seqtab, method = "consensus", multithread = threads, verbose = verbose)
}

# -------------------------------
# Banner / Params
# -------------------------------
start_time0 <- Sys.time()
cat("\n           DADA2 Amplicon Analysis\n")
cat("==============================================\n")
cat("    Input directory : ", input_dir,  "\n")
cat("    Output directory: ", output_dir, "\n")
cat("    Mode            : ", toupper(mode), "\n")
cat("    Threads         : ", threads,     "\n")
cat("    Platform        : ", platform,    "\n")

ensure_dir(output_dir)

# Discover inputs
disc <- discover_inputs(input_dir, reads1_suffix, reads2_suffix, mode)
fastqFs      <- disc$fastqFs
fastqRs      <- disc$fastqRs
sample_names <- disc$sample_names
sample_count <- length(sample_names)

cat("    Sample count    : ", sample_count, "\n")
cat("==============================================\n")

# Filtered file paths
filtpath <- file.path(output_dir, "dada2_filtered")
ensure_dir(filtpath)
filtFs <- file.path(filtpath, sub(".gz$", "", basename(fastqFs)))
names(filtFs) <- sample_names
if (mode == "pe") {
  filtRs <- file.path(filtpath, sub(".gz$", "", basename(fastqRs)))
  names(filtRs) <- sample_names
}
failed_sample_lst <- file.path(output_dir, "filterAndTrim_failed_samples.tsv")

# -------------------------------
# Main pipeline
# -------------------------------
seqtab.nochim <- NULL
track <- NULL

if (mode == "se") {
  # ---- SE pipeline ----
  log_step("1: FilterAndTrim (SE)")
  t <- Sys.time()
  ft <- filter_and_trim_se(fastqFs, filtFs, platform, threads, failed_sample_lst)
  filtFs     <- ft$filtFs
  filter_out <- ft$filter_out
  if (nrow(filter_out) == 0) stop("No reads passed the filter.")
  cat("    Elapsed time: ", elapsed_time(t), "\n")

  log_step(" 2: learnErrors (F)")
  t <- Sys.time()
  errF <- learn_err_safe(filtFs, threads)
  cat("    Elapsed time: ", elapsed_time(t), "\n")

  log_step(" 3: derep & dada (SE)")
  t <- Sys.time()
  ddFs <- denoise_se(filtFs, names(filtFs), errF, threads, platform)
  cat("    Elapsed time: ", elapsed_time(t), "\n")

  log_step(" 4: makeSequenceTable & removeBimeraDenovo")
  t <- Sys.time()
  seqtab <- makeSequenceTable(ddFs)
  seqtab.nochim <- chimera_remove(seqtab, threads)
  cat("    Elapsed time: ", elapsed_time(t), "\n")

  track <- build_track_se(filter_out, ddFs, seqtab.nochim, names(filtFs))

} else {
  # ---- PE pipeline ----
  log_step(" 1: FilterAndTrim (PE)")
  t <- Sys.time()
  ft <- filter_and_trim_pe(fastqFs, fastqRs, filtFs, filtRs, threads,
                           reads1_suffix, reads2_suffix, failed_sample_lst)
  filtFs     <- ft$filtFs
  filtRs     <- ft$filtRs
  filter_out <- ft$filter_out
  if (nrow(filter_out) == 0) stop("No reads passed the filter.")
  cat("    Elapsed time: ", elapsed_time(t), "\n")

  log_step(" 2: learnErrors (F & R)")
  t <- Sys.time()
  errF <- learn_err_safe(filtFs, threads)
  errR <- learn_err_safe(filtRs, threads)
  cat("    Elapsed time: ", elapsed_time(t), "\n")

  log_step(" 3: derep, dada & mergePairs (PE)")
  t <- Sys.time()
  dpe <- denoise_merge_pe(filtFs, filtRs, names(filtFs), errF, errR, threads)
  mergers          <- dpe$mergers
  denoisedF_counts <- dpe$denoisedF_counts
  denoisedR_counts <- dpe$denoisedR_counts
  cat("     Denoise_merge_pe elapsed time: ", elapsed_time(t), "\n")

  log_step(" 4: makeSequenceTable & removeBimeraDenovo")
  t <- Sys.time()
  seqtab <- makeSequenceTable(mergers)
  seqtab.nochim <- chimera_remove(seqtab, threads, verbose = TRUE)
  cat("    Elapsed time: ", elapsed_time(t), "\n")

  track <- build_track_pe(filter_out, denoisedF_counts, denoisedR_counts, mergers, seqtab.nochim, names(filtFs))

  # ---- Suggest SE if many reads failed to merge ----
  need_se <- suggest_pe2se(track, out_dir = output_dir,
                           ratio_thr = 0.50, frac_thr = 0.25,
                           write_note = TRUE)

  if (isTRUE(need_se)) {
    log_warn("Switching to SE analysis on filtered F reads (reuse filtFs).")
    # Re-run SE denoising and chimera removal on F only, reusing errF
    log_step("SE fallback: derep & dada (F only)")
    t <- Sys.time()
    ddFs <- denoise_se(filtFs, names(filtFs), errF, threads, platform)
    cat("    Elapsed time: ", elapsed_time(t), "\n")

    log_step("SE fallback: makeSequenceTable & removeBimeraDenovo")
    t <- Sys.time()
    seqtab <- makeSequenceTable(ddFs)
    seqtab.nochim <- chimera_remove(seqtab, threads)
    cat("    Elapsed time: ", elapsed_time(t), "\n")

    track <- build_track_se(filter_out, ddFs, seqtab.nochim, names(filtFs))
  }
}

# -------------------------------
# Save outputs
# -------------------------------
saveRDS(seqtab.nochim, file = file.path(output_dir, "seqtab.nochim.rds"))
log_info("Saved ASV table: %s", file.path(output_dir, "seqtab.nochim.rds"))

write.table(track, file = file.path(output_dir, "track.summary.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
log_info("Saved summary: %s", file.path(output_dir, "track.summary.tsv"))

# Optional taxonomy
if (!is.null(opt$classifier) && file.exists(opt$classifier)) {
  log_step(" 5: Assign taxonomy")
  t <- Sys.time()
  tax <- assignTaxonomy(seqtab.nochim, opt$classifier, multithread = threads, verbose = TRUE)
  write.table(tax, file = file.path(output_dir, "taxonomy.tsv"),
              sep = "\t", quote = FALSE, col.names = NA)
  log_info("Saved taxonomy: %s (Elapsed: %s)", file.path(output_dir, "taxonomy.tsv"), elapsed_time(t))
}

# Cleanup
log_step("Cleanup: removing filtered FASTQ files (%s)", filtpath)
unlink(filtpath, recursive = TRUE, force = TRUE)

log_info("DADA2 finished. Total elapsed time: %s", elapsed_time(start_time0))