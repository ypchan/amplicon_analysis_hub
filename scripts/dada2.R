#!/usr/bin/env Rscript

# DADA2 engine for amplicon_analysis_hub.
#
# Defaults are selected from a marker/platform profile, but every biologically
# important filter and inference setting can be overridden. Resolved settings
# are always recorded in effective_parameters.tsv.

suppressPackageStartupMessages({
  library(getopt)
})

VERSION <- "2.0.0"
SUPPORTED_MARKERS <- c("16s", "its", "other")
SUPPORTED_PLATFORMS <- c("illumina", "mgi", "element", "aviti", "iontorrent",
                         "454", "pacbio_ccs", "nanopore")

usage <- function(status = 0L) {
  cat("amplicon_analysis_hub: DADA2 ASV inference\n\n")
  cat("Usage:\n  dada2.R -i DIR -o DIR [options]\n\n")
  cat("Required:\n")
  cat("  -i, --input_dir DIR       Primer-free FASTQ directory\n")
  cat("  -o, --output_dir DIR      Output directory\n\n")
  cat("Profile selection:\n")
  cat("  -m, --mode pe|se          Read layout (default: pe)\n")
  cat("  -M, --marker NAME         16s|its|other (default: 16s)\n")
  cat("  -P, --platform NAME       illumina|mgi|element|aviti|iontorrent|454|\n")
  cat("                            pacbio_ccs|nanopore (default: illumina)\n")
  cat("      --print_profile       Print resolved defaults and exit\n\n")
  cat("Input and compute:\n")
  cat("  -1, --reads1_suffix STR   R1/SE suffix (default: _1.fastq.gz)\n")
  cat("  -2, --reads2_suffix STR   R2 suffix (default: _2.fastq.gz)\n")
  cat("  -t, --threads INT         CPU threads (default: 4)\n")
  cat("      --seed INT            Error-learning seed (default: 100)\n")
  cat("      --learn_nbases NUM    Error-learning bases (default: 100000000)\n")
  cat("      --pool MODE           independent|pseudo|true (default: independent)\n")
  cat("                            pseudo raises rare-ASV sensitivity at ~2x denoise time\n\n")
  cat("Filtering overrides (profile defaults shown by --print_profile):\n")
  cat("  -f, --trunc_len_f INT     Fixed R1/SE truncation; 0 keeps full length (default: 0)\n")
  cat("  -r, --trunc_len_r INT     Fixed R2 truncation; 0 keeps full length (default: 0)\n")
  cat("      --trim_left INT       Remove leading bases (Ion Torrent default: 15)\n")
  cat("      --max_ee_f NUM        Maximum expected errors for R1/SE\n")
  cat("      --max_ee_r NUM        Maximum expected errors for R2\n")
  cat("      --trunc_q INT         Truncate at first quality <= value\n")
  cat("      --min_q INT           Reject reads containing quality below value\n")
  cat("      --min_len INT         Minimum retained length\n")
  cat("      --max_len INT         Maximum retained length; 0 disables\n\n")
  cat("Merging, chimera and taxonomy:\n")
  cat("      --min_overlap INT     Minimum PE overlap (default: 12)\n")
  cat("      --max_mismatch INT    Maximum overlap mismatches (default: 0)\n")
  cat("      --chimera METHOD      consensus|pooled|per-sample|none (default: consensus)\n")
  cat("  -c, --classifier FASTA    DADA2 taxonomy training FASTA (optional)\n")
  cat("      --min_boot INT        Taxonomy bootstrap cutoff (default: 50)\n")
  cat("      --no_try_rc           Do not classify reverse complements\n")
  cat("      --keep_filtered       Keep generated filtered FASTQ files\n")
  cat("  -h, --help                Show this help\n")
  cat("  -V, --version             Show version\n\n")
  cat("Ecological defaults:\n")
  cat("  short-read 16S: maxEE=2/2, truncQ=2, minLen=100, no fixed truncation\n")
  cat("  short-read ITS: maxEE=2/2, truncQ=2, minLen=50, no fixed truncation\n")
  cat("  PacBio CCS 16S: maxEE=3, minQ=3, length=1000..1800, PacBioErrfun\n")
  cat("  PacBio CCS ITS: maxEE=5, minQ=3, length=100..3000, PacBioErrfun\n")
  cat("  Nanopore: experimental SE mode; quality ignored by error model\n\n")
  cat("Examples:\n")
  cat("  dada2.R -i 02_cutadapt -o 03_dada2 -m pe -M 16s -P illumina\n")
  cat("  dada2.R -i 02_cutadapt -o 03_dada2 -m pe -M its -P mgi --pool pseudo\n")
  cat("  dada2.R -i 02_cutadapt -o 03_dada2 -m se -M 16s -P pacbio_ccs -1 .fastq.gz\n")
  cat("  dada2.R -M its -P pacbio_ccs -m se --print_profile\n")
  quit(status = status)
}

spec <- matrix(c(
  "input_dir",       "i", 1, "character",
  "output_dir",      "o", 1, "character",
  "mode",            "m", 1, "character",
  "marker",          "M", 1, "character",
  "platform",        "P", 1, "character",
  "reads1_suffix",   "1", 1, "character",
  "reads2_suffix",   "2", 1, "character",
  "threads",         "t", 1, "integer",
  "classifier",      "c", 1, "character",
  "trunc_len_f",     "f", 1, "integer",
  "trunc_len_r",     "r", 1, "integer",
  "trim_left",        NA, 1, "integer",
  "max_ee_f",         NA, 1, "double",
  "max_ee_r",         NA, 1, "double",
  "trunc_q",          NA, 1, "integer",
  "min_q",            NA, 1, "integer",
  "min_len",          NA, 1, "integer",
  "max_len",          NA, 1, "integer",
  "learn_nbases",     NA, 1, "double",
  "pool",             NA, 1, "character",
  "seed",             NA, 1, "integer",
  "min_overlap",      NA, 1, "integer",
  "max_mismatch",     NA, 1, "integer",
  "chimera",          NA, 1, "character",
  "min_boot",         NA, 1, "integer",
  "no_try_rc",        NA, 0, "logical",
  "keep_filtered",    NA, 0, "logical",
  "print_profile",    NA, 0, "logical",
  "help",            "h", 0, "logical",
  "version",         "V", 0, "logical"
), byrow = TRUE, ncol = 4)

opt <- getopt(spec)
if (isTRUE(opt$help)) usage(0L)
if (isTRUE(opt$version)) {
  cat("dada2.R ", VERSION, "\n", sep = "")
  quit(status = 0L)
}

normalize_platform <- function(x) {
  x <- tolower(gsub("[- ]", "_", x))
  aliases <- c(bgi = "mgi", bgiseq = "mgi", mgiseq = "mgi", dnbseq = "mgi",
               roche454 = "454", roche_454 = "454", ion_torrent = "iontorrent",
               pacbio = "pacbio_ccs", ccs = "pacbio_ccs", hifi = "pacbio_ccs",
               ont = "nanopore", oxford_nanopore = "nanopore")
  if (x %in% names(aliases)) unname(aliases[[x]]) else x
}

profile_defaults <- function(marker, platform, mode) {
  cfg <- list(
    marker = marker, platform = platform, mode = mode,
    max_ee_f = 2, max_ee_r = 2, trunc_q = 2L, min_q = 0L,
    min_len = if (marker == "16s") 100L else 50L,
    max_len = 0L, trim_left = 0L, trunc_len_f = 0L, trunc_len_r = 0L,
    learn_nbases = 1e8, pool = "independent", seed = 100L,
    min_overlap = 12L, max_mismatch = 0L, chimera = "consensus",
    error_model = "quality", band_size = 16L,
    homopolymer_gap_penalty = NA_real_, self_consist = FALSE,
    rm_phix = platform %in% c("illumina", "mgi", "element", "aviti")
  )
  if (platform == "iontorrent") {
    cfg$trim_left <- 15L
    cfg$band_size <- 32L
    cfg$homopolymer_gap_penalty <- -1
  } else if (platform == "454") {
    cfg$band_size <- 32L
    cfg$homopolymer_gap_penalty <- -1
  } else if (platform == "pacbio_ccs") {
    cfg$max_ee_f <- if (marker == "16s") 3 else 5
    cfg$max_ee_r <- cfg$max_ee_f
    cfg$trunc_q <- 0L
    cfg$min_q <- 3L
    cfg$min_len <- if (marker == "16s") 1000L else 100L
    cfg$max_len <- if (marker == "16s") 1800L else 3000L
    cfg$error_model <- "pacbio"
    cfg$band_size <- 32L
    cfg$self_consist <- TRUE
    cfg$rm_phix <- FALSE
  } else if (platform == "nanopore") {
    cfg$max_ee_f <- Inf
    cfg$max_ee_r <- Inf
    cfg$trunc_q <- 0L
    cfg$min_len <- if (marker == "16s") 1000L else 100L
    cfg$max_len <- if (marker == "16s") 1800L else 3000L
    cfg$error_model <- "noqual"
    cfg$band_size <- 32L
    cfg$homopolymer_gap_penalty <- -1
    cfg$self_consist <- TRUE
    cfg$rm_phix <- FALSE
  }
  cfg
}

marker <- tolower(if (is.null(opt$marker)) "16s" else opt$marker)
platform <- normalize_platform(if (is.null(opt$platform)) "illumina" else opt$platform)
mode <- tolower(if (is.null(opt$mode)) "pe" else opt$mode)
if (!marker %in% SUPPORTED_MARKERS) stop("--marker must be one of: ", paste(SUPPORTED_MARKERS, collapse = ", "))
if (!platform %in% SUPPORTED_PLATFORMS) stop("--platform must be one of: ", paste(SUPPORTED_PLATFORMS, collapse = ", "))
if (!mode %in% c("pe", "se")) stop("--mode must be pe or se")
if (mode == "pe" && platform %in% c("454", "iontorrent", "pacbio_ccs", "nanopore")) {
  stop(platform, " is supported only in SE mode by this workflow")
}

cfg <- profile_defaults(marker, platform, mode)
override <- function(name) if (!is.null(opt[[name]])) cfg[[name]] <<- opt[[name]]
for (name in c("trunc_len_f", "trunc_len_r", "trim_left", "max_ee_f", "max_ee_r",
               "trunc_q", "min_q", "min_len", "max_len", "learn_nbases", "pool",
               "seed", "min_overlap", "max_mismatch", "chimera")) override(name)
cfg$pool <- tolower(as.character(cfg$pool))
cfg$chimera <- tolower(as.character(cfg$chimera))
if (!cfg$pool %in% c("independent", "pseudo", "true")) stop("--pool must be independent, pseudo, or true")
if (!cfg$chimera %in% c("consensus", "pooled", "per-sample", "none")) {
  stop("--chimera must be consensus, pooled, per-sample, or none")
}

numeric_nonnegative <- c("trunc_len_f", "trunc_len_r", "trim_left", "trunc_q", "min_q",
                         "min_len", "max_len", "min_overlap", "max_mismatch")
for (name in numeric_nonnegative) {
  if (!is.finite(cfg[[name]]) || cfg[[name]] < 0) stop("--", name, " must be non-negative")
}
if (cfg$max_len > 0 && cfg$max_len < cfg$min_len) stop("--max_len must be 0 or >= --min_len")
if (!is.finite(cfg$learn_nbases) || cfg$learn_nbases <= 0) stop("--learn_nbases must be > 0")
for (name in c("max_ee_f", "max_ee_r")) {
  if (is.na(cfg[[name]]) || cfg[[name]] <= 0) stop("--", name, " must be > 0 (Inf disables the filter)")
}
seed_integer <- suppressWarnings(as.integer(cfg$seed))
if (!is.finite(cfg$seed) || cfg$seed < 0 || is.na(seed_integer) || cfg$seed != seed_integer) {
  stop("--seed must be a non-negative integer")
}

parameter_frame <- function(config) {
  data.frame(parameter = names(config),
             value = vapply(config, function(x) paste(x, collapse = ","), character(1)),
             stringsAsFactors = FALSE)
}

if (isTRUE(opt$print_profile)) {
  write.table(parameter_frame(cfg), stdout(), sep = "\t", quote = FALSE, row.names = FALSE)
  quit(status = 0L)
}
if (is.null(opt$input_dir) || is.null(opt$output_dir)) usage(2L)
if (!requireNamespace("dada2", quietly = TRUE)) stop("R package 'dada2' is required")

threads <- if (is.null(opt$threads)) 4L else as.integer(opt$threads)
if (!is.finite(threads) || threads < 1L) stop("--threads must be >= 1")
reads1_suffix <- if (is.null(opt$reads1_suffix)) "_1.fastq.gz" else opt$reads1_suffix
reads2_suffix <- if (is.null(opt$reads2_suffix)) "_2.fastq.gz" else opt$reads2_suffix
min_boot <- if (is.null(opt$min_boot)) 50L else as.integer(opt$min_boot)
if (!is.finite(min_boot) || min_boot < 0L || min_boot > 100L) {
  stop("--min_boot must be between 0 and 100")
}
input_dir <- normalizePath(opt$input_dir, mustWork = TRUE)
output_dir <- normalizePath(opt$output_dir, mustWork = FALSE)
if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE)) stop("Cannot create: ", output_dir)

timestamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
log_message <- function(level, ...) message("[", timestamp(), "] [", level, "] ", paste0(..., collapse = ""))
elapsed <- function(start) sprintf("%.1fs", as.numeric(difftime(Sys.time(), start, units = "secs")))

discover_inputs <- function(path, suffix_f, suffix_r, layout) {
  entries <- sort(list.files(path, full.names = TRUE, recursive = FALSE))
  entry_info <- file.info(entries)
  entries <- entries[!is.na(entry_info$isdir) & !entry_info$isdir]
  forward <- entries[endsWith(entries, suffix_f)]
  if (!length(forward)) stop("No files end with --reads1_suffix '", suffix_f, "' in ", path)
  samples <- substr(basename(forward), 1L, nchar(basename(forward)) - nchar(suffix_f))
  if (any(!nzchar(samples)) || anyDuplicated(samples)) stop("FASTQ suffixes do not produce unique non-empty sample names")
  names(forward) <- samples
  reverse <- NULL
  if (layout == "pe") {
    reverse <- file.path(path, paste0(samples, suffix_r))
    names(reverse) <- samples
    missing <- reverse[!file.exists(reverse)]
    if (length(missing)) stop("Missing R2 for: ", paste(names(missing), collapse = ", "))
  }
  list(forward = forward, reverse = reverse, samples = samples)
}

inputs <- discover_inputs(input_dir, reads1_suffix, reads2_suffix, mode)
filtered_dir <- file.path(output_dir, "dada2_filtered")
if (!dir.exists(filtered_dir) && !dir.create(filtered_dir, recursive = TRUE)) stop("Cannot create: ", filtered_dir)
filt_f <- setNames(file.path(filtered_dir, basename(inputs$forward)), inputs$samples)
filt_r <- if (mode == "pe") setNames(file.path(filtered_dir, basename(inputs$reverse)), inputs$samples) else NULL

cfg$threads <- threads
cfg$reads1_suffix <- reads1_suffix
cfg$reads2_suffix <- if (mode == "pe") reads2_suffix else "NA"
cfg$min_boot <- min_boot
cfg$dada2_version <- as.character(packageVersion("dada2"))
cfg$workflow_version <- VERSION
write.table(parameter_frame(cfg), file.path(output_dir, "effective_parameters.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

if (platform == "nanopore") {
  log_message("WARN", "Nanopore support is experimental. Validate ASVs with a mock community; ",
              "the official DADA2 long-read model is PacBio CCS-specific.")
}
if (marker == "its" && (cfg$trunc_len_f > 0 || cfg$trunc_len_r > 0)) {
  log_message("WARN", "Fixed truncation can remove real ITS length variants; use only after inspecting quality profiles.")
}

filter_args <- list(maxN = 0,
                    maxEE = if (mode == "pe") c(cfg$max_ee_f, cfg$max_ee_r) else cfg$max_ee_f,
                    truncQ = cfg$trunc_q, minQ = cfg$min_q, minLen = cfg$min_len,
                    trimLeft = if (mode == "pe") c(cfg$trim_left, cfg$trim_left) else cfg$trim_left,
                    truncLen = if (mode == "pe") c(cfg$trunc_len_f, cfg$trunc_len_r) else cfg$trunc_len_f,
                    rm.phix = cfg$rm_phix, compress = TRUE,
                    multithread = threads, verbose = TRUE)
if (cfg$max_len > 0) filter_args$maxLen <- cfg$max_len

start_all <- Sys.time()
log_message("INFO", "Filtering ", length(inputs$samples), " sample(s); marker=", marker,
            ", platform=", platform, ", mode=", mode)
start <- Sys.time()
if (mode == "pe") {
  filter_out <- do.call(dada2::filterAndTrim,
                        c(list(fwd = inputs$forward, filt = filt_f,
                               rev = inputs$reverse, filt.rev = filt_r), filter_args))
} else {
  filter_out <- do.call(dada2::filterAndTrim, c(list(fwd = inputs$forward, filt = filt_f), filter_args))
}
log_message("INFO", "Filtering completed in ", elapsed(start))

keep <- filter_out[, "reads.out"] > 0
if (!any(keep)) stop("No reads passed filtering")
if (any(!keep)) {
  failed <- inputs$samples[!keep]
  writeLines(failed, file.path(output_dir, "filter_failed_samples.txt"))
  log_message("WARN", "No reads remained for: ", paste(failed, collapse = ", "))
}
filter_out <- filter_out[keep, , drop = FALSE]
filt_f <- filt_f[keep]
if (mode == "pe") filt_r <- filt_r[keep]
samples <- names(filt_f)

set.seed(as.integer(cfg$seed))
error_fun <- switch(cfg$error_model,
                    pacbio = dada2::PacBioErrfun,
                    noqual = dada2::noqualErrfun,
                    quality = dada2::loessErrfun)
learn_one <- function(files, label) {
  log_message("INFO", "Learning ", label, " error model from ",
              format(cfg$learn_nbases, scientific = FALSE), " bases")
  start <- Sys.time()
  learn_args <- list(files, nbases = cfg$learn_nbases, randomize = TRUE,
                     multithread = threads, errorEstimationFunction = error_fun,
                     verbose = TRUE)
  if (cfg$band_size != 16L) learn_args$BAND_SIZE <- cfg$band_size
  if (!is.na(cfg$homopolymer_gap_penalty)) {
    learn_args$HOMOPOLYMER_GAP_PENALTY <- cfg$homopolymer_gap_penalty
  }
  ans <- do.call(dada2::learnErrors, learn_args)
  log_message("INFO", label, " error learning completed in ", elapsed(start))
  ans
}
err_f <- learn_one(filt_f, "R1/SE")
err_r <- if (mode == "pe") learn_one(filt_r, "R2") else NULL
saveRDS(err_f, file.path(output_dir, "error_model_r1.rds"))
if (mode == "pe") {
  saveRDS(err_r, file.path(output_dir, "error_model_r2.rds"))
} else {
  unlink(file.path(output_dir, "error_model_r2.rds"), force = TRUE)
}

pool_arg <- switch(cfg$pool, independent = FALSE, pseudo = "pseudo", true = TRUE)
dada_args <- list(multithread = threads, pool = pool_arg, verbose = TRUE,
                  BAND_SIZE = cfg$band_size, selfConsist = cfg$self_consist)
if (cfg$self_consist) dada_args$errorEstimationFunction <- error_fun
if (!is.na(cfg$homopolymer_gap_penalty)) {
  dada_args$HOMOPOLYMER_GAP_PENALTY <- cfg$homopolymer_gap_penalty
}

log_message("INFO", "Dereplicating and denoising with pool=", cfg$pool)
start <- Sys.time()
derep_f <- dada2::derepFastq(filt_f, verbose = TRUE)
names(derep_f) <- samples
dd_f <- do.call(dada2::dada, c(list(derep_f, err = err_f), dada_args))
if (mode == "pe") {
  derep_r <- dada2::derepFastq(filt_r, verbose = TRUE)
  names(derep_r) <- samples
  dd_r <- do.call(dada2::dada, c(list(derep_r, err = err_r), dada_args))
  mergers <- dada2::mergePairs(dd_f, derep_f, dd_r, derep_r,
                               minOverlap = cfg$min_overlap,
                               maxMismatch = cfg$max_mismatch, verbose = TRUE)
  seqtab <- dada2::makeSequenceTable(mergers)
} else {
  dd_r <- NULL
  mergers <- NULL
  seqtab <- dada2::makeSequenceTable(dd_f)
}
if (!ncol(seqtab)) stop("Denoising produced no ASVs")
log_message("INFO", "Denoising/merging completed in ", elapsed(start))

if (cfg$chimera == "none") {
  seqtab_nochim <- seqtab
} else {
  log_message("INFO", "Removing chimeras with method=", cfg$chimera)
  seqtab_nochim <- dada2::removeBimeraDenovo(seqtab, method = cfg$chimera,
                                             multithread = threads, verbose = TRUE)
}
if (!ncol(seqtab_nochim)) stop("Chimera removal left no ASVs")

get_n <- function(x) sum(dada2::getUniques(x))
track <- data.frame(input = filter_out[, "reads.in"], filtered = filter_out[, "reads.out"],
                    denoisedF = vapply(dd_f, get_n, numeric(1)), row.names = samples,
                    check.names = FALSE)
if (mode == "pe") {
  track$denoisedR <- vapply(dd_r, get_n, numeric(1))
  track$merged <- vapply(mergers, get_n, numeric(1))
}
track$nonchim <- rowSums(seqtab_nochim[samples, , drop = FALSE])
track$retained_pct <- round(100 * track$nonchim / pmax(track$input, 1), 3)

saveRDS(seqtab_nochim, file.path(output_dir, "seqtab.nochim.rds"), compress = "xz")
write.table(track, file.path(output_dir, "track.summary.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
write.table(seqtab_nochim, file.path(output_dir, "seqtab.nochim.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
asv_seq <- colnames(seqtab_nochim)
writeLines(as.vector(rbind(paste0(">ASV", seq_along(asv_seq)), asv_seq)),
           file.path(output_dir, "ASVs.fasta"))

if (mode == "pe") {
  unlink(file.path(output_dir, c("suggestion.pe2se.note", "suggestion.is_pe.note")), force = TRUE)
  merge_ratio <- track$merged / pmax(track$input, 1)
  low_fraction <- mean(merge_ratio < 0.5)
  note <- c(sprintf("Samples with merged/input < 50%%: %d/%d (%.1f%%)",
                    sum(merge_ratio < 0.5), length(merge_ratio), 100 * low_fraction),
            sprintf("Median merged/input: %.1f%%", 100 * median(merge_ratio)))
  if (low_fraction >= 0.25) {
    note <- c(note, "Recommendation: inspect overlap and quality; consider forward-only SE analysis.")
    writeLines(note, file.path(output_dir, "suggestion.pe2se.note"))
  } else {
    note <- c(note, "Recommendation: PE retention is acceptable under the configured threshold.")
    writeLines(note, file.path(output_dir, "suggestion.is_pe.note"))
  }
} else {
  unlink(file.path(output_dir, c("suggestion.pe2se.note", "suggestion.is_pe.note")), force = TRUE)
}

if (!is.null(opt$classifier)) {
  classifier <- normalizePath(opt$classifier, mustWork = TRUE)
  log_message("INFO", "Assigning taxonomy with ", classifier)
  taxa <- dada2::assignTaxonomy(seqtab_nochim, classifier, minBoot = min_boot,
                                tryRC = !isTRUE(opt$no_try_rc),
                                multithread = threads, verbose = TRUE)
  write.table(taxa, file.path(output_dir, "taxonomy.tsv"),
              sep = "\t", quote = FALSE, col.names = NA)
  saveRDS(taxa, file.path(output_dir, "taxonomy.rds"), compress = "xz")
} else {
  unlink(file.path(output_dir, c("taxonomy.tsv", "taxonomy.rds")), force = TRUE)
}

capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
if (!isTRUE(opt$keep_filtered)) unlink(filtered_dir, recursive = TRUE, force = TRUE)
log_message("INFO", "Finished: ", nrow(seqtab_nochim), " samples, ",
            ncol(seqtab_nochim), " ASVs, elapsed ", elapsed(start_all))
