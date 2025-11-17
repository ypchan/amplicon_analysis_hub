#!/usr/bin/env Rscript

## ------------------------------------------------------------
## Infer amplicon region (V1..V9, V1-2, V3-4, V1-8, etc.)
## for many seqtab.nochim.rds files by:
##   - sampling ASVs
##   - aligning to a full-length 16S reference (E. coli numbering)
##   - using median start/end coordinates and V1–V9 ranges
##
## Output (tab-delimited):
##   identifier   majority_hit_percent   region
##
## Parallelism:
##   --threads          : how many R workers (seqtab-level) in parallel
##   --vsearch_threads  : threads per vsearch call
##
## Coordinate model:
##   E. coli 16S variable regions:
##     V1:  69–99
##     V2:  137–242
##     V3:  433–497
##     V4:  576–682
##     V5:  822–879
##     V6:  986–1043
##     V7:  1117–1173
##     V8:  1243–1294
##     V9:  1435–1465
## ------------------------------------------------------------

suppressPackageStartupMessages({
  library(getopt)
  library(dada2)
  library(readr)
  library(dplyr)
  library(purrr)
  library(stringr)
  library(parallel)
})

## =========================
## 0. Command line options
## =========================

spec <- matrix(c(
  "help",            "h", 0, "logical",   "Show this help message",
  "ref",             "f", 1, "character", "Full-length 16S reference FASTA (bacteria + archaea, >1200bp, E. coli coordinates)",
  "rds_list",        "l", 1, "character", "Text file with one seqtab.nochim.rds path per line",
  "rds",             "r", 1, "character", "One or more RDS files (comma-separated)",
  "out",             "o", 1, "character", "Output TSV file (default: amplicon_regions.tsv)",
  "sample_size",     "n", 1, "integer",   "Number of ASVs to sample per seqtab (default: 1000)",
  "id_threshold",    "i", 1, "double",    "Minimum identity threshold for vsearch (default: 0.8)",
  "window_bp",       "w", 1, "integer",   "Window (bp) around median start/end for majority region (default: 30)",
  "vsearch_bin",     "v", 1, "character", "Path to vsearch executable (default: vsearch)",
  "threads",         "t", 1, "integer",   "Number of R worker processes (seqtab-level) (default: 1)",
  "vsearch_threads", "T", 1, "integer",   "Number of threads per vsearch call (default: 1)"
), byrow = TRUE, ncol = 5)

opt <- getopt(spec)

print_help_and_quit <- function() {
  cat("\nUsage:\n")
  cat("  Rscript infer_amplicon_region.R --ref full_length_16S.fasta \\\n")
  cat("      --rds_list seqtab_list.txt --out amplicon_regions.tsv [options]\n\n")
  cat("Options:\n")
  print(getopt(spec, usage = TRUE))
  cat("\nNotes:\n")
  cat("  --rds can be a comma-separated list of extra RDS files.\n")
  cat("  Either --rds_list, --rds, or both must be provided.\n")
  cat("  --threads controls how many R workers run in parallel (seqtab-level).\n")
  cat("  --vsearch_threads controls how many threads each vsearch call uses.\n\n")
  quit(status = 1)
}

if (!is.null(opt$help)) {
  print_help_and_quit()
}

## Required: reference FASTA
if (is.null(opt$ref)) {
  cat("ERROR: --ref is required.\n")
  print_help_and_quit()
}

## Optional options with defaults
ref_fasta       <- opt$ref
rds_list_file   <- if (!is.null(opt$rds_list)) opt$rds_list else NA_character_
output_tsv      <- if (!is.null(opt$out)) opt$out else "amplicon_regions.tsv"
sample_size     <- if (!is.null(opt$sample_size)) opt$sample_size else 1000L
id_threshold    <- if (!is.null(opt$id_threshold)) opt$id_threshold else 0.8
window_bp       <- if (!is.null(opt$window_bp)) opt$window_bp else 30L
vsearch_bin     <- if (!is.null(opt$vsearch_bin)) opt$vsearch_bin else "vsearch"
threads         <- if (!is.null(opt$threads)) opt$threads else 1L
vsearch_threads <- if (!is.null(opt$vsearch_threads)) opt$vsearch_threads else 1L

if (threads < 1) threads <- 1L
if (vsearch_threads < 1) vsearch_threads <- 1L

## Parse extra RDS files (comma-separated)
rds_files_extra <- character(0)
if (!is.null(opt$rds)) {
  rds_files_extra <- unlist(strsplit(opt$rds, ",", fixed = TRUE))
  rds_files_extra <- trimws(rds_files_extra)
  rds_files_extra <- rds_files_extra[nzchar(rds_files_extra)]
}

## =========================
## 1. Variable region model
## =========================

## Variable region coordinates (E. coli 16S numbering)
## V1–V9 ranges:
## V1:  69–99
## V2:  137–242
## V3:  433–497
## V4:  576–682
## V5:  822–879
## V6:  986–1043
## V7:  1117–1173
## V8:  1243–1294
## V9:  1435–1465

var_regions <- data.frame(
  name  = paste0("V", 1:9),
  start = c(  69,  137,  433,  576,  822,   986,  1117,  1243,  1435),
  end   = c(  99,  242,  497,  682,  879,  1043,  1173,  1294,  1465),
  stringsAsFactors = FALSE
)

## =========================
## 2. Helper functions
## =========================

run_system <- function(cmd, args = character(), ...) {
  status <- system2(cmd, args = args, ...)
  if (!is.null(status) && status != 0) {
    stop("Command failed: ", cmd, " ", paste(args, collapse = " "))
  }
}

read_rds_list_file <- function(path) {
  if (is.na(path) || !file.exists(path)) {
    return(character(0))
  }
  lines <- readr::read_lines(path)
  lines <- lines[nzchar(lines)]
  return(lines)
}

sample_asvs_from_seqtab <- function(rds_path, out_fasta, sample_size = 1000L) {
  message("  Reading seqtab: ", rds_path)
  seqtab <- readRDS(rds_path)

  if (is.null(colnames(seqtab))) {
    stop("No ASV sequences (colnames) found in seqtab: ", rds_path)
  }

  asv_seqs <- colnames(seqtab)
  abund    <- colSums(seqtab)

  probs <- abund / sum(abund)

  n_asv    <- length(asv_seqs)
  n_sample <- min(sample_size, n_asv)

  set.seed(123)  # for reproducibility
  idx <- sample(seq_len(n_asv), size = n_sample, replace = FALSE, prob = probs)

  asv_seqs_sample <- asv_seqs[idx]

  con <- file(out_fasta, "w")
  on.exit(close(con), add = TRUE)

  for (i in seq_along(asv_seqs_sample)) {
    cat(">", "ASV_", i, "\n", asv_seqs_sample[i], "\n", file = con, sep = "")
  }
}

run_vsearch_alignment <- function(query_fasta, ref_fasta, out_tsv,
                                  vsearch = "vsearch",
                                  id_threshold = 0.8,
                                  vsearch_threads = 1L) {
  args <- c(
    "--usearch_global", query_fasta,
    "--db", ref_fasta,
    "--id", as.character(id_threshold),
    "--strand", "both",
    "--maxaccepts", "1",
    "--maxhits", "1",
    "--top_hits_only",
    "--threads", as.character(vsearch_threads),
    "--blast6out", out_tsv
  )
  message("  Running vsearch: ", vsearch, " ", paste(args, collapse = " "))
  run_system(vsearch, args = args)
}

parse_blast6_regions <- function(blast6_path) {
  if (!file.exists(blast6_path)) {
    stop("BLAST6 file not found: ", blast6_path)
  }

  if (file.size(blast6_path) == 0) {
    return(tibble(
      qseqid = character(0),
      sseqid = character(0),
      sstart = integer(0),
      send = integer(0),
      region_start = integer(0),
      region_end = integer(0)
    ))
  }

  tbl <- readr::read_tsv(
    blast6_path,
    col_names = FALSE,
    show_col_types = FALSE,
    progress = FALSE
  )
  if (nrow(tbl) == 0) {
    return(tibble(
      qseqid = character(0),
      sseqid = character(0),
      sstart = integer(0),
      send = integer(0),
      region_start = integer(0),
      region_end = integer(0)
    ))
  }

  colnames(tbl)[1:10] <- c(
    "qseqid", "sseqid", "pident", "length", "mismatch",
    "gapopen", "qstart", "qend", "sstart", "send"
  )

  tbl <- tbl %>%
    mutate(
      region_start = pmin(sstart, send),
      region_end   = pmax(sstart, send)
    ) %>%
    select(qseqid, sseqid, sstart, send, region_start, region_end)

  return(tbl)
}

infer_majority_region <- function(region_df, window_bp = 30L) {
  if (nrow(region_df) == 0) {
    return(list(
      med_start = NA_real_,
      med_end = NA_real_,
      majority_fraction = 0.0,
      amplicon_len = NA_real_
    ))
  }

  starts <- region_df$region_start
  ends   <- region_df$region_end

  med_start    <- stats::median(starts)
  med_end      <- stats::median(ends)
  amplicon_len <- med_end - med_start + 1

  in_majority <- (abs(starts - med_start) <= window_bp) &
                 (abs(ends   - med_end)   <= window_bp)

  majority_fraction <- sum(in_majority) / length(in_majority)

  return(list(
    med_start = med_start,
    med_end = med_end,
    majority_fraction = majority_fraction,
    amplicon_len = amplicon_len
  ))
}

classify_region_by_coords <- function(med_start, med_end, overlap_min_bp = 10L) {
  # Classify amplicon region based on median start/end coordinates
  # relative to E. coli 16S variable regions (V1–V9).
  #
  # Returns examples:
  #   "V4"
  #   "V3-4"
  #   "V1-3"
  #   "V1-8"
  #   "mixed_V1_V3"          (non-contiguous coverage)
  #   "no_V_region_match"
  #   "coords_unknown"

  if (is.na(med_start) || is.na(med_end)) {
    return("coords_unknown")
  }

  s <- min(med_start, med_end)
  e <- max(med_start, med_end)

  covered <- logical(nrow(var_regions))

  for (i in seq_len(nrow(var_regions))) {
    vs <- var_regions$start[i]
    ve <- var_regions$end[i]
    ov_start <- max(s, vs)
    ov_end   <- min(e, ve)
    ov_len   <- ov_end - ov_start + 1
    covered[i] <- ov_len >= overlap_min_bp
  }

  if (!any(covered)) {
    return("no_V_region_match")
  }

  covered_idx   <- which(covered)
  covered_names <- var_regions$name[covered_idx]
  covered_nums  <- as.integer(sub("^V", "", covered_names))
  covered_nums  <- sort(covered_nums)

  if (length(covered_nums) > 1 &&
      any(diff(covered_nums) != 1)) {
    return(paste0("mixed_", paste0("V", covered_nums, collapse = "_")))
  }

  if (length(covered_nums) == 1) {
    return(paste0("V", covered_nums))
  } else {
    return(paste0("V", min(covered_nums), "-", max(covered_nums)))
  }
}

process_one_seqtab <- function(rds_path,
                               ref_fasta,
                               sample_size = 1000L,
                               vsearch = "vsearch",
                               id_threshold = 0.8,
                               window_bp = 30L,
                               vsearch_threads = 1L) {
  rds_path <- normalizePath(rds_path)
  message("Processing: ", rds_path)

  identifier <- basename(rds_path)
  identifier <- sub("\\.rds$", "", identifier)

  tmpdir <- tempfile(pattern = "amplicon_tmp_")
  dir.create(tmpdir, showWarnings = FALSE)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  asv_fasta  <- file.path(tmpdir, "asv_sample.fasta")
  blast6_out <- file.path(tmpdir, "hits.tsv")

  sample_asvs_from_seqtab(rds_path, asv_fasta, sample_size = sample_size)

  run_vsearch_alignment(
    query_fasta      = asv_fasta,
    ref_fasta        = ref_fasta,
    out_tsv          = blast6_out,
    vsearch          = vsearch,
    id_threshold     = id_threshold,
    vsearch_threads  = vsearch_threads
  )

  region_df <- parse_blast6_regions(blast6_out)

  if (nrow(region_df) == 0) {
    return(tibble(
      identifier = identifier,
      majority_hit_percent = 0,
      region = "no_hits"
    ))
  }

  region_info      <- infer_majority_region(region_df, window_bp = window_bp)
  majority_percent <- region_info$majority_fraction * 100

  region_label <- classify_region_by_coords(
    med_start = region_info$med_start,
    med_end   = region_info$med_end,
    overlap_min_bp = 10L
  )

  tibble(
    identifier = identifier,
    majority_hit_percent = majority_percent,
    region = region_label
  )
}

## =========================
## 3. Collect input RDS list
## =========================

rds_from_file <- read_rds_list_file(rds_list_file)
rds_all <- unique(c(rds_from_file, rds_files_extra))

if (length(rds_all) == 0) {
  cat("ERROR: No RDS files provided. Use --rds_list or --rds.\n")
  print_help_and_quit()
}

message("Total RDS files to process: ", length(rds_all))
message("R workers (seqtab-level): ", threads)
message("vsearch threads per worker: ", vsearch_threads)
message("Total logical CPUs needed (approx): ", threads * vsearch_threads)

## =========================
## 4. Sanity checks
## =========================

if (!file.exists(ref_fasta)) {
  stop("Reference FASTA not found: ", ref_fasta)
}

vsearch_check <- suppressWarnings(
  system2(vsearch_bin, args = "--version", stdout = TRUE, stderr = TRUE)
)
if (!is.null(attr(vsearch_check, "status")) &&
    attr(vsearch_check, "status") != 0) {
  stop("Cannot run vsearch. Please check --vsearch_bin: ", vsearch_bin)
}

## =========================
## 5. Main loop (parallel)
## =========================

worker_fun <- function(rds) {
  tryCatch(
    process_one_seqtab(
      rds_path         = rds,
      ref_fasta        = ref_fasta,
      sample_size      = sample_size,
      vsearch          = vsearch_bin,
      id_threshold     = id_threshold,
      window_bp        = window_bp,
      vsearch_threads  = vsearch_threads
    ),
    error = function(e) {
      message("Error processing ", rds, ": ", conditionMessage(e))
      tibble(
        identifier = basename(rds),
        majority_hit_percent = 0,
        region = "ERROR"
      )
    }
  )
}

if (.Platform$OS.type == "windows" || threads == 1L) {
  results_list <- lapply(rds_all, worker_fun)
} else {
  results_list <- mclapply(rds_all, worker_fun, mc.cores = threads)
}

results <- bind_rows(results_list)

## =========================
## 6. Write output
## =========================

readr::write_tsv(results, output_tsv)
message("Done. Result written to: ", output_tsv)
