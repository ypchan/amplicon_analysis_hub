#!/usr/bin/env Rscript

# Infer approximate 16S V-region coverage by aligning sampled ASVs to a
# coordinate-compatible full-length 16S reference FASTA.

suppressPackageStartupMessages({
  library(getopt)
})

VERSION <- "2.0.0"
spec <- matrix(c(
  "ref",             "f", 1, "character",
  "rds_list",        "l", 1, "character",
  "rds",             "r", 1, "character",
  "out",             "o", 1, "character",
  "sample_size",     "n", 1, "integer",
  "id_threshold",    "i", 1, "double",
  "window_bp",       "w", 1, "integer",
  "overlap_min",      NA, 1, "integer",
  "vsearch_bin",     "v", 1, "character",
  "threads",         "t", 1, "integer",
  "vsearch_threads", "T", 1, "integer",
  "help",            "h", 0, "logical",
  "version",         "V", 0, "logical"
), byrow = TRUE, ncol = 4)

usage <- function(status = 0L) {
  cat("Infer approximate 16S variable-region coverage from DADA2 sequence tables.\n\n")
  cat("Usage:\n  infer_16s_regions.R --ref full_length_16S.fasta (--rds_list FILE | --rds FILES) [options]\n\n")
  cat("Options:\n")
  cat("  -f, --ref FASTA             Full-length, coordinate-compatible 16S references (required)\n")
  cat("  -l, --rds_list FILE         One seqtab.nochim.rds path per line\n")
  cat("  -r, --rds FILES             Comma-separated RDS paths\n")
  cat("  -o, --out TSV               Output (default: 16s_regions.tsv)\n")
  cat("  -n, --sample_size INT       Abundance-weighted ASVs per table (default: 1000)\n")
  cat("  -i, --id_threshold NUM      vsearch identity fraction (default: 0.80)\n")
  cat("  -w, --window_bp INT         Coordinate cluster half-window (default: 30)\n")
  cat("      --overlap_min INT       V-region overlap needed in bp (default: 10)\n")
  cat("  -v, --vsearch_bin FILE      vsearch executable (default: vsearch)\n")
  cat("  -t, --threads INT           Concurrent R workers (default: 1)\n")
  cat("  -T, --vsearch_threads INT   Threads per alignment (default: 1)\n")
  cat("  -h, --help                  Show help\n")
  cat("  -V, --version               Show version\n\n")
  cat("CPU use is approximately --threads × --vsearch_threads. Coordinates follow\n")
  cat("the E. coli 16S convention and are approximate for divergent taxa. This tool\n")
  cat("does not classify ITS1/ITS2, whose boundaries require different references.\n")
  quit(status = status)
}

opt <- getopt(spec)
if (isTRUE(opt$help)) usage(0L)
if (isTRUE(opt$version)) { cat("infer_16s_regions.R ", VERSION, "\n", sep = ""); quit(status = 0L) }
if (is.null(opt$ref) || (is.null(opt$rds_list) && is.null(opt$rds))) usage(2L)

value_or <- function(value, default) if (is.null(value)) default else value
output <- value_or(opt$out, "16s_regions.tsv")
sample_size <- as.integer(value_or(opt$sample_size, 1000L))
identity <- as.numeric(value_or(opt$id_threshold, 0.80))
window_bp <- as.integer(value_or(opt$window_bp, 30L))
overlap_min <- as.integer(value_or(opt$overlap_min, 10L))
vsearch <- value_or(opt$vsearch_bin, "vsearch")
workers <- as.integer(value_or(opt$threads, 1L))
vsearch_threads <- as.integer(value_or(opt$vsearch_threads, 1L))
integer_values <- c(sample_size, window_bp, overlap_min, workers, vsearch_threads)
if (any(!is.finite(integer_values)) || sample_size < 1L || window_bp < 0L ||
    overlap_min < 1L || workers < 1L || vsearch_threads < 1L) {
  stop("sample/thread/overlap values are outside their valid range")
}
if (!is.finite(identity) || identity <= 0 || identity > 1) stop("--id_threshold must be in (0,1]")
reference <- normalizePath(opt$ref, mustWork = TRUE)
if (!nzchar(Sys.which(vsearch))) stop("vsearch executable not found: ", vsearch)

regions <- data.frame(
  name = paste0("V", 1:9),
  start = c(69L, 137L, 433L, 576L, 822L, 986L, 1117L, 1243L, 1435L),
  end = c(99L, 242L, 497L, 682L, 879L, 1043L, 1173L, 1294L, 1465L)
)

parse_paths <- function() {
  paths <- character()
  if (!is.null(opt$rds_list)) {
    lines <- trimws(readLines(opt$rds_list, warn = FALSE))
    paths <- c(paths, lines[nzchar(lines) & !startsWith(lines, "#")])
  }
  if (!is.null(opt$rds)) paths <- c(paths, trimws(strsplit(opt$rds, ",", fixed = TRUE)[[1]]))
  unique(paths[nzchar(paths)])
}

write_asv_sample <- function(rds_path, fasta_path) {
  table <- readRDS(rds_path)
  if (!is.matrix(table) || is.null(colnames(table))) stop("Not a sequence-table matrix: ", rds_path)
  abundance <- colSums(table)
  valid <- abundance > 0 & nzchar(colnames(table))
  sequences <- colnames(table)[valid]
  abundance <- abundance[valid]
  if (!length(sequences)) stop("No positive-abundance ASVs: ", rds_path)
  count <- min(sample_size, length(sequences))
  seed <- sum(utf8ToInt(normalizePath(rds_path))) %% .Machine$integer.max
  set.seed(seed)
  chosen <- sample(seq_along(sequences), count, replace = FALSE, prob = abundance)
  lines <- as.vector(rbind(paste0(">ASV", seq_along(chosen)), sequences[chosen]))
  writeLines(lines, fasta_path)
  count
}

region_label <- function(start, end) {
  if (!is.finite(start) || !is.finite(end)) return("coordinates_unknown")
  overlap <- pmax(0L, pmin(end, regions$end) - pmax(start, regions$start) + 1L)
  covered <- which(overlap >= overlap_min)
  if (!length(covered)) return("no_V_region_match")
  if (length(covered) > 1L && any(diff(covered) != 1L)) {
    return(paste0("mixed_", paste(regions$name[covered], collapse = "_")))
  }
  if (length(covered) == 1L) regions$name[covered] else paste0(regions$name[min(covered)], "-", sub("^V", "", regions$name[max(covered)]))
}

process_one <- function(path) {
  full_path <- normalizePath(path, mustWork = TRUE)
  temp_dir <- tempfile("infer-16s-")
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE, force = TRUE), add = TRUE)
  query <- file.path(temp_dir, "asv.fasta")
  hits_path <- file.path(temp_dir, "hits.tsv")
  sampled <- write_asv_sample(full_path, query)
  args <- c("--usearch_global", query, "--db", reference, "--id", identity,
            "--strand", "both", "--maxaccepts", "1", "--maxrejects", "32",
            "--top_hits_only", "--threads", vsearch_threads,
            "--userout", hits_path,
            "--userfields", "query+target+id+alnlen+qlo+qhi+tlo+thi")
  status <- suppressWarnings(system2(vsearch, args = as.character(args), stdout = FALSE, stderr = FALSE))
  if (!identical(status, 0L)) stop("vsearch failed for ", full_path)
  identifier <- paste(basename(dirname(full_path)), sub("\\.rds$", "", basename(full_path)), sep = "/")
  if (!file.exists(hits_path) || file.size(hits_path) == 0L) {
    return(data.frame(identifier, sampled_asvs = sampled, aligned_asvs = 0L,
                      aligned_percent = 0, majority_cluster_percent = 0,
                      median_start = NA, median_end = NA, region = "no_hits"))
  }
  hits <- read.delim(hits_path, header = FALSE,
                     col.names = c("query", "target", "identity", "aln_len", "qlo", "qhi", "tlo", "thi"))
  starts <- pmin(hits$tlo, hits$thi)
  ends <- pmax(hits$tlo, hits$thi)
  med_start <- median(starts)
  med_end <- median(ends)
  majority <- mean(abs(starts - med_start) <= window_bp & abs(ends - med_end) <= window_bp)
  data.frame(identifier, sampled_asvs = sampled, aligned_asvs = nrow(hits),
             aligned_percent = round(100 * nrow(hits) / sampled, 3),
             majority_cluster_percent = round(100 * majority, 3),
             median_start = med_start, median_end = med_end,
             region = region_label(med_start, med_end))
}

paths <- parse_paths()
if (!length(paths)) stop("No RDS paths were provided")
missing <- paths[!file.exists(paths)]
if (length(missing)) stop("RDS path(s) not found: ", paste(missing, collapse = ", "))
message("Files: ", length(paths), "; approximate CPUs: ", workers * vsearch_threads)
safe_process <- function(path) tryCatch(process_one(path), error = function(error) {
  data.frame(identifier = path, sampled_asvs = NA, aligned_asvs = NA,
             aligned_percent = NA, majority_cluster_percent = NA,
             median_start = NA, median_end = NA,
             region = paste0("ERROR: ", conditionMessage(error)))
})
if (.Platform$OS.type == "windows" || workers == 1L) {
  result <- lapply(paths, safe_process)
} else {
  result <- parallel::mclapply(paths, safe_process, mc.cores = workers)
}
output_dir <- dirname(output)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
write.table(do.call(rbind, result), output, sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote: ", output)
