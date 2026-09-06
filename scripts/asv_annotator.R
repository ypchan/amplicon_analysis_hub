#!/usr/bin/env Rscript

# Assign taxonomy to DADA2 ASVs and write an ASV-by-sample count table.

suppressPackageStartupMessages({
  library(getopt)
})

VERSION <- "2.0.0"
spec <- matrix(c(
  "seqtab_rds",    "s", 1, "character",
  "output",        "o", 1, "character",
  "train_fasta",   "t", 1, "character",
  "species_fasta", "p", 1, "character",
  "marker",        "M", 1, "character",
  "min_boot",      "b", 1, "integer",
  "threads",       "n", 1, "integer",
  "tax_delim",     "d", 1, "character",
  "no_try_rc",      NA, 0, "logical",
  "help",          "h", 0, "logical",
  "version",       "V", 0, "logical"
), byrow = TRUE, ncol = 4)

usage <- function(status = 0L) {
  cat("Assign taxonomy and combine it with a DADA2 ASV table.\n\n")
  cat("Usage:\n  asv_annotator.R -s seqtab.nochim.rds -t TRAINING.fa.gz -o ASVs.tsv [options]\n\n")
  cat("Required:\n")
  cat("  -s, --seqtab_rds RDS      DADA2 sample-by-sequence matrix\n")
  cat("  -t, --train_fasta FASTA   assignTaxonomy training FASTA (GTDB/SILVA/UNITE)\n")
  cat("  -o, --output FILE         Output .tsv or .csv\n\n")
  cat("Options:\n")
  cat("  -p, --species_fasta FILE  Exact species reference for addSpecies (usually 16S)\n")
  cat("  -M, --marker NAME         16s|its|other; documentation metadata (default: 16s)\n")
  cat("  -b, --min_boot INT        Minimum taxonomy bootstrap, 0..100 (default: 50)\n")
  cat("  -n, --threads INT         Taxonomy threads (default: min(8, detected cores))\n")
  cat("  -d, --tax_delim STR       Collapsed-taxonomy delimiter (default: ;)\n")
  cat("      --no_try_rc           Disable reverse-complement classification\n")
  cat("  -h, --help                Show help\n")
  cat("  -V, --version             Show version\n\n")
  cat("Use a marker-compatible training set: GTDB/SILVA/RDP for 16S and UNITE for\n")
  cat("fungal ITS. Database release and sequence orientation are recorded in your\n")
  cat("workflow metadata, not inferred by this command.\n")
  quit(status = status)
}

opt <- getopt(spec)
if (isTRUE(opt$help)) usage(0L)
if (isTRUE(opt$version)) { cat("asv_annotator.R ", VERSION, "\n", sep = ""); quit(status = 0L) }
if (is.null(opt$seqtab_rds) || is.null(opt$train_fasta) || is.null(opt$output)) usage(2L)
if (!requireNamespace("dada2", quietly = TRUE)) stop("R package 'dada2' is required")

marker <- tolower(if (is.null(opt$marker)) "16s" else opt$marker)
if (!marker %in% c("16s", "its", "other")) stop("--marker must be 16s, its, or other")
min_boot <- as.integer(if (is.null(opt$min_boot)) 50L else opt$min_boot)
detected <- suppressWarnings(parallel::detectCores(logical = TRUE))
if (is.na(detected)) detected <- 1L
threads <- as.integer(if (is.null(opt$threads)) min(8L, detected) else opt$threads)
delimiter <- if (is.null(opt$tax_delim)) ";" else opt$tax_delim
if (!is.finite(min_boot) || min_boot < 0L || min_boot > 100L) stop("--min_boot must be in [0,100]")
if (!is.finite(threads) || threads < 1L) stop("--threads must be >= 1")
if (!nzchar(delimiter)) stop("--tax_delim cannot be empty")

seqtab_path <- normalizePath(opt$seqtab_rds, mustWork = TRUE)
training_path <- normalizePath(opt$train_fasta, mustWork = TRUE)
seqtab <- readRDS(seqtab_path)
if (!is.matrix(seqtab) || is.null(rownames(seqtab)) || is.null(colnames(seqtab))) {
  stop("seqtab RDS must be a sample-by-ASV matrix with row and column names")
}
if (anyDuplicated(rownames(seqtab)) || anyDuplicated(colnames(seqtab))) stop("seqtab names must be unique")
if (any(!is.finite(seqtab)) || any(seqtab < 0)) stop("seqtab contains invalid counts")

message("Assigning taxonomy to ", ncol(seqtab), " ASVs with ", threads, " thread(s)")
taxa <- dada2::assignTaxonomy(colnames(seqtab), training_path, minBoot = min_boot,
                              tryRC = !isTRUE(opt$no_try_rc), multithread = threads,
                              outputBootstraps = FALSE, verbose = TRUE)
if (!is.null(opt$species_fasta)) {
  if (marker == "its") warning("addSpecies exact matching is usually designed for 16S species references; verify ITS suitability")
  species_path <- normalizePath(opt$species_fasta, mustWork = TRUE)
  taxa <- dada2::addSpecies(taxa, species_path, tryRC = !isTRUE(opt$no_try_rc))
}

tax_character <- as.matrix(taxa)
tax_character[is.na(tax_character)] <- ""
collapsed <- apply(tax_character, 1L, function(row) paste(row[nzchar(row)], collapse = delimiter))
counts <- as.data.frame(t(seqtab), check.names = FALSE)
output <- data.frame(ASV = paste0("ASV", seq_len(nrow(counts))),
                     Sequence = rownames(counts), Taxonomy = unname(collapsed),
                     taxa, counts, check.names = FALSE, stringsAsFactors = FALSE)

output_path <- opt$output
output_dir <- dirname(output_path)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
separator <- if (grepl("\\.csv$", output_path, ignore.case = TRUE)) "," else "\t"
write.table(output, output_path, sep = separator, quote = separator == ",",
            row.names = FALSE, na = "")
message("Wrote ", nrow(output), " ASVs × ", ncol(seqtab), " samples: ", output_path)
