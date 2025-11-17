#!/usr/bin/env Rscript

# asv_annotator.R
# Build an ASV count matrix with taxonomy for DADA2 results.
#
# Output columns:
#   1) ASV           (ASV sequence; one unique sequence per row)
#   2) Taxonomy      (joined ranks: Kingdom;Phylum;Class;Order;Family;Genus;Species)
#   3..N) sample counts (one column per sample; preserved in the same order as input seqtab)
#
# Notes:
#   - This script ONLY annotates using reference FASTA via assignTaxonomy()
#     (and optional addSpecies()). It does NOT merge an existing taxa.rds.
#   - Compatible with dada2 variants where assignTaxonomy takes either
#     'refFasta' or 'trainingSet'; same for addSpecies().
#   - All comments are in English by request.
#
# Example:
#   Rscript asv_annotator.R \
#     --seqtab_rds seqtab.nochim.rds \
#     --output_csv asv_counts.csv \
#     --train_fasta silva_nr99_v138.2_toSpecies_trainset.fa.gz \
#     --species_fasta silva_v138.2_assignSpecies.fa.gz \
#     --minBoot 50 --threads 8 --tax_delim ";"
#
# Required flags:
#   --seqtab_rds    Path to seqtab.nochim RDS (rows=samples, cols=ASVs)
#   --output_csv    Path to write the output CSV
#   --train_fasta   Reference training FASTA for assignTaxonomy()
#
# Optional flags:
#   --species_fasta Path to species FASTA for addSpecies()
#   --minBoot       Integer; minimum bootstrap for assignTaxonomy() [default: 50]
#   --threads       Integer; threads for assignTaxonomy() [default: detectCores-1]
#   --tax_delim     Character; taxonomy rank delimiter [default: ";"]
#   --help          Show help and exit
#
# Exit codes:
#   0 success; 1 usage error; >1 runtime errors.

suppressPackageStartupMessages({
    library(getopt)
    library(dada2)
    library(dplyr)
    library(readr)
    library(tibble)
})

# --------------------------- CLI (getopt) --------------------------------------
spec <- matrix(c(
    "seqtab_rds",    "s", 1, "character",
    "output_csv",    "o", 1, "character",
    "train_fasta",   "t", 1, "character",
    "species_fasta", "p", 2, "character",
    "minBoot",       "b", 2, "integer",
    "threads",       "n", 2, "integer",
    "tax_delim",     "d", 2, "character",
    "help",          "h", 0, "logical"
), byrow = TRUE, ncol = 4)

opt <- getopt(spec)

print_usage <- function() {
    cat("
asv_annotator.R

Required:
  --seqtab_rds    <path>  DADA2 seqtab.nochim RDS (rows=samples, cols=ASVs)
  --output_csv    <path>  Output CSV path
  --train_fasta   <path>  Reference training FASTA (e.g., SILVA) for assignTaxonomy()

Optional:
  --species_fasta <path>  Species FASTA for addSpecies()
  --minBoot       <int>   Minimum bootstrap for assignTaxonomy() [default: 50]
  --threads       <int>   Threads for assignTaxonomy() [default: detectCores-1]
  --tax_delim     <char>  Delimiter used to join taxonomy ranks [default: ';']
  --help                  Show this help

Example:
  Rscript asv_annotator.R \\
    --seqtab_rds seqtab.nochim.rds \\
    --output_csv asv_counts.csv \\
    --train_fasta silva_nr99_v138.2_toSpecies_trainset.fa.gz \\
    --species_fasta silva_v138.2_assignSpecies.fa.gz \\
    --minBoot 50 --threads 8 --tax_delim ';'
\n")
}

if (!is.null(opt$help)) { print_usage(); quit(status = 0) }

req_missing <- c(
    if (is.null(opt$seqtab_rds)) "--seqtab_rds" else NULL,
    if (is.null(opt$output_csv)) "--output_csv" else NULL,
    if (is.null(opt$train_fasta)) "--train_fasta" else NULL
)
if (length(req_missing) > 0) {
    cat("Missing required options: ", paste(req_missing, collapse = ", "), "\n\n", sep = "")
    print_usage(); quit(status = 1)
}

minBoot  <- if (!is.null(opt$minBoot)) as.integer(opt$minBoot) else 50L
threads  <- if (!is.null(opt$threads)) as.integer(opt$threads) else max(1L, parallel::detectCores() - 1L)
tax_delim <- if (!is.null(opt$tax_delim)) as.character(opt$tax_delim) else ";"

# --------------------------- Load seqtab ---------------------------------------
if (!file.exists(opt$seqtab_rds)) stop("seqtab RDS not found: ", opt$seqtab_rds)
seqtab <- readRDS(opt$seqtab_rds)
if (!is.matrix(seqtab)) seqtab <- as.matrix(seqtab)

# Validate structure: rows = samples, cols = ASVs (sequences as colnames)
if (is.null(colnames(seqtab)) || any(!nzchar(colnames(seqtab)))) {
    stop("seqtab has no valid column names; columns must be ASV sequences.")
}

# Ensure integer-like counts (avoid scientific notation in output)
storage.mode(seqtab) <- "integer"

# Transpose to ASV rows × sample columns (preferred output orientation)
count_mat <- t(seqtab)  # rows = ASV, cols = samples

# Remember sample column order (exactly as in input)
sample_cols <- colnames(seqtab)

# Defensive checks on ASV keys
asv_keys <- rownames(count_mat)
if (anyDuplicated(asv_keys) > 0) {
    dupn <- sum(duplicated(asv_keys))
    stop("Detected duplicated ASV sequences (rownames after transpose): ", dupn,
         ". Each ASV sequence must be unique.")
}

# --------------------------- Taxonomy ------------------------------------------
expected_ranks <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")

if (!file.exists(opt$train_fasta)) stop("train_fasta not found: ", opt$train_fasta)
asv_vec <- asv_keys

message("Running assignTaxonomy() [threads=", threads, ", minBoot=", minBoot, "] ...")

# Version-robust call to assignTaxonomy: accept 'refFasta' or 'trainingSet'
at_formals <- names(formals(dada2::assignTaxonomy))
use_refFasta    <- "refFasta"    %in% at_formals
use_trainingSet <- "trainingSet" %in% at_formals

if (use_refFasta) {
    taxa <- assignTaxonomy(asv_vec,
                           refFasta = opt$train_fasta,
                           tryRC = TRUE,
                           multithread = threads,
                           minBoot = minBoot)
} else if (use_trainingSet) {
    taxa <- assignTaxonomy(asv_vec,
                           trainingSet = opt$train_fasta,
                           tryRC = TRUE,
                           multithread = threads,
                           minBoot = minBoot)
} else {
    # Fallback: positional second argument
    taxa <- assignTaxonomy(asv_vec,
                           opt$train_fasta,
                           tryRC = TRUE,
                           multithread = threads,
                           minBoot = minBoot)
}

# Optional species refinement
if (!is.null(opt$species_fasta)) {
    if (!file.exists(opt$species_fasta)) stop("species_fasta not found: ", opt$species_fasta)
    message("Running addSpecies() ...")
    
    as_formals <- names(formals(dada2::addSpecies))
    if ("refFasta" %in% as_formals) {
        taxa <- addSpecies(taxa, refFasta = opt$species_fasta, tryRC = TRUE)
    } else if ("trainingSet" %in% as_formals) {
        taxa <- addSpecies(taxa, trainingSet = opt$species_fasta, tryRC = TRUE)
    } else {
        taxa <- addSpecies(taxa, opt$species_fasta, tryRC = TRUE)
    }
}

# Keep expected ranks only (if present)
keep <- intersect(expected_ranks, colnames(taxa))
if (length(keep) == 0) stop("No expected taxonomy ranks present in taxonomy result.")

# Build a data.frame: ASV + collapsed Taxonomy
tax_df <- as.data.frame(taxa[, keep, drop = FALSE], stringsAsFactors = FALSE)
tax_df$ASV <- rownames(taxa)

# Fast vectorized collapse of ranks per ASV; remove empty/NA levels before joining
collapse_tax <- function(m, delim = ";") {
    # m is a character matrix of ranks, rows = ASV, cols = ranks
    m[is.na(m) | m == ""] <- NA
    # Apply over rows without creating row-wise groups (faster than rowwise())
    vapply(seq_len(nrow(m)), function(i) {
        x <- m[i, ]
        x <- x[!is.na(x)]
        if (length(x) == 0) "" else paste(x, collapse = delim)
    }, FUN.VALUE = character(1))
}
tax_df$Taxonomy <- collapse_tax(as.matrix(tax_df[, keep, drop = FALSE]), delim = tax_delim)

tax_df <- tax_df[, c("ASV", "Taxonomy"), drop = FALSE]

# --------------------------- Assemble output -----------------------------------
# Convert count matrix (ASV rows × sample cols) to data.frame
asv_df <- as.data.frame(count_mat, check.names = FALSE) |>
    rownames_to_column(var = "ASV")

# Right-join to keep all ASVs even if taxonomy is empty
out_df <- tax_df |>
    right_join(asv_df, by = "ASV")

# Ensure 'Taxonomy' exists and is character; replace NA with empty string
if (!"Taxonomy" %in% names(out_df)) out_df$Taxonomy <- ""
out_df$Taxonomy[is.na(out_df$Taxonomy)] <- ""

# Final column order: ASV, Taxonomy, <samples in original order>
out_df <- out_df |>
    select(c("ASV", "Taxonomy", sample_cols))

# --------------------------- Write ---------------------------------------------
# Avoid scientific notation; keep large strings intact
options(scipen = 999)

# readr::write_csv preserves column order and UTF-8; counts will be written as integers
readr::write_csv(out_df, opt$output_csv)

message("Done. Wrote: ", opt$output_csv)
