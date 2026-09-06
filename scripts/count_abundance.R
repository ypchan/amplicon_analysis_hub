#!/usr/bin/env Rscript

# Build ASV and taxon abundance tables without expanding a dense sequence table
# into a potentially enormous sample-ASV long table.

suppressPackageStartupMessages(library(getopt))
VERSION <- "2.0.0"
RANKS <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

spec <- matrix(c(
  "seqtab",            "s", 1, "character",
  "taxonomy",          "t", 1, "character",
  "rank",              "r", 1, "character",
  "key",               "k", 1, "character",
  "outdir",            "o", 1, "character",
  "taxon_col",          NA, 1, "character",
  "tax_sep",            NA, 1, "character",
  "prefix",             NA, 1, "character",
  "unclassified",       NA, 1, "character",
  "write_asv",          NA, 0, "logical",
  "help",              "h", 0, "logical",
  "version",           "V", 0, "logical"
), byrow = TRUE, ncol = 4)

usage <- function(status = 0L) {
  cat("Aggregate a DADA2 sequence table at one taxonomic rank.\n\n")
  cat("Usage:\n  count_abundance.R -s seqtab.nochim.rds -t taxonomy.tsv [options]\n\n")
  cat("Required:\n")
  cat("  -s, --seqtab RDS          DADA2 sample-by-sequence matrix\n")
  cat("  -t, --taxonomy TSV        Taxonomy table keyed by sequence\n\n")
  cat("Options:\n")
  cat("  -r, --rank NAME           Kingdom|Phylum|Class|Order|Family|Genus|Species\n")
  cat("                            (default: Genus)\n")
  cat("  -k, --key COLUMN          Sequence-key column; auto-detected by default\n")
  cat("  -o, --outdir DIR          Output directory (default: .)\n")
  cat("      --taxon_col COLUMN    Column containing semicolon taxonomy strings\n")
  cat("      --tax_sep STR         Taxonomy-string separator (default: ;)\n")
  cat("      --prefix STR          Output prefix (default: abundance)\n")
  cat("      --unclassified MODE   lineage|collapse|drop (default: lineage)\n")
  cat("                            lineage keeps parent-specific unknown groups;\n")
  cat("                            collapse merges all as Unclassified; drop removes them\n")
  cat("      --write_asv           Also write a sample-by-ASV table and ASV map\n")
  cat("  -h, --help                Show help\n")
  cat("  -V, --version             Show version\n")
  quit(status = status)
}

opt <- getopt(spec)
if (isTRUE(opt$help)) usage(0L)
if (isTRUE(opt$version)) { cat("count_abundance.R ", VERSION, "\n", sep = ""); quit(status = 0L) }
if (is.null(opt$seqtab) || is.null(opt$taxonomy)) usage(2L)
rank <- if (is.null(opt$rank)) "Genus" else opt$rank
outdir <- if (is.null(opt$outdir)) "." else opt$outdir
prefix <- if (is.null(opt$prefix)) "abundance" else opt$prefix
tax_sep <- if (is.null(opt$tax_sep)) ";" else opt$tax_sep
unclassified <- if (is.null(opt$unclassified)) "lineage" else tolower(opt$unclassified)
if (!rank %in% RANKS) stop("--rank must be one of: ", paste(RANKS, collapse = ", "))
if (!unclassified %in% c("lineage", "collapse", "drop")) stop("--unclassified must be lineage, collapse, or drop")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

seqtab <- readRDS(normalizePath(opt$seqtab, mustWork = TRUE))
if (!is.matrix(seqtab) || is.null(rownames(seqtab)) || is.null(colnames(seqtab))) {
  stop("seqtab must be a named sample-by-sequence matrix")
}
if (any(!is.finite(seqtab)) || any(seqtab < 0)) stop("seqtab contains invalid counts")

taxonomy_path <- normalizePath(opt$taxonomy, mustWork = TRUE)
first_line <- readLines(taxonomy_path, n = 1L, warn = FALSE)
leading_key <- startsWith(first_line, "\t")
taxonomy <- read.delim(taxonomy_path, header = TRUE, check.names = FALSE,
                       stringsAsFactors = FALSE, row.names = if (leading_key) 1L else NULL,
                       quote = "", comment.char = "")

if (!is.null(opt$taxon_col)) {
  if (!opt$taxon_col %in% names(taxonomy)) stop("--taxon_col not found: ", opt$taxon_col)
  split <- strsplit(as.character(taxonomy[[opt$taxon_col]]), tax_sep, fixed = TRUE)
  matrix <- t(vapply(split, function(values) {
    values <- trimws(sub("^[a-zA-Z]__", "", values))
    length(values) <- length(RANKS)
    values
  }, character(length(RANKS))))
  colnames(matrix) <- RANKS
  for (name in RANKS) if (!name %in% names(taxonomy)) taxonomy[[name]] <- matrix[, name]
}

key_candidates <- c("Sequence", "sequence", "ASV", "Feature.ID", "FeatureID", "#OTU ID", "ID")
key <- if (!is.null(opt$key)) opt$key else key_candidates[key_candidates %in% names(taxonomy)][1]
if (is.null(key) || !length(key) || is.na(key)) {
  if (!leading_key) stop("No sequence key found; provide --key")
  taxonomy$Sequence <- rownames(taxonomy)
} else {
  if (!key %in% names(taxonomy)) stop("--key column not found: ", key)
  taxonomy$Sequence <- as.character(taxonomy[[key]])
}
if (anyDuplicated(taxonomy$Sequence)) {
  warning("Duplicate taxonomy keys: retaining the first occurrence")
  taxonomy <- taxonomy[!duplicated(taxonomy$Sequence), , drop = FALSE]
}
for (name in RANKS) if (!name %in% names(taxonomy)) taxonomy[[name]] <- NA_character_

sequence <- colnames(seqtab)
matched <- match(sequence, taxonomy$Sequence)
if (!any(!is.na(matched))) {
  stop("No sequence-table ASVs match the taxonomy keys; use a sequence-keyed taxonomy table")
}
tax <- taxonomy[matched, RANKS, drop = FALSE]
rownames(tax) <- sequence
tax[] <- lapply(tax, function(values) {
  values <- trimws(as.character(values))
  values[!nzchar(values) | is.na(values)] <- NA_character_
  values
})

labels <- tax[[rank]]
missing <- is.na(labels)
if (any(missing) && unclassified == "lineage") {
  rank_index <- match(rank, RANKS)
  parent <- rep("root", nrow(tax))
  if (rank_index > 1L) {
    for (index in seq_len(rank_index - 1L)) {
      known <- !is.na(tax[[RANKS[index]]])
      parent[known] <- tax[[RANKS[index]]][known]
    }
  }
  labels[missing] <- paste0("Unclassified_", gsub("[^A-Za-z0-9_.-]+", "_", parent[missing]))
} else if (any(missing) && unclassified == "collapse") {
  labels[missing] <- "Unclassified"
}

keep <- if (unclassified == "drop") !missing else rep(TRUE, length(labels))
if (!any(keep)) stop("No ASVs remain after unclassified filtering")
taxon_by_sample <- rowsum(t(seqtab[, keep, drop = FALSE]), group = labels[keep], reorder = TRUE)
sample_by_taxon <- t(taxon_by_sample)
totals <- rowSums(sample_by_taxon)
relative <- sample_by_taxon / ifelse(totals > 0, totals, 1)

write_table <- function(value, name, row_label) {
  frame <- data.frame(setNames(list(rownames(value)), row_label), value,
                      check.names = FALSE, row.names = NULL)
  write.table(frame, file.path(outdir, name), sep = "\t", quote = FALSE,
              row.names = FALSE, na = "")
}
write_table(sample_by_taxon, paste0(prefix, "_", rank, "_counts.tsv"), "Sample")
write_table(relative, paste0(prefix, "_", rank, "_relative.tsv"), "Sample")

rank_feature <- data.frame(Taxon = rownames(taxon_by_sample), taxon_by_sample,
                           check.names = FALSE, row.names = NULL)
write.table(rank_feature, file.path(outdir, paste0(prefix, "_", rank, "_counts_by_feature.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

asv_map <- data.frame(ASV = paste0("ASV", seq_along(sequence)), Sequence = sequence,
                      tax, t(seqtab), check.names = FALSE, row.names = NULL)
write.table(asv_map, file.path(outdir, paste0(prefix, "_ASV_counts_with_taxonomy.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "")
if (isTRUE(opt$write_asv)) {
  asv_counts <- seqtab
  colnames(asv_counts) <- asv_map$ASV
  write_table(asv_counts, paste0(prefix, "_ASV_counts.tsv"), "Sample")
  write.table(asv_map[, c("ASV", "Sequence", RANKS)], file.path(outdir, paste0(prefix, "_ASV_taxonomy.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}
message("Wrote abundance tables: ", outdir)
