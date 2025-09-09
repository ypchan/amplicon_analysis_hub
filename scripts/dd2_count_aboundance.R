#!/usr/bin/env Rscript
# Build abundance tables from DADA2 seqtab.nochim.rds and taxonomy.tsv
# - CLI parsing via getopt
# - Outputs:
#   1) <prefix>_ASV_counts_with_taxonomy.tsv            (features rows: ASV+Sequence+taxonomy+samples)
#   2) <prefix>_<Rank>_counts_with_taxonomy.tsv         (features rows: Taxon+taxonomy+samples)
#   3) <prefix>_<Rank>_counts.tsv / <prefix>_<Rank>_relative.tsv  (samples rows)
# Comments in English only

suppressPackageStartupMessages({
  library(getopt)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
})

# ---------- CLI (getopt) ----------
spec <- matrix(c(
  "seqtab",      "s", 1, "character", "Path to seqtab.nochim.rds (required)",
  "taxonomy",    "t", 1, "character", "Path to taxonomy.tsv (required)",
  "rank",        "r", 1, "character", "Aggregation rank: Kingdom|Phylum|Class|Order|Family|Genus|Species [Genus]",
  "key",         "k", 1, "character", "Column in taxonomy holding the sequence key (e.g. Sequence/ASV/Feature.ID)",
  "outdir",      "o", 1, "character", "Output directory [.]",
  "taxon-col",    NA, 1, "character", "If taxonomy has a single 'Taxon' column, provide its name",
  "tax-sep",      NA, 1, "character", "Separator used in the Taxon column [';']",
  "prefix",       NA, 1, "character", "Output filename prefix [abundance]",
  "write-asv",    NA, 0, "logical",   "Also write samples×ASV wide table (optional)",
  "help",        "h", 0, "logical",   "Show help and exit"
), byrow = TRUE, ncol = 5)

opt <- getopt(spec)

usage_text <- function() {
  cat("Usage:\n  Rscript make_abundance_getopt.R --seqtab seqtab.nochim.rds --taxonomy taxonomy.tsv [options]\n\nOptions:\n")
  apply(spec, 1, function(row) {
    cat(sprintf("  --%-12s %s  %s\n",
                row[1],
                ifelse(is.na(row[2]), "   ", paste0("(-", row[2], ")")),
                row[5]))
  })
  cat("\nExamples:\n",
      "  dd2_count_aboundace.R -s seqtab.nochim.rds -t taxonomy.tsv -r Genus -o dd2_count_aboundace\n",
      "  dd2_count_aboundace.R -s seqtab.nochim.rds -t taxonomy.tsv -r Species -o out --taxon-col Taxon --tax-sep ';'\n",
      "  dd2_count_aboundace.R -s seqtab.nochim.rds -t taxonomy.tsv -r Phylum -o out --key Feature.ID\n",
      sep = "")
}

if (isTRUE(opt$help) || is.null(opt$seqtab) || is.null(opt$taxonomy)) {
  usage_text(); quit(status = if (isTRUE(opt$help)) 0 else 2)
}

# Defaults
if (is.null(opt$rank))      opt$rank   <- "Genus"
if (is.null(opt$outdir))    opt$outdir <- "."
if (is.null(opt$`tax-sep`)) opt$`tax-sep` <- ";"
if (is.null(opt$prefix))    opt$prefix <- "abundance"
if (!dir.exists(opt$outdir)) dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

ranks_allowed <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
if (!(opt$rank %in% ranks_allowed)) {
  stop("--rank must be one of: ", paste(ranks_allowed, collapse = ", "))
}

message("[1/6] Loading seqtab: ", opt$seqtab)
seqtab <- readRDS(opt$seqtab)
if (!is.matrix(seqtab)) stop("seqtab.nochim.rds must be a matrix")
if (is.null(rownames(seqtab)) || is.null(colnames(seqtab)))
  stop("seqtab must have rownames (samples) and colnames (ASV sequences)")

message("[2/6] Loading taxonomy: ", opt$taxonomy)
tax <- read_tsv(opt$taxonomy, show_col_types = FALSE, progress = FALSE, comment = "")
tax <- as.data.frame(tax, stringsAsFactors = FALSE, check.names = FALSE)

# ---------- Helpers ----------
split_taxon_column <- function(df, col, sep = ";") {
  # Split a single 'Taxon' column into 7 ranks (K..S). Supports "k__Bacteria; p__Firmicutes; ..."
  if (!col %in% names(df)) return(df)
  x <- df[[col]]
  parts <- str_split(x, pattern = fixed(sep), n = 7, simplify = TRUE)
  parts <- trimws(parts)
  strip_pref <- function(s) sub("^[a-z]__", "", s, perl = TRUE)  # remove k__/p__/...
  parts <- apply(parts, 2, strip_pref)
  colnames(parts) <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")[seq_len(ncol(parts))]
  bind_cols(df, as.data.frame(parts, stringsAsFactors = FALSE))
}

detect_key <- function(df) {
  candidates <- c("Sequence","sequence","Seq","seq","ASV","Feature.ID","FeatureID","feature","Feature",
                  "#OTU ID","OTU","ID")
  for (c in candidates) if (c %in% names(df)) return(c)
  NA_character_
}

fill_na_rank <- function(x, label = "Unclassified") ifelse(is.na(x) | x == "", label, x)

pick1 <- function(x) {
  # Prefer the first non-missing and non-"Unclassified" value
  ux <- unique(x)
  ux <- ux[!(is.na(ux) | ux == "")]
  pr <- ux[ux != "Unclassified"]
  if (length(pr)) pr[1] else if (length(ux)) ux[1] else NA_character_
}

# ---------- Prepare taxonomy ----------
if (!is.null(opt$`taxon-col`)) {
  tax <- split_taxon_column(tax, opt$`taxon-col`, sep = opt$`tax-sep`)
}

key_col <- if (!is.null(opt$key)) opt$key else detect_key(tax)
if (is.na(key_col)) {
  if (!is.null(rownames(tax)) && all(nchar(rownames(tax)) > 0)) {
    tax$Sequence <- rownames(tax); key_col <- "Sequence"
  } else {
    stop("Cannot find sequence key column in taxonomy.tsv; supply one via --key")
  }
}
if (!key_col %in% names(tax)) stop("Key column '", key_col, "' not found in taxonomy.tsv")
names(tax)[names(tax) == key_col] <- "Sequence"

# Ensure all rank columns exist
for (rk in ranks_allowed) if (!rk %in% names(tax)) tax[[rk]] <- NA_character_

# Deduplicate taxonomy by Sequence if needed
if (any(duplicated(tax$Sequence))) {
  warning("Duplicate Sequence keys in taxonomy; keeping the first for each Sequence.")
  tax <- tax |> group_by(Sequence) |> summarise(across(everything(), pick1), .groups = "drop")
}

# ---------- Build ASV map and long table ----------
message("[3/6] Building ASV map and long table")
seqs   <- colnames(seqtab)
asv_id <- paste0("ASV", seq_along(seqs))
map_asv <- tibble(ASV = asv_id, Sequence = seqs)

counts_long <- as.data.frame(seqtab, check.names = FALSE) |>
  rownames_to_column("Sample") |>
  as_tibble() |>
  pivot_longer(-Sample, names_to = "Sequence", values_to = "Count") |>
  filter(Count > 0L)

asv_long <- counts_long |>
  left_join(map_asv, by = "Sequence") |>
  relocate(ASV, .after = Sequence) |>
  left_join(tax, by = "Sequence")

for (rk in ranks_allowed) if (rk %in% names(asv_long)) asv_long[[rk]] <- fill_na_rank(asv_long[[rk]])

# ---------- ASV-level merged table (features as rows) ----------
message("[4/6] Writing ASV-level merged table with taxonomy")
asv_counts_feature <-
  as.data.frame(seqtab, check.names = FALSE) |>
  t() |>
  as.data.frame(check.names = FALSE) |>
  tibble::rownames_to_column("Sequence") |>
  left_join(tibble(Sequence = colnames(seqtab), ASV = paste0("ASV", seq_along(colnames(seqtab)))), by = "Sequence") |>
  left_join(tax, by = "Sequence") |>
  relocate(ASV, Sequence, Kingdom, Phylum, Class, Order, Family, Genus, Species)

write_tsv(asv_counts_feature, file.path(opt$outdir, paste0(opt$prefix, "_ASV_counts_with_taxonomy.tsv")))

# Optionally also write samples×ASV wide table (samples rows) without taxonomy
if (isTRUE(opt$`write-asv`)) {
  asv_wide <- asv_long |>
    select(Sample, ASV, Count) |>
    pivot_wider(names_from = ASV, values_from = Count, values_fill = 0L) |>
    arrange(Sample)
  write_tsv(asv_wide, file.path(opt$outdir, paste0(opt$prefix, "_ASV_counts.tsv")))
  write_tsv(map_asv |> left_join(tax, by = "Sequence"),
            file.path(opt$outdir, "ASV_taxonomy.tsv"))
}

# ---------- Aggregate to chosen rank ----------
message("[5/6] Aggregating at rank: ", opt$rank)
rank <- opt$rank
label_vec <- if (rank == "Species" && "Genus" %in% names(asv_long)) {
  starts_with_genus <- grepl(paste0("^\\Q", asv_long$Genus, "\\E\\b"), asv_long$Species, perl = TRUE)
  ifelse(starts_with_genus, asv_long$Species, paste(asv_long$Genus, asv_long$Species))
} else {
  asv_long[[rank]]
}

tab <- asv_long |>
  mutate(Taxon = label_vec) |>
  group_by(Sample, Taxon) |>
  summarise(Abundance = sum(Count), .groups = "drop")

# Samples×Taxon wide (counts)
wide_counts <- tab |>
  pivot_wider(names_from = Taxon, values_from = Abundance, values_fill = 0L) |>
  arrange(Sample)

# Samples×Taxon wide (relative)
wide_rel <- tab |>
  group_by(Sample) |>
  mutate(RelAbund = Abundance / sum(Abundance)) |>
  ungroup() |>
  select(Sample, Taxon, RelAbund) |>
  pivot_wider(names_from = Taxon, values_from = RelAbund, values_fill = 0) |>
  arrange(Sample)

write_tsv(wide_counts, file.path(opt$outdir, paste0(opt$prefix, "_", rank, "_counts.tsv")))
write_tsv(wide_rel,    file.path(opt$outdir, paste0(opt$prefix, "_", rank, "_relative.tsv")))

# ---------- Rank-level merged table (features as rows) ----------
message("[6/6] Writing rank-level merged table with taxonomy")
tax_by_rank <- asv_long |>
  mutate(Taxon = label_vec) |>
  select(Taxon, Kingdom, Phylum, Class, Order, Family, Genus, Species) |>
  group_by(Taxon) |>
  summarise(
    Kingdom = pick1(Kingdom),
    Phylum  = pick1(Phylum),
    Class   = pick1(Class),
    Order   = pick1(Order),
    Family  = pick1(Family),
    Genus   = pick1(Genus),
    Species = pick1(Species),
    .groups = "drop"
  )

rank_counts_feature <- tab |>
  pivot_wider(names_from = Sample, values_from = Abundance, values_fill = 0L) |>
  left_join(tax_by_rank, by = "Taxon") |>
  relocate(Taxon, Kingdom, Phylum, Class, Order, Family, Genus, Species)

write_tsv(rank_counts_feature,
          file.path(opt$outdir, paste0(opt$prefix, "_", rank, "_counts_with_taxonomy.tsv")))

cat("Done.\n")
