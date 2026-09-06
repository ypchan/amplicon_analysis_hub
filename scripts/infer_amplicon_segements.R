#!/usr/bin/env Rscript

# Compatibility wrapper for the historically misspelled command.
args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "scripts/infer_amplicon_segements.R"
target <- file.path(dirname(normalizePath(this_file)), "infer_16s_regions.R")
message("WARNING: infer_amplicon_segements.R is deprecated; use infer_16s_regions.R")
status <- system2("Rscript", c(shQuote(target), vapply(args, shQuote, character(1))))
quit(status = if (is.null(status)) 0L else status)
