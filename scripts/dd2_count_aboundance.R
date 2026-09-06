#!/usr/bin/env Rscript

# Compatibility wrapper for the historically misspelled command.
args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "scripts/dd2_count_aboundance.R"
target <- file.path(dirname(normalizePath(this_file)), "count_abundance.R")
message("WARNING: dd2_count_aboundance.R is deprecated; use count_abundance.R")
status <- system2("Rscript", c(shQuote(target), vapply(args, shQuote, character(1))))
quit(status = if (is.null(status)) 0L else status)
