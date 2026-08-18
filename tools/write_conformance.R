#!/usr/bin/env Rscript
## Write a .cytome from R for the reverse conformance check: read the Python
## fixture into a Seurat object, write it back out, hand the result to
## tools/reverse_conformance.py.
args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1] else "/tmp/from_r.cytome"
suppressMessages(library(cytome))
ref <- system.file("extdata", "reference.cytome", package = "cytome")
if (!nzchar(ref)) stop("reference.cytome not installed")
so <- read_cytome(ref, as = "Seurat")
cat("read:", paste(names(so@assays), collapse = ", "), "\n")
write_cytome(so, out)
cat("wrote:", out, "\n")
