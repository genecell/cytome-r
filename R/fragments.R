#' Query ATAC fragments overlapping a genomic interval
#'
#' Fragments are stored per chromosome in `fragments_<chrom>` tables (columns
#' `start, end_, cell_idx, dup_count`) with an SQLite R*-tree index, so interval queries
#' are fast and native — no Python.
#'
#' @param x A `cytome` handle.
#' @param chrom Chromosome, e.g. `"chr1"`.
#' @param start,end Optional 0-based half-open interval `[start, end)`. If both `NULL`, returns
#'   all fragments on `chrom`.
#' @param with_barcode If `TRUE` (default), join `cell_idx` to the cell `barcode`.
#' @return A data.frame: chrom, start, end, cell_idx, [barcode], dup_count.
#' @export
cytome_fragments <- function(x, chrom, start = NULL, end = NULL, with_barcode = TRUE) {
  tbl <- paste0("fragments_", chrom)
  if (!.has_table(x, tbl)) stop("cytome: no fragment table for '", chrom, "'")
  q <- sprintf("SELECT start, end_, cell_idx, dup_count FROM %s", DBI::dbQuoteIdentifier(x$con, tbl))
  params <- list()
  if (!is.null(start) && !is.null(end)) {
    q <- paste0(q, " WHERE end_ > ? AND start < ?")     # half-open overlap
    params <- list(start, end)
  }
  df <- DBI::dbGetQuery(x$con, q, params = params)
  if (!nrow(df)) return(data.frame(chrom = character(0), start = integer(0), end = integer(0),
                                   cell_idx = integer(0), dup_count = integer(0)))
  names(df)[names(df) == "end_"] <- "end"
  df <- cbind(chrom = chrom, df)
  if (with_barcode && .has_table(x, "cells")) {
    bc <- DBI::dbGetQuery(x$con, "SELECT cell_idx, barcode FROM cells")
    df$barcode <- bc$barcode[match(df$cell_idx, bc$cell_idx)]
  }
  df
}
