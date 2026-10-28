#' Query ATAC fragments on a chromosome, optionally within an interval
#'
#' Fragments are stored in one of two forms, and both are read:
#' * compressed chunks (`fragment_chunks`), written by PIASO's importer and the Python
#'   package's streaming import, and the only form `subset()` and `filter_cells()` keep;
#' * one row per fragment in `fragments_<chrom>` tables (columns
#'   `start, end_, cell_idx, dup_count`), written by the older Python import.
#'
#' A chunk holds about half a million fragments sorted by start. An interval query decodes
#' only the chunks that start before its end, one at a time, as the Python reader does.
#'
#' @param x A `cytome` handle.
#' @param chrom Chromosome, e.g. `"chr1"`.
#' @param start,end Optional 0-based half-open interval `[start, end)`. Give both, or neither
#'   for every fragment on `chrom`.
#' @param with_barcode If `TRUE` (default), join `cell_idx` to the cell `barcode`.
#' @return A data.frame: chrom, start, end, cell_idx, dup_count and, with `with_barcode`,
#'   barcode, ordered by start, end and cell. `dup_count` is `NA` for fragments stored in
#'   chunks, which do not keep it.
#' @export
cytome_fragments <- function(x, chrom, start = NULL, end = NULL, with_barcode = TRUE) {
  if (xor(is.null(start), is.null(end)))
    stop("cytome: give both start and end, or neither")
  if (!is.null(start) && !(is.numeric(start) && is.numeric(end) && length(start) == 1L &&
                           length(end) == 1L && !is.na(start) && !is.na(end) && end >= start))
    stop("cytome: start and end must be two numbers with start <= end")
  tbl <- paste0("fragments_", chrom)
  qt <- DBI::dbQuoteIdentifier(x$con, tbl)
  has_rows <- .has_table(x, tbl) &&
    nrow(DBI::dbGetQuery(x$con, sprintf("SELECT 1 FROM %s LIMIT 1", qt))) > 0
  has_chunks <- .has_table(x, "fragment_chunks") &&
    nrow(DBI::dbGetQuery(x$con, "SELECT 1 FROM fragment_chunks WHERE chrom = ? LIMIT 1",
                         params = list(chrom))) > 0
  if (has_rows) {
    q <- sprintf("SELECT start, end_ AS \"end\", cell_idx, dup_count FROM %s", qt)
    df <- if (is.null(start)) DBI::dbGetQuery(x$con, q) else
      DBI::dbGetQuery(x$con, paste0(q, " WHERE end_ > ? AND start < ?"),   # half-open overlap
                      params = list(start, end))
  } else if (has_chunks) {
    df <- .fragments_from_chunks(x, chrom, start, end)
  } else {
    if (!.has_table(x, tbl)) {
      have <- .fragment_chroms(x)
      stop("cytome: no fragments for '", chrom, "'",
           if (length(have)) paste0("; this store has fragments on ", paste(have, collapse = ", ")))
    }
    df <- data.frame(start = integer(0), end = integer(0), cell_idx = integer(0), dup_count = integer(0))
  }
  df <- data.frame(chrom = rep(chrom, nrow(df)), start = df$start, end = df$end,
                   cell_idx = df$cell_idx, dup_count = df$dup_count, stringsAsFactors = FALSE)
  df <- df[order(df$start, df$end, df$cell_idx), , drop = FALSE]
  rownames(df) <- NULL
  if (with_barcode && .has_table(x, "cells")) {
    bc <- DBI::dbGetQuery(x$con, "SELECT cell_idx, barcode FROM cells")
    m <- match(df$cell_idx, bc$cell_idx)
    if (anyNA(m))
      warning("cytome: ", sum(is.na(m)), " fragments on ", chrom,
              " name cells that are not in the cells table; their barcode is NA")
    df$barcode <- as.character(bc$barcode)[m]
  }
  df
}

# Fragments of `chrom` from the compressed chunks, decoded one chunk at a time. Encoding 0
# stores absolute starts and ends; encoding 1 the starts as differences from the previous one
# and the ends as lengths. A chunk that starts at or after the interval's end is skipped.
.fragments_from_chunks <- function(x, chrom, start = NULL, end = NULL) {
  cols <- .table_columns(x, "fragment_chunks")
  where <- "chrom = ?"
  params <- list(chrom)
  if (!is.null(end) && "min_start" %in% cols) {
    where <- paste(where, "AND (min_start IS NULL OR min_start < ?)")
    params <- c(params, list(end))
  }
  ids <- DBI::dbGetQuery(x$con, paste("SELECT chunk_idx FROM fragment_chunks WHERE", where,
                                      "ORDER BY chunk_idx"), params = params)$chunk_idx
  sql <- paste0("SELECT starts_blob, ends_blob, cell_idx_blob, compression, ",
                if ("encoding" %in% cols) "COALESCE(encoding, 0)" else "0",
                " AS encoding FROM fragment_chunks WHERE chrom = ? AND chunk_idx = ?")
  parts <- lapply(ids, function(k) {
    r <- DBI::dbGetQuery(x$con, sql, params = list(chrom, k))
    f <- cytome_decode_fragments(r$starts_blob[[1]], r$ends_blob[[1]], r$cell_idx_blob[[1]],
                                 .comp(r$compression[1]), as.integer(r$encoding[1]),
                                 sprintf("fragments %s chunk %d", chrom, as.integer(k)))
    if (!is.null(start)) {
      keep <- f$start < end & f$end > start
      f <- lapply(f, `[`, keep)
    }
    f
  })
  s <- as.integer(unlist(lapply(parts, `[[`, "start"), use.names = FALSE))
  data.frame(start = s,
             end = as.integer(unlist(lapply(parts, `[[`, "end"), use.names = FALSE)),
             cell_idx = as.integer(unlist(lapply(parts, `[[`, "cell_idx"), use.names = FALSE)),
             dup_count = rep(NA_integer_, length(s)))
}

# Chromosomes that hold fragments, in either form (for messages).
.fragment_chroms <- function(x) {
  out <- character(0)
  if (.has_table(x, "fragment_chunks"))
    out <- DBI::dbGetQuery(x$con, "SELECT DISTINCT chrom FROM fragment_chunks")$chrom
  if (.has_table(x, "fragment_meta"))
    out <- c(out, DBI::dbGetQuery(x$con,
      "SELECT chrom FROM fragment_meta WHERE n_fragments > 0")$chrom)
  sort(unique(as.character(out)))
}
