#' Graphs stored in a cytome
#'
#' Lists the cell-by-cell graphs of a store, in either form the Python package
#' writes: compressed row chunks (`graph_meta` / `graph_chunks`, cytome 0.3.6
#' and later), and one row per edge (`graph_edges`, earlier versions and this
#' package's own writer).
#'
#' @param x A `cytome` handle.
#' @return A character vector of graph names.
#' @export
cytome_graphs <- function(x) {
  names <- character(0)
  if (.has_table(x, "graph_meta"))
    names <- DBI::dbGetQuery(x$con,
      "SELECT graph_name FROM graph_meta WHERE axis = 'obs'")$graph_name
  if (.has_table(x, "graph_edges"))
    names <- c(names, DBI::dbGetQuery(x$con,
      "SELECT DISTINCT graph_name FROM graph_edges")$graph_name)
  sort(unique(as.character(names)))
}

#' Read one graph as a sparse matrix
#'
#' @param x A `cytome` handle.
#' @param name Graph name, from [cytome_graphs()].
#' @param n Square dimension; defaults to the stored shape, or the cell count
#'   for a graph stored one row per edge.
#' @return A [Matrix::dgCMatrix-class], cells x cells.
#' @export
cytome_graph <- function(x, name, n = NULL) {
  meta <- .graph_meta(x, name)
  if (!is.null(meta)) return(.read_graph_chunks(x, name, meta, n))
  if (!.has_table(x, "graph_edges"))
    stop("cytome: graph '", name, "' not found; have: ",
         paste(cytome_graphs(x), collapse = ", "))
  e <- DBI::dbGetQuery(x$con,
    "SELECT row_idx, col_idx, value FROM graph_edges WHERE graph_name = ?",
    params = list(name))
  if (!nrow(e)) stop("cytome: graph '", name, "' not found; have: ",
                     paste(cytome_graphs(x), collapse = ", "))
  if (is.null(n)) n <- nrow(cytome_obs(x, "cells"))
  ## stored 0-based, as Python writes them
  Matrix::sparseMatrix(i = e$row_idx + 1L, j = e$col_idx + 1L, x = e$value,
                       dims = c(n, n))
}

## The chunked form's shape, entry count and value type, or NULL when the
## graph is not stored that way (or the store predates the tables).
.graph_meta <- function(x, name) {
  if (!.has_table(x, "graph_meta")) return(NULL)
  m <- DBI::dbGetQuery(x$con,
    "SELECT n_rows, n_cols, n_nonzero, dtype, n_chunks FROM graph_meta
      WHERE graph_name = ? AND axis = 'obs'", params = list(name))
  if (!nrow(m)) NULL else m
}

## A graph from its row chunks: each holds a run of rows as CSR, data in the
## graph's float type, indices int32 and a chunk-local indptr int64, each blob
## compressed as the chunk's `compression` says (lz4 by default).
.read_graph_chunks <- function(x, name, meta, n = NULL) {
  n_rows <- as.integer(meta$n_rows[1]); n_cols <- as.integer(meta$n_cols[1])
  if (!is.null(n)) {
    if (n < n_rows || n < n_cols)
      stop("cytome: graph '", name, "' is ", n_rows, " x ", n_cols,
           ", larger than n = ", n)
    n_rows <- n_cols <- as.integer(n)
  }
  ch <- DBI::dbGetQuery(x$con,
    "SELECT row_start, row_end, data_blob, indices_blob, indptr_blob, compression
       FROM graph_chunks WHERE graph_name = ? AND axis = 'obs' ORDER BY chunk_idx",
    params = list(name))
  dtype <- meta$dtype[1]
  parts <- lapply(seq_len(nrow(ch)), function(k) {
    comp <- ch$compression[k]
    ptr <- .read_ints(cytome_decompress(ch$indptr_blob[[k]], comp), 8L)
    rows <- seq_len(ch$row_end[k] - ch$row_start[k]) + ch$row_start[k]   # 1-based
    list(i = rep.int(rows, diff(ptr)),
         j = .read_ints(cytome_decompress(ch$indices_blob[[k]], comp)) + 1L,
         x = .read_floats(cytome_decompress(ch$data_blob[[k]], comp), dtype))
  })
  i <- unlist(lapply(parts, `[[`, "i"), use.names = FALSE)
  j <- unlist(lapply(parts, `[[`, "j"), use.names = FALSE)
  v <- unlist(lapply(parts, `[[`, "x"), use.names = FALSE)
  if (length(v) != meta$n_nonzero[1])
    stop("cytome: graph '", name, "': ", length(v), " entries in its chunks, ",
         meta$n_nonzero[1], " recorded")
  Matrix::sparseMatrix(i = as.integer(i), j = as.integer(j), x = as.numeric(v),
                       dims = c(n_rows, n_cols))
}
