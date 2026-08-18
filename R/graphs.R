#' Graphs stored in a cytome
#'
#' Reads the `graph_edges` table -- the same table the Python package writes,
#' so a KNN graph built on either side is readable from the other.
#'
#' @param x A `cytome` handle.
#' @return A character vector of graph names.
#' @export
cytome_graphs <- function(x) {
  if (!.has_table(x, "graph_edges")) return(character(0))
  as.character(DBI::dbGetQuery(x$con,
    "SELECT DISTINCT graph_name FROM graph_edges ORDER BY graph_name")$graph_name)
}

#' Read one graph as a sparse matrix
#'
#' @param x A `cytome` handle.
#' @param name Graph name, from [cytome_graphs()].
#' @param n Square dimension; defaults to the cell count.
#' @return A [Matrix::dgCMatrix-class], cells x cells.
#' @export
cytome_graph <- function(x, name, n = NULL) {
  if (!.has_table(x, "graph_edges"))
    stop("cytome: this file has no graph_edges table")
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
