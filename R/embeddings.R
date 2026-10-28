#' List embeddings stored in a cytome
#' @param x A `cytome` handle.
#' @return Character vector of embedding/array names (from `embedding_meta`).
#' @export
cytome_embeddings <- function(x) {
  if (!.has_table(x, "embedding_meta")) return(character(0))
  DBI::dbGetQuery(x$con, "SELECT array_name FROM embedding_meta")$array_name
}

#' Read a dense embedding (e.g. X_umap, X_svd) into a matrix
#'
#' Embeddings are stored as dense, row-major, compressed chunks (`dense_chunks`), in any
#' numeric dtype (a Visium `spatial` embedding is int64).
#' @param x A `cytome` handle.
#' @param name Embedding name, e.g. `"X_umap"` (see [cytome_embeddings()]).
#' @return A numeric matrix (rows x n_cols) with the entity's identifiers (cell barcodes for
#'   a cell embedding) as row names when they are available.
#' @export
cytome_embedding <- function(x, name) {
  cols <- .table_columns(x, "embedding_meta")
  meta <- DBI::dbGetQuery(x$con, paste0(
    "SELECT n_rows, n_cols, dtype", if ("entity" %in% cols) ", entity" else "",
    " FROM embedding_meta WHERE array_name = ?"), params = list(name))
  if (!nrow(meta)) stop("cytome: embedding '", name, "' not found")
  n_rows <- as.integer(meta$n_rows[1]); n_cols <- as.integer(meta$n_cols[1]); dtype <- meta$dtype[1]
  ch <- DBI::dbGetQuery(x$con,
    "SELECT chunk_idx, row_start, row_end, data_blob, dtype, compression
       FROM dense_chunks WHERE array_name = ? ORDER BY chunk_idx", params = list(name))
  out <- matrix(0.0, nrow = n_rows, ncol = n_cols)
  covered <- 0
  for (k in seq_len(nrow(ch))) {
    r0 <- as.integer(ch$row_start[k]); r1 <- as.integer(ch$row_end[k]); nr <- r1 - r0
    if (r0 < 0L || r1 > n_rows || nr < 0L)
      stop("cytome: embedding '", name, "' chunk ", ch$chunk_idx[k], " covers rows ", r0, "-", r1,
           " of ", n_rows)
    vals <- cytome_decode_values(ch$data_blob[[k]], .comp(ch$compression[k]),
                                 .chunk_dtype(ch$dtype[k], dtype), as.numeric(nr) * n_cols,
                                 sprintf("embedding %s chunk %d", name, as.integer(ch$chunk_idx[k])))
    if (nr) out[(r0 + 1L):r1, ] <- matrix(vals, nrow = nr, ncol = n_cols, byrow = TRUE)   # row-major
    covered <- covered + nr
  }
  if (covered != n_rows)
    warning("cytome: embedding '", name, "' has chunks for ", covered, " of its ", n_rows, " rows")
  entity <- if ("entity" %in% names(meta)) meta$entity[1] else "cells"
  ids <- tryCatch(.entity_ids(x, if (is.na(entity)) "cells" else entity), error = function(e) NULL)
  if (length(ids) == n_rows) rownames(out) <- ids
  out
}
