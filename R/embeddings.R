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
#' Embeddings are stored as dense, row-major, compressed chunks (`dense_chunks`).
#' @param x A `cytome` handle.
#' @param name Embedding name, e.g. `"X_umap"` (see [cytome_embeddings()]).
#' @return A numeric matrix (cells x n_cols) with cell barcodes as row names when available.
#' @export
cytome_embedding <- function(x, name) {
  meta <- DBI::dbGetQuery(x$con,
    "SELECT n_rows, n_cols, dtype FROM embedding_meta WHERE array_name = ?", params = list(name))
  if (!nrow(meta)) stop("cytome: embedding '", name, "' not found")
  n_rows <- meta$n_rows[1]; n_cols <- meta$n_cols[1]; dtype <- meta$dtype[1]
  ch <- DBI::dbGetQuery(x$con,
    "SELECT chunk_idx, row_start, row_end, data_blob, dtype, compression
       FROM dense_chunks WHERE array_name = ? ORDER BY chunk_idx", params = list(name))
  out <- matrix(0.0, nrow = n_rows, ncol = n_cols)
  for (k in seq_len(nrow(ch))) {
    vals <- .read_floats(cytome_decompress(ch$data_blob[[k]], ch$compression[k]), dtype)
    nr <- ch$row_end[k] - ch$row_start[k]
    # row-major: fill by row (byrow = TRUE)
    out[(ch$row_start[k] + 1L):ch$row_end[k], ] <- matrix(vals, nrow = nr, ncol = n_cols, byrow = TRUE)
  }
  out
}
