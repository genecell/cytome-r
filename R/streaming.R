#' Stream a cytome matrix chunk-by-chunk (cytome-native, out-of-core)
#'
#' Iterates the matrix's stored chunks, decoding **one chunk at a time** into a sparse
#' sub-matrix and passing it to `FUN` — so you can compute over a matrix far larger than RAM
#' without ever materialising the whole thing, directly from the cytome format (no BPCells or
#' other backend required). This is the foundation for R-side streaming single-cell analysis.
#'
#' @param x A `cytome` handle.
#' @param name Matrix name, e.g. `"RNA_counts"`.
#' @param FUN A function called per chunk as `FUN(mat, row_start, row_end, chunk_idx, ...)`,
#'   where `mat` is the chunk's sparse sub-matrix.
#' @param features_x_cells If `TRUE` each chunk is features x cells; else cells x features (default).
#' @param ... Extra arguments passed to `FUN`.
#' @return A list of the per-chunk return values of `FUN` (reduce as you wish).
#' @examples
#' \dontrun{
#' x <- cytome_open("data.cytome")
#' # per-cell total counts without loading the whole matrix:
#' totals <- cytome_stream(x, "RNA_counts", function(m, s, e, i) Matrix::rowSums(m))
#' totals <- unlist(totals)
#' }
#' @export
cytome_stream <- function(x, name, FUN, features_x_cells = FALSE, ...) {
  meta <- DBI::dbGetQuery(x$con,
    "SELECT n_cols, dtype FROM matrix_meta WHERE matrix_name = ?", params = list(name))
  if (!nrow(meta)) stop("cytome: matrix '", name, "' not found")
  n_cols <- meta$n_cols[1]; dtype <- meta$dtype[1]
  cm <- .chunk_meta(x, name)
  out <- vector("list", nrow(cm))
  for (k in seq_len(nrow(cm))) {
    tp <- .read_chunk_triplets(x, name, cm$chunk_idx[k], dtype)
    local_i <- tp$i - cm$row_start[k]                 # rows local to this chunk (1-based)
    sub <- Matrix::sparseMatrix(i = local_i, j = tp$j, x = tp$x,
                                dims = c(tp$nrow, n_cols), index1 = TRUE)   # cells_chunk x features
    if (features_x_cells) sub <- Matrix::t(sub)
    out[[k]] <- FUN(sub, cm$row_start[k], cm$row_end[k], cm$chunk_idx[k], ...)
  }
  out
}
