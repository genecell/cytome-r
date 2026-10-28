#' Stream a cytome matrix chunk-by-chunk (cytome-native, out-of-core)
#'
#' Iterates the matrix's stored chunks, decoding **one chunk at a time** into a sparse
#' sub-matrix and passing it to `FUN` — so you can compute over a matrix far larger than RAM
#' without ever materialising the whole thing, directly from the cytome format (no BPCells or
#' other backend required). This is the foundation for R-side streaming single-cell analysis.
#'
#' A matrix is stored either by cells (row chunks, the layout of the Python writers) or by
#' features (column chunks, the default of PIASO's importer and quantifier for ATAC peaks and
#' tiles). `along` says which chunks to stream: `"cells"` passes blocks of cells with every
#' feature, `"features"` blocks of features with every cell. A matrix streams only along the
#' chunks it has; asking for the other axis is an error that says so, never an empty result.
#'
#' @param x A `cytome` handle.
#' @param name Matrix name, e.g. `"RNA_counts"`.
#' @param FUN A function called per chunk as `FUN(mat, start, end, chunk_idx, ...)`, where
#'   `mat` is the chunk's sparse sub-matrix and `start`, `end` its 0-based, end-exclusive
#'   range of cells (`along = "cells"`) or of features (`along = "features"`).
#' @param features_x_cells If `TRUE` each chunk is features x cells; else cells x features (default).
#' @param ... Extra arguments passed to `FUN`.
#' @param along `"cells"` (default) or `"features"`: the axis the stored chunks run along.
#' @return A list of the per-chunk return values of `FUN` (reduce as you wish).
#' @examples
#' ref <- system.file("extdata", "reference.cytome", package = "cytome")
#' if (nzchar(ref)) {
#'   x <- cytome_open(ref)
#'   # per-cell total counts without loading the whole matrix:
#'   totals <- unlist(cytome_stream(x, "RNA_counts", function(m, s, e, i) Matrix::rowSums(m)))
#'   cytome_close(x)
#' }
#' @export
cytome_stream <- function(x, name, FUN, features_x_cells = FALSE, ...,
                          along = c("cells", "features")) {
  along <- match.arg(along)
  info <- .matrix_info(x, name)
  layout <- if (along == "cells") "csr" else "csc"
  if (along == "cells" && info$layout == "csc")
    stop("cytome: '", name, "' is stored by features (column chunks), so it has no chunks of ",
         "cells to stream; use along = \"features\", or read it whole with read_cytome_matrix()")
  if (along == "features" && !info$has_csc)
    stop("cytome: '", name, "' is stored by cells (row chunks) and has no column chunks; ",
         "use along = \"cells\"")
  ch <- .chunk_index(x, info, layout)
  .check_chunks_present(info, ch, layout)
  out <- vector("list", nrow(ch))
  for (k in seq_len(nrow(ch))) {
    d <- .decode_chunks(x, info, ch$chunk_idx[k], layout)[[1L]]
    n <- d$end - d$start
    local <- rep.int(seq_len(n), diff(d$p))
    sub <- if (layout == "csr")
      Matrix::sparseMatrix(i = local, j = d$i + 1L, x = d$x, dims = c(n, info$n_cols))   # cells_chunk x features
    else
      Matrix::sparseMatrix(i = d$i + 1L, j = local, x = d$x, dims = c(info$n_rows, n))   # cells x features_chunk
    if (features_x_cells) sub <- Matrix::t(sub)
    out[[k]] <- FUN(sub, d$start, d$end, ch$chunk_idx[k], ...)
  }
  out
}
