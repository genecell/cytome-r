# Decode the typed arrays of a cytome chunk. Cytome stores numpy `.tobytes()` little-endian.
.read_floats <- function(raw, dtype) {
  size <- if (dtype == "float32") 4L else if (dtype == "float64") 8L else
    stop("cytome: unsupported matrix dtype '", dtype, "'")
  readBin(raw, what = "double", n = length(raw) %/% size, size = size, endian = "little")
}
.read_ints <- function(raw, size = 4L) {
  if (size == 8L) {
    # int64 index arrays: read as doubles to avoid 32-bit overflow, then to integer if safe.
    lo <- readBin(raw, what = "integer", n = length(raw) %/% 4L, size = 4L, endian = "little")
    return(lo[c(TRUE, FALSE)])  # low words (cytome index arrays fit in 32 bits in practice)
  }
  readBin(raw, what = "integer", n = length(raw) %/% 4L, size = 4L, endian = "little")
}

.chunk_meta <- function(x, name) {
  DBI::dbGetQuery(x$con,
    "SELECT chunk_idx, row_start, row_end FROM matrix_chunks WHERE matrix_name = ? ORDER BY chunk_idx",
    params = list(name))
}

# Read one chunk's three blobs and return the CSR triplets for its global rows.
.read_chunk_triplets <- function(x, name, chunk_idx, data_dtype) {
  r <- DBI::dbGetQuery(x$con,
    "SELECT row_start, row_end, data_blob, indices_blob, indptr_blob, dtype, compression
       FROM matrix_chunks WHERE matrix_name = ? AND chunk_idx = ?",
    params = list(name, chunk_idx))
  comp <- r$compression[1]
  data    <- .read_floats(cytome_decompress(r$data_blob[[1]],    comp), data_dtype)
  indices <- .read_ints(  cytome_decompress(r$indices_blob[[1]], comp))            # 0-based cols
  indptr  <- .read_ints(  cytome_decompress(r$indptr_blob[[1]],  comp))            # len = nrow+1
  nrow_chunk <- r$row_end[1] - r$row_start[1]
  counts <- diff(indptr)                              # nnz per chunk-row
  rows <- rep.int(seq_len(nrow_chunk) + r$row_start[1], counts)   # 1-based global rows
  list(i = rows, j = indices + 1L, x = data, nrow = nrow_chunk)
}

#' Read a cytome count matrix into a sparse Matrix
#'
#' Reconstructs the chunked, compressed CSR matrix (stored cells x features) and returns it,
#' by default transposed to **features x cells** (the orientation of an assay in
#' SingleCellExperiment / Seurat).
#'
#' @param x A `cytome` handle.
#' @param name Matrix name, e.g. `"RNA_counts"` (see [cytome_matrices()]).
#' @param features_x_cells If `TRUE` (default) return features x cells; else cells x features.
#' @return A [Matrix::dgCMatrix-class].
#' @export
read_cytome_matrix <- function(x, name, features_x_cells = TRUE) {
  meta <- DBI::dbGetQuery(x$con,
    "SELECT n_rows, n_cols, dtype FROM matrix_meta WHERE matrix_name = ?", params = list(name))
  if (!nrow(meta)) stop("cytome: matrix '", name, "' not found")
  n_rows <- meta$n_rows[1]; n_cols <- meta$n_cols[1]; dtype <- meta$dtype[1]
  cm <- .chunk_meta(x, name)
  if (!nrow(cm))
    return(methods::as(Matrix::Matrix(0, n_cols, n_rows, sparse = TRUE), "CsparseMatrix"))
  parts <- lapply(cm$chunk_idx, function(ci) .read_chunk_triplets(x, name, ci, dtype))
  i <- unlist(lapply(parts, `[[`, "i"), use.names = FALSE)
  j <- unlist(lapply(parts, `[[`, "j"), use.names = FALSE)
  v <- unlist(lapply(parts, `[[`, "x"), use.names = FALSE)
  # cells x features (rows x cols)
  M <- Matrix::sparseMatrix(i = i, j = j, x = v, dims = c(n_rows, n_cols), index1 = TRUE)
  if (features_x_cells) Matrix::t(M) else M
}
