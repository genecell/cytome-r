# Matrix storage.
#
# A matrix is stored cells x features, in chunks, in one of two layouts, decided as the
# Python package decides it (MeasurementLayer.is_column_major):
#   - row-major: CSR chunks of consecutive cells in `matrix_chunks`, the layout of the Python
#     writers;
#   - column-major: CSC chunks of consecutive features in `matrix_csc_chunks`, every cell in
#     each (`has_csc = 1` and no row chunks), the default of PIASO's importer and quantifier
#     for ATAC peaks and tiles.
# A row-major matrix may also carry column chunks (a feature index); its row chunks are read.
# `matrix_meta.layout`, where a writer recorded it, must name one of the two.
#
# A chunk's arrays are decoded in C++ (src/decode.cpp): the values by the chunk's own dtype
# column, since a merge of an int32 and a float32 store writes both into one matrix, and the
# indices as int32 or int64 by the blob's length. Each array is checked against the extent
# the chunk records, so a corrupt chunk is an error, not a matrix with values out of place.

# The chunk's own dtype, or the matrix's when the chunk does not record one.
.chunk_dtype <- function(chunk_dtype, matrix_dtype) {
  if (length(chunk_dtype) && !is.na(chunk_dtype) && nzchar(chunk_dtype)) chunk_dtype else matrix_dtype
}

# A missing compression label means the bytes are stored as they are, as in Python.
.comp <- function(method) {
  if (!length(method) || is.na(method)) "none" else method
}

# What matrix_meta says about `name`, with its layout decided.
.matrix_info <- function(x, name) {
  want <- c("n_rows", "n_cols", "n_nonzero", "dtype", "row_entity", "col_entity", "chunk_size",
            "n_chunks", "has_csc", "csc_chunk_size", "layout")
  have <- intersect(want, .table_columns(x, "matrix_meta"))
  m <- DBI::dbGetQuery(x$con, sprintf("SELECT %s FROM matrix_meta WHERE matrix_name = ?",
                                      paste(have, collapse = ", ")), params = list(name))
  if (!nrow(m)) stop("cytome: matrix '", name, "' not found")
  get <- function(k) if (k %in% names(m)) m[[k]][1] else NA
  label <- get("layout")
  if (!is.na(label) && !label %in% c("csr", "csc"))
    stop("cytome: matrix '", name, "' is stored in layout '", label,
         "', which this version of cytome cannot read")
  has_csc <- isTRUE(as.numeric(get("has_csc")) == 1)
  n_chunks <- as.numeric(get("n_chunks"))
  col_major <- has_csc && (is.na(n_chunks) || n_chunks == 0)
  list(name = name, n_rows = as.integer(get("n_rows")), n_cols = as.integer(get("n_cols")),
       n_nonzero = as.numeric(get("n_nonzero")), dtype = get("dtype"),
       row_entity = get("row_entity"), col_entity = get("col_entity"),
       layout = if (col_major) "csc" else "csr", has_csc = has_csc,
       chunk_size = as.integer(if (col_major) get("csc_chunk_size") else get("chunk_size")))
}

.layout_table <- function(layout) {
  if (layout == "csc") list(table = "matrix_csc_chunks", start = "col_start", end = "col_end", axis = "column")
  else list(table = "matrix_chunks", start = "row_start", end = "row_end", axis = "row")
}

# The chunks of a matrix in `layout` (its own by default): chunk_idx, start and end (0-based,
# end exclusive), in order along the axis, checked to lie within it without overlapping.
.chunk_index <- function(x, info, layout = info$layout) {
  lt <- .layout_table(layout)
  if (!.has_table(x, lt$table))
    return(data.frame(chunk_idx = integer(0), start = integer(0), end = integer(0)))
  ch <- DBI::dbGetQuery(x$con, sprintf(
    "SELECT chunk_idx, %s AS start, %s AS \"end\" FROM %s WHERE matrix_name = ? ORDER BY %s, %s",
    lt$start, lt$end, lt$table, lt$start, lt$end), params = list(info$name))
  extent <- if (layout == "csc") info$n_cols else info$n_rows
  bad <- ch$start < 0 | ch$end > extent | ch$end < ch$start
  if (any(bad))
    stop("cytome: matrix '", info$name, "': ", lt$axis, " chunk ", ch$chunk_idx[which(bad)[1]],
         " covers ", ch$start[which(bad)[1]], "-", ch$end[which(bad)[1]], " of ", extent)
  n <- nrow(ch)
  if (n > 1 && any(ch$start[-1] < ch$end[-n]))
    stop("cytome: matrix '", info$name, "': its ", lt$axis, " chunks overlap")
  ch
}

# A matrix that records nonzeros must have chunks to hold them.
.check_chunks_present <- function(info, ch, layout = info$layout) {
  if (!nrow(ch) && !is.na(info$n_nonzero) && info$n_nonzero > 0)
    stop("cytome: matrix '", info$name, "' records ", format(info$n_nonzero, big.mark = ","),
         " nonzeros but has no ", .layout_table(layout)$axis, " chunks")
}

# Decode chunks `ids` of a matrix in `layout`: a list of list(start, end, p, i, x), with `p`
# the chunk-local pointer over its rows (CSR) or columns (CSC) and `i` the 0-based feature
# (CSR) or cell (CSC) of each value.
.decode_chunks <- function(x, info, ids, layout = info$layout) {
  lt <- .layout_table(layout)
  minor <- if (layout == "csc") info$n_rows else info$n_cols
  sql <- sprintf(paste("SELECT chunk_idx, %s AS start, %s AS \"end\", data_blob, indices_blob,",
                       "indptr_blob, dtype, compression FROM %s WHERE matrix_name = ? AND chunk_idx = ?"),
                 lt$start, lt$end, lt$table)
  lapply(ids, function(ci) {
    r <- DBI::dbGetQuery(x$con, sql, params = list(info$name, ci))
    d <- cytome_decode_sparse(r$data_blob[[1]], r$indices_blob[[1]], r$indptr_blob[[1]],
                              .comp(r$compression[1]), .chunk_dtype(r$dtype[1], info$dtype),
                              as.integer(r$end[1] - r$start[1]), minor,
                              sprintf("%s %s chunk %d", info$name, lt$axis, as.integer(ci)))
    d$start <- as.integer(r$start[1])
    d$end <- as.integer(r$end[1])
    d
  })
}

# The decoded chunks as triplets in stored coordinates (1-based): `major` is the cell (CSR) or
# feature (CSC) of each value, `minor` the other.
.triplets <- function(parts) {
  # typed even when no chunk was read: DelayedArray asks for an empty block to learn the type
  list(major = as.integer(unlist(lapply(parts, function(d)
         rep.int(seq.int(d$start + 1L, length.out = d$end - d$start), diff(d$p))), use.names = FALSE)),
       minor = as.integer(unlist(lapply(parts, `[[`, "i"), use.names = FALSE)) + 1L,
       x = as.numeric(unlist(lapply(parts, `[[`, "x"), use.names = FALSE)))
}

.check_nnz <- function(info, nnz, whole = TRUE) {
  if (nnz > .Machine$integer.max)
    stop("cytome: matrix '", info$name, "' has ", format(nnz, big.mark = ","),
         " nonzeros, past what one R sparse matrix holds (2^31 - 1); stream it with ",
         "cytome_stream() or back it with cytome_delayed()")
  if (whole && !is.na(info$n_nonzero) && nnz != info$n_nonzero)
    warning("cytome: matrix '", info$name, "': matrix_meta records ",
            format(info$n_nonzero, big.mark = ","), " nonzeros but its chunks hold ",
            format(nnz, big.mark = ","), "; the chunks are used")
}

#' Read a cytome count matrix into a sparse Matrix
#'
#' Reconstructs the chunked, compressed matrix (stored cells x features) and returns it, by
#' default transposed to **features x cells** (the orientation of an assay in
#' SingleCellExperiment / Seurat). Both storage layouts are read: row-major chunks of cells
#' (the Python writers) and column-major chunks of features (PIASO's importer and quantifier,
#' for ATAC peaks and tiles).
#'
#' @param x A `cytome` handle.
#' @param name Matrix name, e.g. `"RNA_counts"` (see [cytome_matrices()]).
#' @param features_x_cells If `TRUE` (default) return features x cells; else cells x features.
#' @return A [Matrix::dgCMatrix-class].
#' @export
read_cytome_matrix <- function(x, name, features_x_cells = TRUE) {
  info <- .matrix_info(x, name)
  ch <- .chunk_index(x, info)
  .check_chunks_present(info, ch)
  t3 <- .triplets(.decode_chunks(x, info, ch$chunk_idx))
  .check_nnz(info, length(t3$x))
  # columns = the stored major axis: CSR gives features x cells, CSC cells x features
  csr <- info$layout == "csr"
  M <- Matrix::sparseMatrix(i = t3$minor, j = t3$major, x = t3$x,
                            dims = if (csr) c(info$n_cols, info$n_rows) else c(info$n_rows, info$n_cols))
  if (csr == features_x_cells) M else Matrix::t(M)
}

# The features x cells block of `cells` and `features` (1-based stored indices, in any order,
# repeats allowed; NULL = all), reading only the chunks that hold them: chunks of cells for a
# row-major matrix, of features for a column-major one. Linear in the values read.
.extract_block <- function(x, info, cells = NULL, features = NULL) {
  uc <- if (is.null(cells)) NULL else sort(unique(as.integer(cells)))
  uf <- if (is.null(features)) NULL else sort(unique(as.integer(features)))
  ch <- .chunk_index(x, info)
  .check_chunks_present(info, ch)
  csr <- info$layout == "csr"
  along <- if (csr) uc else uf
  ids <- if (is.null(along)) ch$chunk_idx else .cytome_chunks_for(ch, along)
  t3 <- .triplets(.decode_chunks(x, info, ids))
  cell <- if (csr) t3$major else t3$minor
  feat <- if (csr) t3$minor else t3$major
  ci <- if (is.null(uc)) cell else match(cell, uc)
  fi <- if (is.null(uf)) feat else match(feat, uf)
  keep <- !is.na(ci) & !is.na(fi)
  .check_nnz(info, sum(keep), whole = FALSE)
  M <- Matrix::sparseMatrix(i = fi[keep], j = ci[keep], x = t3$x[keep],
                            dims = c(if (is.null(uf)) info$n_cols else length(uf),
                                     if (is.null(uc)) info$n_rows else length(uc)))
  if (!is.null(cells) && !identical(as.integer(cells), uc))
    M <- M[, match(as.integer(cells), uc), drop = FALSE]
  if (!is.null(features) && !identical(as.integer(features), uf))
    M <- M[match(as.integer(features), uf), , drop = FALSE]
  M
}

## The chunks (chunk_idx) whose span holds at least one of `want`, 1-based stored indices
## along the chunked axis; `cm` is a chunk index (chunk_idx, start, end; 0-based, end
## exclusive, ordered by start). By interval search, not a scan of every chunk per index.
.cytome_chunks_for <- function(cm, want) {
  if (!nrow(cm) || !length(want)) return(integer(0))
  w0 <- as.numeric(want) - 1
  k <- findInterval(w0, cm$start)          # the last chunk that starts at or before w0
  ok <- k > 0
  k <- k[ok]
  k <- k[w0[ok] < cm$end[k]]
  cm$chunk_idx[sort(unique(k))]
}
