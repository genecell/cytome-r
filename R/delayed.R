## Out-of-core backend: a DelayedArray seed that pulls from the .cytome on
## demand. Once the assay is a DelayedArray the rest of Bioconductor operates
## on it unmodified, which is the difference between a converter and a backend.
##
## DelayedArray / S4Arrays / SparseArray are in Suggests, not Imports: someone
## who only wants the converter should not have to install the Bioconductor
## array stack. That means the class and its methods CANNOT be declared at
## build time -- they are registered in .onLoad(), and only if the packages are
## actually present. A NAMESPACE importFrom() on a Suggests package is exactly
## the R CMD check note this avoids.

.cytome_delayed_ready <- function() {
  all(vapply(c("DelayedArray", "S4Arrays", "SparseArray"),
             requireNamespace, logical(1), quietly = TRUE))
}

## Chunk-aligned extraction.
##
## The stored matrix is cells x features and is chunked along CELLS; the seed
## presents features x cells, so a viewport's COLUMN range is what selects
## chunks. Read only the chunks that overlap it -- the previous version read
## the whole matrix and then subset, which made peak memory the full matrix
## plus the block, i.e. strictly worse than loading it once. That is the
## difference between a lazy wrapper and an out-of-core backend.
##
## Chunk k covers 1-based stored rows (row_start + 1):row_end.

.cytome_seed_extract <- function(x, index) {
  i <- index[[1L]]                                   # features (stored columns)
  j <- index[[2L]]                                   # cells    (stored rows)
  n_features <- x@dim[1L]; n_cells <- x@dim[2L]
  want <- if (is.null(j)) seq_len(n_cells) else as.integer(j)

  h <- cytome_open(x@path)
  on.exit(cytome_close(h))
  dtype <- DBI::dbGetQuery(h$con,
    "SELECT dtype FROM matrix_meta WHERE matrix_name = ?",
    params = list(x@matrix_name))$dtype[1]
  cm <- .chunk_meta(h, x@matrix_name)

  empty <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                                dims = c(n_features, length(want)))
  if (!nrow(cm) || !length(want)) return(.subset_features(empty, i))

  ## Only the chunks whose row span intersects the requested cells. Pulled out
  ## so the selectivity itself can be asserted in a test -- "we read fewer
  ## chunks" is the whole claim, and it is not observable from the values.
  hit <- .cytome_chunks_for(cm, want)
  if (!length(hit)) return(.subset_features(empty, i))

  parts <- lapply(hit, function(ci)
    .read_chunk_triplets(h, x@matrix_name, ci, dtype))
  gi <- unlist(lapply(parts, `[[`, "i"), use.names = FALSE)   # global cell, 1-based
  gj <- unlist(lapply(parts, `[[`, "j"), use.names = FALSE)   # feature, 1-based
  gx <- unlist(lapply(parts, `[[`, "x"), use.names = FALSE)

  ## Map global cell -> position in `want`. A cell can be asked for more than
  ## once (DelayedArray permits repeated subscripts), so build the output by
  ## expanding rather than by a single match().
  keep <- gi %in% want
  gi <- gi[keep]; gj <- gj[keep]; gx <- gx[keep]
  idx_by_cell <- split(seq_along(gi), gi)
  out_i <- integer(0); out_j <- integer(0); out_x <- numeric(0)
  for (n in seq_along(want)) {
    sel <- idx_by_cell[[as.character(want[n])]]
    if (is.null(sel)) next
    out_i <- c(out_i, gj[sel])                       # feature -> seed row
    out_j <- c(out_j, rep.int(n, length(sel)))       # position -> seed column
    out_x <- c(out_x, gx[sel])
  }
  M <- Matrix::sparseMatrix(i = out_i, j = out_j, x = out_x,
                            dims = c(n_features, length(want)))
  .subset_features(M, i)
}

.subset_features <- function(M, i) if (is.null(i)) M else M[i, , drop = FALSE]

## Chunk k covers 1-based stored rows (row_start + 1):row_end. Returns the
## chunk_idx values whose span contains at least one requested cell.
.cytome_chunks_for <- function(cm, want) {
  if (!nrow(cm) || !length(want)) return(integer(0))
  cm$chunk_idx[vapply(seq_len(nrow(cm)), function(k)
    any(want > cm$row_start[k] & want <= cm$row_end[k]), logical(1))]
}

.register_cytome_delayed <- function(pkg_env) {
  if (!.cytome_delayed_ready()) return(invisible(FALSE))

  methods::setClass("CytomeArraySeed",
    slots = c(path = "character", matrix_name = "character",
              dim = "integer", dimnames = "list", chunk_size = "integer"),
    where = pkg_env)

  methods::setMethod("dim", "CytomeArraySeed", function(x) x@dim, where = pkg_env)
  methods::setMethod("dimnames", "CytomeArraySeed", function(x) x@dimnames, where = pkg_env)
  methods::setMethod(S4Arrays::is_sparse, "CytomeArraySeed",
                     function(x) TRUE, where = pkg_env)
  ## Without chunkdim() the default block grid straddles storage chunks, so a
  ## single block pulls two chunks and every chunk is read twice.
  methods::setMethod(DelayedArray::chunkdim, "CytomeArraySeed",
                     ## clamped: a cytome written with chunk_size 2000 but
                     ## holding 12 cells would otherwise report a chunk larger
                     ## than the array, which DelayedArray rejects.
                     function(x) c(x@dim[1L], min(x@chunk_size, x@dim[2L])),
                     where = pkg_env)
  methods::setMethod(S4Arrays::extract_array, "CytomeArraySeed",
                     function(x, index) as.matrix(.cytome_seed_extract(x, index)),
                     where = pkg_env)
  methods::setMethod(SparseArray::extract_sparse_array, "CytomeArraySeed",
                     function(x, index)
                       methods::as(.cytome_seed_extract(x, index), "SVT_SparseArray"),
                     where = pkg_env)

  ## DelayedArray is mid-transition between two sparse-block generics: the
  ## block machinery in 0.28 still calls OLD_extract_sparse_array() (returning
  ## a SparseArraySeed) and only newer versions call extract_sparse_array()
  ## (returning an SVT_SparseArray). Declaring is_sparse() TRUE while
  ## implementing only the new one means colSums() dies inside
  ## read_sparse_block(). Register whichever the installed version has.
  if ("OLD_extract_sparse_array" %in% getNamespaceExports("DelayedArray")) {
    methods::setMethod(DelayedArray::OLD_extract_sparse_array, "CytomeArraySeed",
                       function(x, index)
                         methods::as(.cytome_seed_extract(x, index), "SparseArraySeed"),
                       where = pkg_env)
  }
  invisible(TRUE)
}

.onLoad <- function(libname, pkgname) {
  .register_cytome_delayed(asNamespace(pkgname))
}

#' A `.cytome` matrix as a DelayedArray
#'
#' Returns a [DelayedArray::DelayedArray] backed by the file: nothing is read
#' until something asks for values, and block-processing code
#' (`DelayedMatrixStats`, `scuttle`, `scran`, ...) then works on a matrix that
#' never has to be held whole.
#'
#' Extraction is **chunk-aligned**: only the storage chunks overlapping the
#' requested cells are read, and `chunkdim()` reports the storage geometry so
#' DelayedArray's block grid lands on chunk boundaries rather than straddling
#' them. Peak memory is therefore set by the block size, not by the size of
#' the matrix.
#'
#' For data that fits comfortably in memory, `delayed = FALSE` is faster --
#' the reason to use this is data that does not fit, and code you did not
#' write (`scran`, `scuttle`, `DelayedMatrixStats`) that must run on it.
#'
#' Requires the suggested `DelayedArray`, `S4Arrays` and `SparseArray`
#' packages.
#'
#' @param path Path to a `.cytome` file.
#' @param matrix_name Matrix to wrap, e.g. `"RNA_counts"`.
#' @return A `DelayedArray`, features x cells.
#' @seealso [read_cytome()] with `delayed = TRUE`, [cytome_stream()]
#' @examples
#' ref <- system.file("extdata", "reference.cytome", package = "cytome")
#' if (nzchar(ref) && requireNamespace("DelayedArray", quietly = TRUE)) {
#'   m <- cytome_delayed(ref, "RNA_counts")
#'   dim(m)
#' }
#' @export
cytome_delayed <- function(path, matrix_name) {
  if (!.cytome_delayed_ready())
    stop("cytome: install 'DelayedArray', 'S4Arrays' and 'SparseArray' to use ",
         "cytome_delayed()")
  x <- cytome_open(path); on.exit(cytome_close(x))
  meta <- cytome_matrices(x)
  if (!matrix_name %in% meta$matrix_name)
    stop("cytome: matrix '", matrix_name, "' not found; have: ",
         paste(meta$matrix_name, collapse = ", "))
  row <- meta[meta$matrix_name == matrix_name, , drop = FALSE]
  ## cytome_matrices() does not select chunk_size, and as.integer(NULL) is
  ## integer(0) -- which made chunkdim() return a length-1 vector and
  ## DelayedArray reject the seed. Ask matrix_meta directly.
  chunk_size <- as.integer(DBI::dbGetQuery(x$con,
    "SELECT chunk_size FROM matrix_meta WHERE matrix_name = ?",
    params = list(matrix_name))$chunk_size[1])
  if (is.na(chunk_size) || !length(chunk_size) || chunk_size < 1L)
    chunk_size <- as.integer(row$n_rows[1])
  var <- cytome_var(x, matrix_name)
  obs <- cytome_obs(x, "cells")
  ## Stored cells x features; the seed presents features x cells.
  seed <- methods::new("CytomeArraySeed",
                       path = normalizePath(path), matrix_name = matrix_name,
                       dim = c(as.integer(row$n_cols[1]), as.integer(row$n_rows[1])),
                       dimnames = list(.feature_ids(var), .cell_ids(obs)),
                       chunk_size = chunk_size)
  DelayedArray::DelayedArray(seed)
}
