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
## The stored matrix is cells x features; the seed presents features x cells.
## A row-major matrix is chunked along cells, so a viewport's column range
## selects its chunks; a column-major one (PIASO's ATAC and tiles matrices) is
## chunked along features, so the row range does. Only the chunks overlapping
## the request are read (.extract_block() in read_matrix.R), and the block is
## assembled in time linear in the values read: the previous version grew its
## output with c() once per requested cell, quadratic in the cells of a block.

.cytome_seed_extract <- function(x, index) {
  h <- cytome_open(x@path)
  on.exit(cytome_close(h))
  .extract_block(h, .matrix_info(h, x@matrix_name), cells = index[[2L]], features = index[[1L]])
}

## Seeds saved before the layout slot existed hold row-major matrices.
.seed_layout <- function(x) {
  if (methods::.hasSlot(x, "layout") && length(x@layout)) x@layout else "csr"
}

.register_cytome_delayed <- function(pkg_env) {
  if (!.cytome_delayed_ready()) return(invisible(FALSE))

  methods::setClass("CytomeArraySeed",
    slots = c(path = "character", matrix_name = "character",
              dim = "integer", dimnames = "list", chunk_size = "integer",
              layout = "character"),
    where = pkg_env)

  methods::setMethod("dim", "CytomeArraySeed", function(x) x@dim, where = pkg_env)
  methods::setMethod("dimnames", "CytomeArraySeed", function(x) x@dimnames, where = pkg_env)
  methods::setMethod(S4Arrays::is_sparse, "CytomeArraySeed",
                     function(x) TRUE, where = pkg_env)
  ## Without chunkdim() the default block grid straddles storage chunks, so a
  ## single block pulls two chunks and every chunk is read twice. A row-major
  ## matrix is chunked along cells (seed columns), a column-major one along
  ## features (seed rows).
  methods::setMethod(DelayedArray::chunkdim, "CytomeArraySeed",
                     ## clamped: a cytome written with chunk_size 2000 but
                     ## holding 12 cells would otherwise report a chunk larger
                     ## than the array, which DelayedArray rejects.
                     function(x) {
                       if (.seed_layout(x) == "csc") c(min(x@chunk_size, x@dim[1L]), x@dim[2L])
                       else c(x@dim[1L], min(x@chunk_size, x@dim[2L]))
                     },
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
  info <- .matrix_info(x, matrix_name)
  .check_chunks_present(info, .chunk_index(x, info))
  ## The extent of a storage chunk along its axis: cells for a row-major
  ## matrix, features for a column-major one. as.integer(NULL) is integer(0),
  ## which made chunkdim() return a length-1 vector and DelayedArray reject the
  ## seed, so a missing size falls back to the whole axis.
  extent <- if (info$layout == "csc") info$n_cols else info$n_rows
  chunk_size <- info$chunk_size
  if (!length(chunk_size) || is.na(chunk_size) || chunk_size < 1L) chunk_size <- as.integer(extent)
  var <- cytome_var(x, matrix_name)
  obs <- cytome_obs(x, "cells")
  ## Stored cells x features; the seed presents features x cells.
  seed <- methods::new("CytomeArraySeed",
                       path = normalizePath(path), matrix_name = matrix_name,
                       dim = c(info$n_cols, info$n_rows),
                       dimnames = list(.feature_ids(var), .cell_names(x, obs)),
                       chunk_size = chunk_size, layout = info$layout)
  DelayedArray::DelayedArray(seed)
}
