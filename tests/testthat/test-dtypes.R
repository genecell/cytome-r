# Values are decoded in the dtype each chunk records, indices as int32 or int64 by their
# length. Only float32 and float64 used to decode: an int32 store (from_10x_h5,
# from_cellranger) failed entirely, an int64 embedding (a Visium spatial one) vanished
# from read_cytome(), and the chunks of a merge of int32 and float32 stores could not be
# read.

dec <- function(bytes, dtype, ...) cytome:::cytome_decode_values(bytes, "none", dtype, ...)

test_that("every numpy dtype decodes, in either byte order", {
  expect_equal(dec(writeBin(c(-1L, 2L, 300L), raw(), size = 4, endian = "little"), "int32"), c(-1, 2, 300))
  expect_equal(dec(writeBin(c(-1L, 2L, 300L), raw(), size = 4, endian = "big"), ">i4"), c(-1, 2, 300))
  expect_equal(dec(writeBin(c(-7L, 9L), raw(), size = 2, endian = "little"), "int16"), c(-7, 9))
  expect_equal(dec(as.raw(c(0xff, 0x01)), "int8"), c(-1, 1))
  expect_equal(dec(as.raw(c(0xff, 0x01)), "uint8"), c(255, 1))
  expect_equal(dec(as.raw(c(0xff, 0xff)), "uint16"), 65535)
  expect_equal(dec(as.raw(c(0, 1, 2)), "bool"), c(0, 1, 1))
  expect_equal(dec(writeBin(c(0.5, -2), raw(), size = 4, endian = "little"), "<f4"), c(0.5, -2))
  expect_equal(dec(writeBin(c(0.1, -2), raw(), size = 8, endian = "big"), ">f8"), c(0.1, -2))
  # int64: low word first; 2^40 + 5 and -3
  i64 <- c(as.raw(c(5, 0, 0, 0, 0, 1, 0, 0)), as.raw(c(0xfd, rep(0xff, 7))))
  expect_equal(dec(i64, "int64"), c(2^40 + 5, -3))
  # float16: 1, -2, 65504, Inf
  f16 <- as.raw(c(0x00, 0x3c, 0x00, 0xc0, 0xff, 0x7b, 0x00, 0x7c))
  expect_equal(dec(f16, "float16"), c(1, -2, 65504, Inf))
  expect_error(dec(raw(3), "int32"), "whole number")
  expect_error(dec(raw(8), "int32", n = 3), "recorded")
  expect_error(dec(raw(4), "complex64"), "unsupported dtype")
  expect_error(dec(as.raw(c(rep(0, 7), 0x40)), "int64"), "2\\^53")      # 2^62
})

test_that("int32 matrices and an int64 embedding written by from_cellranger read", {
  p <- fixture_store("py_10x"); e <- fixture_expected("py_10x")
  x <- cytome_open(p); on.exit(cytome_close(x))
  for (nm in c("RNA_counts", "ATAC_counts")) {
    expect_equal(cytome_matrices(x)$dtype[cytome_matrices(x)$matrix_name == nm], "int32")
    M <- read_cytome_matrix(x, nm, features_x_cells = FALSE)
    expect_equal(unname(as.matrix(M)), as.matrix(expected_matrix(e$matrices[[nm]])), info = nm)
  }
  E <- cytome_embedding(x, "RNA_spatial")
  expect_equal(unname(E), expected_dense(e$embeddings$RNA_spatial))
  expect_equal(rownames(E), e$barcodes)
})

test_that("read_cytome builds objects from an int32 store", {
  skip_if_not_installed("SingleCellExperiment")
  p <- fixture_store("py_10x"); e <- fixture_expected("py_10x")
  sce <- read_cytome(p)
  expect_equal(dim(sce), c(4L, e$n_cells))
  expect_true("ATAC" %in% SingleCellExperiment::altExpNames(sce))
  expect_true("RNA_spatial" %in% SingleCellExperiment::reducedDimNames(sce))
  expect_equal(unname(as.matrix(SummarizedExperiment::assay(sce))),
               t(as.matrix(expected_matrix(e$matrices$RNA_counts))))
})

test_that("each chunk is decoded in its own dtype", {
  # merge() of an int32 store and a float32 store: one chunk of each, matrix_meta says int32
  p <- fixture_store("py_merge_mixed"); e <- fixture_expected("py_merge_mixed")
  expect_setequal(e$chunk_dtypes, c("float32", "int32"))
  x <- cytome_open(p); on.exit(cytome_close(x))
  M <- read_cytome_matrix(x, "RNA_counts", features_x_cells = FALSE)
  expect_equal(unname(as.matrix(M)), as.matrix(expected_matrix(e$matrices$RNA_counts)))
})

# A store whose index blobs are int64 (scipy keeps int64 indices through slicing, and the
# Python writer stored them as they were): the length tells the width.
.widen_indices <- function(path) {
  con <- DBI::dbConnect(RSQLite::SQLite(), path)
  on.exit(DBI::dbDisconnect(con))
  ch <- DBI::dbGetQuery(con, "SELECT id, indices_blob, indptr_blob, compression FROM matrix_chunks")
  for (k in seq_len(nrow(ch))) {
    w <- function(b) {
      v <- readBin(cytome:::cytome_decompress(b, ch$compression[k]), "integer", n = length(b) * 4, size = 4)
      cytome:::cytome_compress(writeBin(as.vector(rbind(v, 0L)), raw(), size = 4), ch$compression[k])
    }
    DBI::dbExecute(con, "UPDATE matrix_chunks SET indices_blob = ?, indptr_blob = ? WHERE id = ?",
                   params = list(list(w(ch$indices_blob[[k]])), list(w(ch$indptr_blob[[k]])), ch$id[k]))
  }
}

test_that("int64 index blobs are read by their length", {
  p <- fixture_copy("py_10x"); e <- fixture_expected("py_10x")
  .widen_indices(p)
  x <- cytome_open(p); on.exit(cytome_close(x))
  M <- read_cytome_matrix(x, "RNA_counts", features_x_cells = FALSE)
  expect_equal(unname(as.matrix(M)), as.matrix(expected_matrix(e$matrices$RNA_counts)))
})

test_that("a chunk whose arrays disagree is an error, not a matrix", {
  p <- fixture_copy("py_10x")
  con <- DBI::dbConnect(RSQLite::SQLite(), p)
  b <- DBI::dbGetQuery(con, "SELECT data_blob, compression FROM matrix_chunks WHERE matrix_name = 'RNA_counts' AND chunk_idx = 0")
  short <- cytome:::cytome_compress(utils::head(cytome:::cytome_decompress(b$data_blob[[1]], b$compression), -4), b$compression)
  DBI::dbExecute(con, "UPDATE matrix_chunks SET data_blob = ? WHERE matrix_name = 'RNA_counts' AND chunk_idx = 0",
                 params = list(list(short)))
  DBI::dbDisconnect(con)
  x <- cytome_open(p); on.exit(cytome_close(x))
  expect_error(read_cytome_matrix(x, "RNA_counts"), "RNA_counts row chunk 0: .* bytes of data")
})
