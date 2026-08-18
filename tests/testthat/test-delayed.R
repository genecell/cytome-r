# The DelayedArray backend. Point of the tests: a delayed assay must give the
# SAME numbers as the in-memory one -- an out-of-core path that is subtly
# different is worse than none.

test_that("cytome_delayed matches the in-memory matrix", {
  skip_if_not_installed("DelayedArray")
  skip_if_not_installed("SparseArray")
  ref <- system.file("extdata", "reference.cytome", package = "cytome")
  skip_if(ref == "")

  d <- cytome_delayed(ref, "RNA_counts")
  expect_s4_class(d, "DelayedArray")

  x <- cytome_open(ref); on.exit(cytome_close(x))
  M <- read_cytome_matrix(x, "RNA_counts")
  expect_equal(dim(d), dim(M))
  ## read_cytome_matrix() returns the bare matrix; the DelayedArray carries
  ## feature and cell names, so compare values unnamed and the names separately.
  expect_equal(unname(as.matrix(d)), unname(as.matrix(M)), tolerance = 1e-6)
  expect_equal(rownames(d), cytome:::.feature_ids(cytome_var(x, "RNA_counts")))
  expect_equal(colnames(d), cytome:::.cell_ids(cytome_obs(x, "cells")))
  expect_equal(unname(as.numeric(DelayedArray::colSums(d))),
               unname(as.numeric(Matrix::colSums(M))), tolerance = 1e-6)
  ## subsetting must go through the seed's index, not a materialised copy
  expect_equal(unname(as.matrix(d[1:3, 1:4])),
               unname(as.matrix(M[1:3, 1:4])), tolerance = 1e-6)
})

test_that("read_cytome(delayed = TRUE) gives an SCE whose assay is on disk", {
  skip_if_not_installed("DelayedArray")
  skip_if_not_installed("SingleCellExperiment")
  ref <- system.file("extdata", "reference.cytome", package = "cytome")
  skip_if(ref == "")

  sce_d <- read_cytome(ref, delayed = TRUE)
  sce_m <- read_cytome(ref, delayed = FALSE)
  expect_s4_class(SummarizedExperiment::assay(sce_d), "DelayedArray")
  expect_equal(dim(sce_d), dim(sce_m))
  expect_equal(as.matrix(SummarizedExperiment::assay(sce_d)),
               as.matrix(SummarizedExperiment::assay(sce_m)), tolerance = 1e-6)
})

# The point of the chunk-aligned seed is that it reads FEWER chunks. That is
# not observable from the values, so assert it directly on the selection.

test_that("only the chunks overlapping the request are selected", {
  cm <- data.frame(chunk_idx = 0:3,
                   row_start = c(0L, 10L, 20L, 30L),
                   row_end   = c(10L, 20L, 30L, 40L))
  expect_equal(cytome:::.cytome_chunks_for(cm, 1:10), 0L)          # first chunk only
  expect_equal(cytome:::.cytome_chunks_for(cm, 11L), 1L)           # boundary, low side
  expect_equal(cytome:::.cytome_chunks_for(cm, 20L), 1L)           # boundary, high side
  expect_equal(cytome:::.cytome_chunks_for(cm, 21L), 2L)
  expect_equal(cytome:::.cytome_chunks_for(cm, c(5L, 35L)), c(0L, 3L))   # skips 1 and 2
  expect_equal(cytome:::.cytome_chunks_for(cm, 1:40), 0:3)
  expect_equal(cytome:::.cytome_chunks_for(cm, integer(0)), integer(0))
})

test_that("a multi-chunk matrix reads back identically through the seed", {
  skip_if_not_installed("DelayedArray")
  skip_if_not_installed("SingleCellExperiment")

  set.seed(1)
  m <- Matrix::rsparsematrix(20, 57, density = 0.4)      # 20 features x 57 cells
  dimnames(m) <- list(paste0("g", 1:20), paste0("c", 1:57))
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = m))
  p <- tempfile(fileext = ".cytome")
  write_cytome(sce, p, chunk_size = 10L)                 # -> 6 chunks

  x <- cytome_open(p); on.exit(cytome_close(x))
  expect_gt(nrow(cytome:::.chunk_meta(x, "RNA_counts")), 1L)   # really multi-chunk

  d <- cytome_delayed(p, "RNA_counts")
  M <- read_cytome_matrix(x, "RNA_counts")
  expect_equal(unname(as.matrix(d)), unname(as.matrix(M)), tolerance = 1e-6)

  ## a slice wholly inside one chunk, a slice spanning chunks, and a
  ## reordered/repeated one
  for (j in list(3:8, 8:15, c(50L, 2L, 2L, 33L))) {
    expect_equal(unname(as.matrix(d[, j, drop = FALSE])),
                 unname(as.matrix(M[, j, drop = FALSE])), tolerance = 1e-6)
  }
  expect_equal(unname(as.numeric(DelayedArray::colSums(d))),
               unname(as.numeric(Matrix::colSums(M))), tolerance = 1e-6)
})
