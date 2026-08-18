# The API surface itself: one write generic, one read entry point.

test_that("write_cytome dispatches on the object, not on the function name", {
  skip_if_not_installed("SingleCellExperiment")
  m <- Matrix::rsparsematrix(6, 4, density = 0.5)
  dimnames(m) <- list(paste0("g", 1:6), paste0("c", 1:4))
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = m))

  p <- tempfile(fileext = ".cytome")
  write_cytome(sce, p)                       # no _sce suffix needed
  expect_true(file.exists(p))
  back <- read_cytome(p)
  expect_s4_class(back, "SingleCellExperiment")
  expect_equal(dim(back), dim(sce))
})

test_that("write_cytome refuses a class it has no method for, by name", {
  expect_error(write_cytome(1:10, tempfile()), "no write_cytome\\(\\) method")
})

test_that("read_cytome validates 'as' instead of failing later", {
  ref <- system.file("extdata", "reference.cytome", package = "cytome")
  skip_if(ref == "")
  expect_error(read_cytome(ref, as = "seurat"), "should be one of")
})

test_that("as = 'cytome' returns the open handle", {
  ref <- system.file("extdata", "reference.cytome", package = "cytome")
  skip_if(ref == "")
  x <- read_cytome(ref, as = "cytome")
  expect_s3_class(x, "cytome")
  expect_true("RNA_counts" %in% cytome_matrices(x)$matrix_name)
  cytome_close(x)
})

test_that("delayed = TRUE is refused for Seurat rather than silently ignored", {
  ref <- system.file("extdata", "reference.cytome", package = "cytome")
  skip_if(ref == "")
  expect_error(read_cytome(ref, as = "Seurat", delayed = TRUE),
               "only implemented for")
})
