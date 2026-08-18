# Alternative assays used to be dropped silently: a Seurat object with RNA and
# ATAC round-tripped as RNA only, and the file gave no sign half of it was
# missing. These assert the second modality survives.

.fixture <- function() system.file("extdata", "reference.cytome", package = "cytome")

test_that("a Seurat object's extra assays survive the round trip", {
  skip_if_not_installed("SeuratObject")
  ref <- .fixture(); skip_if(ref == "")

  so <- read_cytome(ref, as = "Seurat")
  expect_true(all(c("RNA", "ATAC") %in% names(so@assays)))

  p <- tempfile(fileext = ".cytome")
  write_cytome(so, p)

  x <- cytome_open(p); on.exit(cytome_close(x))
  mats <- cytome_matrices(x)$matrix_name
  expect_true("RNA_counts" %in% mats)
  expect_true("ATAC_counts" %in% mats)          # the regression

  back <- read_cytome(p, as = "Seurat")
  expect_true(all(c("RNA", "ATAC") %in% names(back@assays)))
  a1 <- SeuratObject::LayerData(so, assay = "ATAC", layer = "counts")
  a2 <- SeuratObject::LayerData(back, assay = "ATAC", layer = "counts")
  expect_equal(dim(a1), dim(a2))
  expect_equal(as.numeric(as.matrix(a1)), as.numeric(as.matrix(a2)))
})

test_that("an SCE's altExps survive the round trip", {
  skip_if_not_installed("SingleCellExperiment")
  ref <- .fixture(); skip_if(ref == "")

  sce <- read_cytome(ref)
  expect_true("ATAC" %in% SingleCellExperiment::altExpNames(sce))

  p <- tempfile(fileext = ".cytome")
  write_cytome(sce, p)
  back <- read_cytome(p)
  expect_true("ATAC" %in% SingleCellExperiment::altExpNames(back))
  expect_equal(dim(SingleCellExperiment::altExp(back, "ATAC")),
               dim(SingleCellExperiment::altExp(sce, "ATAC")))
})

test_that("a modality R cannot build a feature table for is an error, not a bad file", {
  skip_if_not_installed("SingleCellExperiment")
  m <- Matrix::rsparsematrix(4, 3, density = 0.6)
  dimnames(m) <- list(paste0("f", 1:4), paste0("c", 1:3))
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = m))
  expect_error(write_cytome(sce, tempfile(), main_modality = "PROTEIN"),
               "cannot write modality")
})

test_that("ATAC feature ids that are not coordinates fail with the offending id", {
  skip_if_not_installed("SingleCellExperiment")
  m <- Matrix::rsparsematrix(3, 3, density = 0.6)
  dimnames(m) <- list(c("chr1:1-2", "notapeak", "chr2:5-9"), paste0("c", 1:3))
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = m))
  expect_error(write_cytome(sce, tempfile(), main_modality = "ATAC"), "notapeak")
})
