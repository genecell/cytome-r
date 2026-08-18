# Layers and graphs. Both were "no" in the transfer table; neither was a
# format limitation. Layers are opt-in because an archive and a handoff want
# opposite things; graphs are on by default because recomputing a KNN graph on
# millions of cells is the expensive part.

## colPair<- wants a SelfHits with integer from/to, not a sparse matrix.
.hits <- function(G, n) {
  tg <- methods::as(G, "TsparseMatrix")
  h <- S4Vectors::SelfHits(from = as.integer(tg@i) + 1L,
                           to = as.integer(tg@j) + 1L, nnode = n)
  S4Vectors::mcols(h)$value <- as.numeric(tg@x)
  h
}

.mini_sce <- function(n_g = 8, n_c = 12) {
  set.seed(2)
  m <- abs(Matrix::rsparsematrix(n_g, n_c, density = 0.5))
  dimnames(m) <- list(paste0("g", seq_len(n_g)), paste0("c", seq_len(n_c)))
  SingleCellExperiment::SingleCellExperiment(assays = list(counts = m))
}

test_that("layers are off by default", {
  skip_if_not_installed("SingleCellExperiment")
  sce <- .mini_sce()
  SummarizedExperiment::assay(sce, "logcounts") <- log1p(
    SummarizedExperiment::assay(sce, "counts"))
  p <- tempfile(fileext = ".cytome")
  write_cytome(sce, p)
  x <- cytome_open(p); on.exit(cytome_close(x))
  expect_equal(cytome_matrices(x)$matrix_name, "RNA_counts")
})

test_that("layers = TRUE carries them, and they come back as assays", {
  skip_if_not_installed("SingleCellExperiment")
  sce <- .mini_sce()
  SummarizedExperiment::assay(sce, "logcounts") <- log1p(
    SummarizedExperiment::assay(sce, "counts"))
  p <- tempfile(fileext = ".cytome")
  write_cytome(sce, p, layers = TRUE)

  x <- cytome_open(p)
  expect_true("RNA_logcounts" %in% cytome_matrices(x)$matrix_name)
  cytome_close(x)

  back <- read_cytome(p)
  expect_true("logcounts" %in% SummarizedExperiment::assayNames(back))
  expect_equal(as.matrix(SummarizedExperiment::assay(back, "logcounts")),
               as.matrix(SummarizedExperiment::assay(sce, "logcounts")),
               tolerance = 1e-6, check.attributes = FALSE)
})

test_that("naming a layer that does not exist says which, and what does", {
  skip_if_not_installed("SingleCellExperiment")
  expect_error(write_cytome(.mini_sce(), tempfile(), layers = "lognorm"),
               "no such layer")
})

test_that("an SCE colPair round-trips through graph_edges", {
  skip_if_not_installed("SingleCellExperiment")
  sce <- .mini_sce()
  n <- ncol(sce)
  set.seed(3)
  G <- Matrix::sparseMatrix(i = c(1L, 2L, 3L), j = c(4L, 5L, 6L),
                            x = c(0.5, 0.25, 1), dims = c(n, n))
  SingleCellExperiment::colPair(sce, "knn") <- .hits(G, n)

  p <- tempfile(fileext = ".cytome")
  write_cytome(sce, p)
  x <- cytome_open(p)
  expect_true("knn" %in% cytome_graphs(x))
  expect_equal(sum(cytome_graph(x, "knn")), sum(G))
  cytome_close(x)

  back <- read_cytome(p)
  expect_true("knn" %in% SingleCellExperiment::colPairNames(back))
})

test_that("graphs = FALSE opts out", {
  skip_if_not_installed("SingleCellExperiment")
  sce <- .mini_sce(); n <- ncol(sce)
  SingleCellExperiment::colPair(sce, "knn") <-
    .hits(Matrix::sparseMatrix(i = 1L, j = 2L, x = 1, dims = c(n, n)), n)
  p <- tempfile(fileext = ".cytome")
  write_cytome(sce, p, graphs = FALSE)
  x <- cytome_open(p); on.exit(cytome_close(x))
  expect_equal(cytome_graphs(x), character(0))
})

test_that("a Seurat graph round-trips", {
  skip_if_not_installed("SeuratObject")
  ref <- system.file("extdata", "reference.cytome", package = "cytome")
  skip_if(ref == "")
  so <- read_cytome(ref, as = "Seurat")
  n <- ncol(so)
  G <- Matrix::sparseMatrix(i = c(1L, 2L), j = c(3L, 4L), x = c(1, 0.5),
                            dims = c(n, n), dimnames = list(colnames(so), colnames(so)))
  so@graphs[["RNA_nn"]] <- SeuratObject::as.Graph(G)

  p <- tempfile(fileext = ".cytome")
  write_cytome(so, p)
  x <- cytome_open(p)
  expect_true("RNA_nn" %in% cytome_graphs(x))
  cytome_close(x)

  back <- read_cytome(p, as = "Seurat")
  expect_true("RNA_nn" %in% names(back@graphs))
})
