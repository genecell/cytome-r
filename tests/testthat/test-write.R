test_that("write_cytome_sce round-trips (R write -> R read) for all codecs", {
  skip_if_not_installed("SingleCellExperiment")
  suppressMessages(library(SingleCellExperiment)); library(Matrix)
  set.seed(1); nc <- 40; ng <- 25
  counts <- as(Matrix::rsparsematrix(ng, nc, 0.2, rand.x = function(n) rpois(n, 3) + 1L), "CsparseMatrix")
  counts@x <- abs(counts@x); rownames(counts) <- paste0("G", 1:ng); colnames(counts) <- paste0("C", 1:nc)
  sce <- SingleCellExperiment(assays = list(counts = counts),
                              colData = DataFrame(grp = sample(c("a", "b"), nc, TRUE)),
                              reducedDims = list(PCA = matrix(rnorm(nc * 4), nc, 4)))
  ref <- as.matrix(Matrix::t(counts))
  for (cp in c("zstd", "lz4", "zlib")) {
    p <- tempfile(fileext = ".cytome")
    write_cytome(sce, p, compression = cp, chunk_size = 15L)
    ds <- cytome_open(p)
    rt <- as.matrix(Matrix::t(read_cytome_matrix(ds, "RNA_counts")))
    expect_equal(unname(rt), unname(ref), info = cp)
    expect_equal(nrow(cytome_obs(ds)), nc)
    cytome_close(ds); unlink(p)
  }
})

test_that("write_cytome_seurat round-trips", {
  skip_if_not_installed("SeuratObject"); library(Matrix)
  set.seed(2); counts <- as(Matrix::rsparsematrix(20, 30, 0.25, rand.x = function(n) rpois(n, 2) + 1L), "CsparseMatrix")
  counts@x <- abs(counts@x); rownames(counts) <- paste0("G", 1:20); colnames(counts) <- paste0("C", 1:30)
  obj <- SeuratObject::CreateSeuratObject(counts)
  p <- tempfile(fileext = ".cytome"); write_cytome(obj, p)
  ds <- cytome_open(p)
  rt <- as.matrix(Matrix::t(read_cytome_matrix(ds, "RNA_counts")))
  expect_equal(unname(rt), unname(as.matrix(Matrix::t(counts))))
  cytome_close(ds); unlink(p)
})

test_that("write_cytome_sce populates symbol from rowData symbol column", {
  skip_if_not_installed("SingleCellExperiment")
  suppressMessages(library(SingleCellExperiment)); library(Matrix)
  set.seed(7); nc <- 12; ng <- 6
  counts <- as(Matrix::rsparsematrix(ng, nc, 0.4, rand.x = function(n) rpois(n, 2) + 1L), "CsparseMatrix")
  counts@x <- abs(counts@x)
  rownames(counts) <- paste0("ENSG", sprintf("%05d", 1:ng))   # ids = Ensembl
  colnames(counts) <- paste0("C", 1:nc)
  sce <- SingleCellExperiment(assays = list(counts = counts))
  SummarizedExperiment::rowData(sce)$Symbol <- paste0("Gene", 1:ng)  # symbols under 'Symbol'
  p <- tempfile(fileext = ".cytome")
  write_cytome(sce, p)
  con <- DBI::dbConnect(RSQLite::SQLite(), p)
  g <- DBI::dbGetQuery(con, "SELECT gene_id, symbol FROM genes ORDER BY gene_idx")
  DBI::dbDisconnect(con); unlink(p)
  expect_equal(g$gene_id, rownames(counts))            # id stays Ensembl
  expect_equal(g$symbol, paste0("Gene", 1:ng))         # symbol from rowData, not a copy of id
})
