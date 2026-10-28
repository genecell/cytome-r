# Both directions of the format, with the Python package itself: the fixtures test what R
# reads from Python; this writes a store from R and has Python read it back. It needs a
# Python with cytome and anndata (set CYTOME_PYTHON to its interpreter), so it is skipped
# where there is none -- on CRAN and r-universe, for instance.

.python_with_cytome <- function() {
  py <- Sys.getenv("CYTOME_PYTHON", unset = "")
  if (!nzchar(py)) py <- unname(Sys.which("python3"))
  if (!nzchar(py)) return(NULL)
  ok <- tryCatch(suppressWarnings(system2(py, c("-c", shQuote("import cytome, anndata")),
                                          stdout = FALSE, stderr = FALSE)) == 0L,
                 error = function(e) FALSE)
  if (isTRUE(ok)) py else NULL
}

test_that("a store written by R reads in Python as R wrote it", {
  skip_if_not_installed("SingleCellExperiment")
  py <- .python_with_cytome()
  skip_if(is.null(py), "no Python with cytome (set CYTOME_PYTHON)")
  set.seed(5)
  n <- 40L
  cnt <- Matrix::rsparsematrix(7, n, 0.4, rand.x = function(k) stats::rpois(k, 3) + 1)
  dimnames(cnt) <- list(paste0("G", 1:7), sprintf("AAAC-%02d", 1:n))
  atac <- Matrix::rsparsematrix(3, n, 0.5, rand.x = function(k) stats::rpois(k, 1) + 1)
  dimnames(atac) <- list(c("chr1:100-200", "chr1:300-450", "chr2:10-90"), colnames(cnt))
  cd <- S4Vectors::DataFrame(
    cluster = factor(sample(c("b", "a"), n, TRUE), levels = c("b", "a", "unused")),
    stage = factor(sample(c("early", "late"), n, TRUE), levels = c("early", "late"), ordered = TRUE),
    is_doublet = sample(c(TRUE, FALSE), n, TRUE),
    n_umi = as.integer(Matrix::colSums(cnt)),
    score = c(NA, stats::runif(n - 1L)),
    label = c(NA, rep("xä", n - 1L)))
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = cnt, logcounts = log1p(cnt)),
                                                    colData = cd)
  SingleCellExperiment::altExp(sce, "ATAC") <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = atac))
  pca <- matrix(stats::rnorm(n * 3), n, 3)
  SingleCellExperiment::reducedDim(sce, "PCA") <- pca
  G <- Matrix::sparseMatrix(i = c(1L, 2L, 40L), j = c(2L, 1L, 7L), x = c(0.5, 0.5, 1), dims = c(n, n))
  SingleCellExperiment::colPair(sce, "knn") <- methods::as(G, "dgCMatrix")
  p <- tempfile(fileext = ".cytome")
  write_cytome(sce, p, layers = TRUE)

  rowmajor <- function(M) as.vector(t(as.matrix(M)))                     # cells x features, by row
  obs <- function(v, kind, ...) list(kind = kind, values = as.list(v), ...)
  e <- list(
    cells = colnames(cnt),
    matrices = list(RNA = list(shape = c(n, 7L), values = rowmajor(Matrix::t(cnt)),
                               features = rownames(cnt), id_col = "gene_id"),
                    ATAC = list(shape = c(n, 3L), values = rowmajor(Matrix::t(atac)),
                                features = rownames(atac), id_col = "peak_id")),
    layers = list(logcounts = rowmajor(Matrix::t(log1p(cnt)))),
    obs = list(cluster = obs(as.character(cd$cluster), "factor", levels = levels(cd$cluster), ordered = FALSE),
               stage = obs(as.character(cd$stage), "factor", levels = levels(cd$stage), ordered = TRUE),
               is_doublet = obs(cd$is_doublet, "logical"),
               n_umi = obs(cd$n_umi, "integer"),
               score = obs(cd$score, "double"),
               label = obs(cd$label, "character")),
    obsm = list(X_pca = as.vector(t(pca))),
    obsp = list(knn = as.vector(t(as.matrix(G)))))
  ej <- tempfile(fileext = ".json")
  writeLines(jsonlite::toJSON(e, auto_unbox = TRUE, digits = NA, na = "null"), ej)
  out <- suppressWarnings(system2(py, c(test_path("python", "check_r_store.py"), p, ej),
                                  stdout = TRUE, stderr = TRUE))
  expect_null(attr(out, "status"), info = paste(out, collapse = "\n"))
  expect_equal(utils::tail(out, 1), "ok", info = paste(out, collapse = "\n"))
})
