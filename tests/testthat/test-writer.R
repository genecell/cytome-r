# The writer (audit part 4, major 23): values by name, never by position; a matrix or graph
# over only some of the features or cells left out, with a warning (decision 4); one transaction in a
# temporary file, renamed into place; column types kept; graphs as Python stores them; the
# writable modalities written and the others named.

.counts <- function(ng = 6, nc = 30, seed = 1) {
  set.seed(seed)
  m <- Matrix::Matrix(stats::rpois(ng * nc, 2), ng, nc, sparse = TRUE)
  dimnames(m) <- list(c("Sox2", "Pax6", "Gfap", "Aqp4", "Olig2", "Snap25")[seq_len(ng)],
                      paste0("cell", seq_len(nc)))
  m
}
.con <- function(p) DBI::dbConnect(RSQLite::SQLite(), p, flags = RSQLite::SQLITE_RO)

# Stored over all the genes, the zeros of the genes a scale.data lacks would read as values
# (0 is the mean there), and no reader checked a record of the genes it covered.
test_that("a scale.data over some of the genes is left out, with a warning that names it", {
  skip_if_not_installed("SeuratObject")
  cnt <- .counts()
  so <- SeuratObject::CreateSeuratObject(counts = cnt)
  so <- SeuratObject::SetAssayData(so, layer = "data", new.data = log1p(cnt))
  sd <- as.matrix(log1p(cnt))[c("Olig2", "Gfap"), ] - 0.5               # what ScaleData keeps
  so <- SeuratObject::SetAssayData(so, layer = "scale.data", new.data = sd)
  p <- tempfile(fileext = ".cytome")
  w <- capture_warnings(write_cytome(so, p, layers = TRUE))
  expect_identical(w, paste("cytome: RNA layer 'scale.data' is not written: it covers only 2 of 6",
                            "features. Recompute it after reading, with ScaleData()."))
  x <- cytome_open(p); on.exit(cytome_close(x))
  expect_setequal(cytome_matrices(x)$matrix_name, c("RNA_counts", "RNA_data"))
  keys <- DBI::dbGetQuery(x$con, "SELECT key FROM _metadata")$key
  expect_false(any(startsWith(keys, "layer_features")))
  expect_equal(cytome_var(x, "RNA_counts")$gene_id, rownames(cnt))
  expect_equal(unname(as.matrix(read_cytome_matrix(x, "RNA_counts"))), unname(as.matrix(cnt)))
  expect_equal(unname(as.matrix(read_cytome_matrix(x, "RNA_data"))), unname(as.matrix(log1p(cnt))),
               tolerance = 1e-6)
})

test_that("a partial layer of an SCTransform assay is recomputed with SCTransform()", {
  cnt <- .counts()
  sd <- as.matrix(log1p(cnt))[c("Olig2", "Gfap"), ]
  w <- capture_warnings(cytome:::.spec_modality("RNA", rownames(cnt), NULL, cnt, list(scale.data = sd),
                                                colnames(cnt), "assay 'RNA'", sct = TRUE))
  expect_identical(w, paste("cytome: RNA layer 'scale.data' is not written: it covers only 2 of 6",
                            "features. Recompute it after reading, with SCTransform(); the SCT model",
                            "is not stored."))
})

# Seurat puts a layer it is given into the assay's gene order, so the writer is handed one
# in another order directly: its values must land under their own genes (Gfap's values
# were once stored under Sox2).
test_that("a layer over all the genes, in another order, is stored under its own genes", {
  cnt <- .counts()
  lg <- log1p(cnt)[rev(rownames(cnt)), ]                                 # Snap25 first, Sox2 last
  expect_no_warning(spec <- cytome:::.spec_modality("RNA", rownames(cnt), NULL, cnt, list(data = lg),
                                                    colnames(cnt), "assay 'RNA'"))
  p <- tempfile(fileext = ".cytome")
  cytome:::.write_store(p, list(barcodes = colnames(cnt), obs = NULL, modalities = list(spec),
                                embeddings = list(), graphs = list()), "float32", 10L, "zstd")
  x <- cytome_open(p); on.exit(cytome_close(x))
  M <- read_cytome_matrix(x, "RNA_data")                                 # genes x cells
  rownames(M) <- cytome_var(x, "RNA_data")$gene_id
  expect_equal(rownames(M), rownames(cnt))
  expect_equal(as.numeric(M["Gfap", ]), as.numeric(log1p(cnt)["Gfap", ]), tolerance = 1e-6)
  expect_equal(unname(as.matrix(M)), unname(as.matrix(log1p(cnt))), tolerance = 1e-6)
})

# The same rule for cells: an assay, layer or graph over some of the cells was stored with 0
# for the others (a cell without counts, or without neighbours), and before that by position.
test_that("an assay, layer or graph over only some of the cells is left out, with a warning", {
  skip_if_not_installed("SeuratObject")
  cnt <- .counts()
  n <- ncol(cnt); half <- colnames(cnt)[16:n]
  so <- SeuratObject::CreateSeuratObject(counts = cnt)
  atac <- Matrix::Matrix(rep(1:3, length.out = 3 * n), 3, n, sparse = TRUE,
                         dimnames = list(c("chr1:1-5", "chr1:9-20", "chr2:3-9"), colnames(cnt)))
  so[["ATAC"]] <- SeuratObject::CreateAssay5Object(counts = atac[, half])   # one batch only
  adt <- Matrix::Matrix(5, 2, n, sparse = TRUE, dimnames = list(c("CD3", "CD19"), colnames(cnt)))
  so[["ADT"]] <- SeuratObject::CreateAssay5Object(counts = adt[, half])
  SeuratObject::LayerData(so, assay = "RNA", layer = "data") <- log1p(cnt)[, half]
  G <- SeuratObject::as.Graph(Matrix::sparseMatrix(i = 1:2, j = 2:1, x = 1, dims = c(15, 15),
                                                   dimnames = list(half, half)))
  SeuratObject::DefaultAssay(G) <- "ATAC"
  so[["ATAC_snn"]] <- G
  p <- tempfile(fileext = ".cytome")
  w <- capture_warnings(write_cytome(so, p, layers = TRUE))
  some <- function(what) paste("cytome:", what, "is not written: it covers only 15 of 30 cells.",
                               "To keep it, subset the object to those cells and write that as a",
                               "separate store.")
  expect_true(some("assay 'ATAC'") %in% w)
  expect_true(some("RNA layer 'data'") %in% w)
  expect_true(some("graph 'ATAC_snn'") %in% w)
  expect_true(any(grepl("no feature table for them: ADT", w)))           # not writable, as before
  expect_length(w, 4L)
  x <- cytome_open(p); on.exit(cytome_close(x))
  expect_equal(cytome_matrices(x)$matrix_name, "RNA_counts")
  expect_length(cytome_graphs(x), 0L)
  expect_equal(cytome_obs(x)$barcode, colnames(cnt))
  expect_equal(unname(as.matrix(read_cytome_matrix(x, "RNA_counts"))), unname(as.matrix(cnt)))
  ## the store is written around the main assay, which must cover every cell
  expect_error(write_cytome(so, tempfile(fileext = ".cytome"), assay = "ATAC"),
               "covers only 15 of the object's 30 cells")
})

# Seurat puts an assay's cells into the object's order, so the writer is handed counts in
# another order directly.
test_that("counts over all the cells, in another order, are stored under their own cells", {
  cnt <- .counts()
  expect_no_warning(spec <- cytome:::.spec_modality("RNA", rownames(cnt), NULL, cnt[, rev(colnames(cnt))],
                                                    list(), colnames(cnt), "assay 'RNA'", main = TRUE))
  p <- tempfile(fileext = ".cytome")
  cytome:::.write_store(p, list(barcodes = colnames(cnt), obs = NULL, modalities = list(spec),
                                embeddings = list(), graphs = list()), "float32", 7L, "zstd")
  x <- cytome_open(p); on.exit(cytome_close(x))
  expect_equal(cytome_obs(x)$barcode, colnames(cnt))
  expect_equal(unname(as.matrix(read_cytome_matrix(x, "RNA_counts"))), unname(as.matrix(cnt)))
})

# assay(withDimnames = FALSE) returns the names an assay had when it was set, and renaming
# the object's genes or cells (Ensembl ids to symbols, as most SCE workflows do) leaves them
# behind: the write failed, or put each row's values under the name the row had before.
test_that("an SCE whose genes and cells were renamed is written under its current names", {
  skip_if_not_installed("SingleCellExperiment")
  cnt <- .counts()
  atac <- Matrix::Matrix(rep(1:4, length.out = 2 * ncol(cnt)), 2, ncol(cnt), sparse = TRUE,
                         dimnames = list(c("chr1:1-5", "chr2:3-9"), colnames(cnt)))
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = cnt, logcounts = log1p(cnt)))
  SingleCellExperiment::altExp(sce, "ATAC") <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = atac))
  rownames(sce) <- rev(rownames(cnt))                                    # relabelled, not reordered
  colnames(sce) <- paste0("bc", seq_len(ncol(sce)))
  rownames(SingleCellExperiment::altExp(sce, "ATAC")) <- c("chr3:1-5", "chr4:3-9")
  p <- tempfile(fileext = ".cytome")
  write_cytome(sce, p, layers = TRUE)
  x <- cytome_open(p); on.exit(cytome_close(x))
  expect_equal(cytome_obs(x)$barcode, colnames(sce))
  expect_equal(cytome_var(x, "RNA_counts")$gene_id, rownames(sce))
  for (a in c("counts", "logcounts"))
    expect_equal(unname(as.matrix(read_cytome_matrix(x, paste0("RNA_", a)))),
                 unname(as.matrix(SummarizedExperiment::assay(sce, a))), tolerance = 1e-6)
  expect_equal(cytome_var(x, "ATAC_counts")$peak_id, c("chr3:1-5", "chr4:3-9"))
  expect_equal(unname(as.matrix(read_cytome_matrix(x, "ATAC_counts"))), unname(as.matrix(atac)))
})

test_that("a write that fails part-way leaves the old file as it was, and nothing else", {
  skip_if_not_installed("SingleCellExperiment")
  cnt <- .counts()
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = cnt))
  atac <- Matrix::Matrix(1, 2, ncol(cnt), sparse = TRUE, dimnames = list(c("chr1:1-5", "notapeak"), colnames(cnt)))
  SingleCellExperiment::altExp(sce, "ATAC") <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = atac))
  d <- tempfile("write_"); dir.create(d)
  p <- file.path(d, "old.cytome")
  writeLines("the old file", p)
  expect_error(write_cytome(sce, p), "notapeak")                          # after the RNA was written
  expect_equal(readLines(p), "the old file")
  expect_equal(list.files(d, all.files = TRUE, no.. = TRUE), "old.cytome")
})

test_that("column types and names survive, and sample_id fills the sample column", {
  skip_if_not_installed("SingleCellExperiment")
  cnt <- .counts()
  n <- ncol(cnt)
  cd <- S4Vectors::DataFrame(
    cluster = factor(rep(c("b", "a"), length.out = n), levels = c("b", "a", "unused")),
    stage = factor(rep(c("late", "early"), length.out = n), levels = c("early", "late"), ordered = TRUE),
    is_doublet = rep(c(TRUE, FALSE, NA), length.out = n),
    n_umi = c(NA, seq_len(n - 1L)),
    score = c(NaN, NA, seq_len(n - 2L) / 7),
    label = c(NA, rep("xä", n - 1L)),
    day = as.Date("2026-10-27") + seq_len(n),
    "a.b" = seq_len(n), "a_b" = -seq_len(n), "col/with slash" = seq_len(n) * 2,
    sample_id = rep(c("donorA", "donorB"), each = n / 2), check.names = FALSE)
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = cnt), colData = cd)
  p <- tempfile(fileext = ".cytome")
  write_cytome(sce, p)
  x <- cytome_open(p); on.exit(cytome_close(x))
  obs <- cytome_obs(x)
  expect_equal(levels(obs$cluster), c("b", "a", "unused"))
  expect_equal(as.character(obs$cluster), as.character(cd$cluster))
  expect_true(is.ordered(obs$stage)); expect_equal(levels(obs$stage), c("early", "late"))
  expect_identical(obs$is_doublet, cd$is_doublet)
  expect_identical(obs$n_umi, cd$n_umi)
  expect_equal(obs$score[-(1:2)], cd$score[-(1:2)]); expect_true(all(is.na(obs$score[1:2])))
  expect_identical(obs$label, cd$label)
  expect_identical(obs$day, as.character(cd$day))
  expect_identical(obs[["a.b"]], cd[["a.b"]]); expect_identical(obs[["a_b"]], cd[["a_b"]])
  expect_identical(obs[["col/with slash"]], cd[["col/with slash"]])
  expect_identical(obs$sample_id, cd$sample_id)
  s <- DBI::dbGetQuery(x$con, "SELECT sample_id, n_cells FROM samples ORDER BY sample_idx")
  expect_equal(s$sample_id, c("donorA", "donorB")); expect_equal(s$n_cells, c(n / 2, n / 2))
  ## SQLite column names ignore case: the second of Batch and batch is renamed, with a warning
  cd2 <- S4Vectors::DataFrame(Batch = rep("B1", n), batch = rep("b1", n))
  sce2 <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = cnt), colData = cd2)
  p2 <- tempfile(fileext = ".cytome")
  expect_warning(write_cytome(sce2, p2), "'batch' is stored as 'batch_'")
})

test_that("modalities cytome has no table for are skipped by default and refused by name", {
  skip_if_not_installed("SingleCellExperiment")
  cnt <- .counts()
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = cnt))
  ercc <- Matrix::Matrix(1, 2, ncol(cnt), sparse = TRUE, dimnames = list(c("ERCC-1", "ERCC-2"), colnames(cnt)))
  SingleCellExperiment::altExp(sce, "ERCC") <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = ercc))
  atac <- Matrix::Matrix(2, 2, ncol(cnt), sparse = TRUE, dimnames = list(c("chr1:1-5", "chr2:3-9"), colnames(cnt)))
  SingleCellExperiment::altExp(sce, "ATAC") <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = atac))
  p <- tempfile(fileext = ".cytome")
  expect_warning(write_cytome(sce, p), "ERCC")
  x <- cytome_open(p); on.exit(cytome_close(x))
  expect_setequal(cytome_matrices(x)$matrix_name, c("RNA_counts", "ATAC_counts"))
  expect_error(write_cytome(sce, tempfile(), alt_modalities = "ERCC"), "cannot write modality 'ERCC'")
  expect_error(write_cytome(sce, tempfile(), alt_modalities = "ADT"), "no such altExp: ADT")
})

test_that("Seurat: an SCT assay is skipped, split layers are joined, the default assay must be writable", {
  skip_if_not_installed("SeuratObject")
  cnt <- .counts()
  so <- SeuratObject::CreateSeuratObject(counts = cnt)
  so[["SCT"]] <- SeuratObject::CreateAssay5Object(counts = cnt)
  p <- tempfile(fileext = ".cytome")
  expect_warning(write_cytome(so, p), "SCT")
  SeuratObject::DefaultAssay(so) <- "SCT"
  expect_error(write_cytome(so, tempfile()), "assay = \"RNA\"")
  sp <- SeuratObject::CreateSeuratObject(counts = cnt)
  sp$grp <- rep(c("A", "B"), length.out = ncol(cnt))
  sp[["RNA"]] <- split(sp[["RNA"]], f = sp$grp)
  p2 <- tempfile(fileext = ".cytome")
  write_cytome(sp, p2)
  x <- cytome_open(p2); on.exit(cytome_close(x))
  M <- read_cytome_matrix(x, "RNA_counts")
  expect_equal(unname(as.matrix(M)), unname(as.matrix(cnt)))
  expect_equal(cytome_obs(x)$barcode, colnames(cnt))
})

test_that("graphs are stored as Python stores them, and Python's to_anndata is told", {
  skip_if_not_installed("SingleCellExperiment")
  cnt <- .counts()
  n <- ncol(cnt)
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = cnt))
  G <- Matrix::sparseMatrix(i = c(1L, 2L, 5L, 30L), j = c(2L, 1L, 9L, 3L), x = c(0.5, 0.5, 1, 0.25), dims = c(n, n))
  SingleCellExperiment::colPair(sce, "knn") <- methods::as(G, "dgCMatrix")
  p <- tempfile(fileext = ".cytome")
  write_cytome(sce, p)
  con <- .con(p); on.exit(DBI::dbDisconnect(con))
  gm <- DBI::dbGetQuery(con, "SELECT graph_name, axis, n_rows, n_cols, n_nonzero, dtype FROM graph_meta")
  expect_equal(gm$graph_name, "knn"); expect_equal(gm$n_rows, n); expect_equal(gm$n_nonzero, 4)
  expect_equal(gm$dtype, "float64")
  expect_equal(DBI::dbGetQuery(con, "SELECT COUNT(*) n FROM graph_edges")$n, 0)
  om <- DBI::dbGetQuery(con, "SELECT value FROM _metadata WHERE key = '_anndata_obsp_map'")$value
  expect_equal(jsonlite::fromJSON(om), list(knn = "knn"))
  x <- cytome_open(p); on.exit(cytome_close(x), add = TRUE)
  expect_equal(as.matrix(cytome_graph(x, "knn")), as.matrix(G))
})

test_that("dtype is checked, int32 holds whole numbers, and is_integer is recorded", {
  skip_if_not_installed("SingleCellExperiment")
  cnt <- .counts()
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = cnt))
  expect_error(write_cytome(sce, tempfile(), dtype = "int64"), "dtype must be")
  p <- tempfile(fileext = ".cytome")
  write_cytome(sce, p, dtype = "int32")
  x <- cytome_open(p); on.exit(cytome_close(x))
  m <- DBI::dbGetQuery(x$con, "SELECT dtype, is_integer FROM matrix_meta")
  expect_equal(m$dtype, "int32"); expect_equal(m$is_integer, 1L)
  expect_equal(unname(as.matrix(read_cytome_matrix(x, "RNA_counts"))), unname(as.matrix(cnt)))
  frac <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = cnt / 3))
  expect_error(write_cytome(frac, tempfile(), dtype = "int32"), "not whole numbers")
})

test_that("the written store is a single file, and an old store's -wal is not replayed over it", {
  skip_if_not_installed("SingleCellExperiment")
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = .counts()))
  d <- tempfile("store_"); dir.create(d)
  p <- file.path(d, "s.cytome")
  for (f in paste0(p, c("-wal", "-shm"))) writeLines("stale", f)
  write_cytome(sce, p)
  expect_equal(list.files(d, all.files = TRUE, no.. = TRUE), "s.cytome")
  expect_equal(as.integer(readBin(p, "raw", 20L)[19:20]), c(1L, 1L))     # not WAL
})

# The writer filled each obs column with one UPDATE per row and inserted each graph edge,
# each row in its own transaction: 13.4 s for 20,000 cells with five columns and a 15-NN
# graph, against 0.3 s without them. It now writes the store in one transaction, and they
# add a fraction of a second. Timed as what they add, with a generous bound, so that only
# the row-by-row form fails it.
test_that("the store is written in one transaction", {
  skip_on_cran()
  skip_if_not_installed("SingleCellExperiment")
  n <- 20000L
  set.seed(3)
  m <- Matrix::rsparsematrix(20, n, density = 0.05, rand.x = function(k) stats::rpois(k, 2) + 1)
  dimnames(m) <- list(paste0("g", 1:20), paste0("c", seq_len(n)))
  cd <- S4Vectors::DataFrame(a = stats::runif(n), b = sample(letters, n, TRUE), c = seq_len(n),
                             d = stats::runif(n) > 0.5, e = factor(sample(c("x", "y"), n, TRUE)))
  bare <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = m))
  full <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = m), colData = cd)
  G <- Matrix::sparseMatrix(i = rep(seq_len(n), each = 15L), j = sample.int(n, 15L * n, TRUE),
                            x = stats::runif(15L * n), dims = c(n, n))
  SingleCellExperiment::colPair(full, "knn") <- methods::as(G, "dgCMatrix")
  write_cytome(bare, tempfile(fileext = ".cytome"))                     # warm up
  t0 <- system.time(write_cytome(bare, tempfile(fileext = ".cytome")))[["elapsed"]]
  t1 <- system.time(write_cytome(full, tempfile(fileext = ".cytome")))[["elapsed"]]
  expect_lt(t1 - t0, 4)
})

test_that("an addNA() level is missing, not a category, and an unnamed column gets a name", {
  skip_if_not_installed("SingleCellExperiment")
  m <- Matrix::Matrix(1, 2, 3, sparse = TRUE, dimnames = list(c("g1", "g2"), c("a", "b", "c")))
  s <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = m), colData = S4Vectors::DataFrame(f = addNA(factor(c("x", NA, "y")))))
  p <- tempfile(fileext = ".cytome")
  write_cytome(s, p)
  x <- cytome_open(p); on.exit(cytome_close(x))
  o <- cytome_obs(x)
  expect_equal(levels(o$f), c("x", "y"))
  expect_equal(as.character(o$f), c("x", NA, "y"))
  df <- data.frame(a = 1:3); names(df) <- ""
  p2 <- tempfile(fileext = ".cytome")
  cytome:::.write_store(p2, list(barcodes = c("a", "b", "c"), obs = df,
                                 modalities = list(cytome:::.spec_modality("RNA", c("g1", "g2"), NULL, m, list(),
                                                                           c("a", "b", "c"), "counts")),
                                 embeddings = list(), graphs = list()), "float32", 10L, "zstd")
  y <- cytome_open(p2); on.exit(cytome_close(y), add = TRUE)
  expect_equal(cytome_obs(y)$V1, 1:3)
})
