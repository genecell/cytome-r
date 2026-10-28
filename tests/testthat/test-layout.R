# Column-major matrices: CSC chunks of features in matrix_csc_chunks, no row chunks. This is
# the default layout of PIASO's importer and quantifier for ATAC peaks and tiles. The R
# package read only matrix_chunks, so every such matrix came back as an empty dgCMatrix with
# no error, streamed as zero chunks, and summed to zero through cytome_delayed() (audit
# part 4, critical 3). The fixture is written by PIASO's Rust importer and quantifier.

test_that("column-major matrices written by PIASO read as Python reads them", {
  p <- fixture_store("piaso_atac"); e <- fixture_expected("piaso_atac")
  x <- cytome_open(p); on.exit(cytome_close(x))
  for (nm in c("ATAC_counts", "tiles_counts")) {
    expect_equal(e$layouts[[nm]], "csc")
    want <- expected_matrix(e$matrices[[nm]])                     # cells x features
    expect_gt(length(want@x), 0)
    M <- read_cytome_matrix(x, nm, features_x_cells = FALSE)
    expect_equal(dim(M), dim(want), info = nm)
    expect_equal(unname(as.matrix(M)), as.matrix(want), info = nm)
    expect_equal(unname(as.matrix(read_cytome_matrix(x, nm))), t(as.matrix(want)), info = nm)
  }
})

test_that("a column-major matrix streams along features, and says so along cells", {
  p <- fixture_store("piaso_atac"); e <- fixture_expected("piaso_atac")
  x <- cytome_open(p); on.exit(cytome_close(x))
  want <- expected_matrix(e$matrices$ATAC_counts)
  expect_error(cytome_stream(x, "ATAC_counts", function(m, s, e, k) sum(m)), "along = \"features\"")
  sums <- cytome_stream(x, "ATAC_counts", function(m, s, e, k) Matrix::colSums(m), along = "features")
  expect_gt(length(sums), 1L)                                     # really several chunks
  expect_equal(unname(unlist(sums)), unname(Matrix::colSums(want)))
  spans <- cytome_stream(x, "ATAC_counts", function(m, s, e, k) c(s, e, ncol(m)), along = "features")
  expect_equal(spans[[1]][1], 0L)
  expect_equal(spans[[length(spans)]][2], ncol(want))
  expect_true(all(vapply(spans, function(v) v[2] - v[1] == v[3], logical(1))))
  tm <- cytome_stream(x, "ATAC_counts", function(m, s, e, k) dim(m), along = "features",
                      features_x_cells = TRUE)
  expect_equal(tm[[1]][2], nrow(want))                            # features x cells
})

test_that("a row-major matrix with a feature index streams along both axes", {
  p <- fixture_store("py_types"); e <- fixture_expected("py_types")
  x <- cytome_open(p); on.exit(cytome_close(x))
  want <- expected_matrix(e$matrices$RNA_log1p)
  by_cells <- cytome_stream(x, "RNA_log1p", function(m, s, e, k) Matrix::rowSums(m))
  by_feats <- cytome_stream(x, "RNA_log1p", function(m, s, e, k) Matrix::colSums(m), along = "features")
  expect_equal(unname(unlist(by_cells)), unname(Matrix::rowSums(want)), tolerance = 1e-6)
  expect_equal(unname(unlist(by_feats)), unname(Matrix::colSums(want)), tolerance = 1e-6)
  expect_error(cytome_stream(x, "RNA_counts", identity, along = "features"), "no column chunks")
})

test_that("a column-major matrix backs a DelayedArray, chunked along features", {
  skip_if_not_installed("DelayedArray")
  skip_if_not_installed("SparseArray")
  p <- fixture_store("piaso_atac"); e <- fixture_expected("piaso_atac")
  want <- t(as.matrix(expected_matrix(e$matrices$tiles_counts)))   # features x cells
  d <- cytome_delayed(p, "tiles_counts")
  expect_equal(unname(as.matrix(d)), unname(want))
  cd <- DelayedArray::chunkdim(d)
  expect_equal(cd[2], ncol(want))                                 # every cell in a chunk
  expect_equal(unname(as.numeric(DelayedArray::colSums(d))), unname(colSums(want)))
  i <- c(9L, 2L, 2L, 14L); j <- c(30L, 1L, 7L, 7L)
  expect_equal(unname(as.matrix(d[i, j])), unname(want[i, j]))
})

test_that("read_cytome builds an object from a PIASO ATAC store", {
  skip_if_not_installed("SingleCellExperiment")
  p <- fixture_store("piaso_atac"); e <- fixture_expected("piaso_atac")
  expect_error(read_cytome(p), "ATAC")                            # no RNA; the message names what there is
  sce <- read_cytome(p, main_modality = "ATAC")
  expect_equal(rownames(sce), e$peaks)
  expect_equal(unname(as.matrix(SummarizedExperiment::assay(sce))),
               t(as.matrix(expected_matrix(e$matrices$ATAC_counts))))
  expect_true("tiles" %in% SingleCellExperiment::altExpNames(sce))
  tl <- SingleCellExperiment::altExp(sce, "tiles")
  expect_equal(unname(as.matrix(SummarizedExperiment::assay(tl))),
               t(as.matrix(expected_matrix(e$matrices$tiles_counts))))
  expect_true(all(grepl("^chr[12]:[0-9]+-[0-9]+$", rownames(tl))))
})

test_that("a layout this version does not know is refused, not read as empty", {
  p <- fixture_copy("piaso_atac")
  con <- DBI::dbConnect(RSQLite::SQLite(), p)
  DBI::dbExecute(con, "UPDATE matrix_meta SET layout = 'coo' WHERE matrix_name = 'ATAC_counts'")
  DBI::dbExecute(con, "DELETE FROM matrix_csc_chunks WHERE matrix_name = 'tiles_counts'")
  DBI::dbDisconnect(con)
  x <- cytome_open(p); on.exit(cytome_close(x))
  expect_error(read_cytome_matrix(x, "ATAC_counts"), "layout 'coo'")
  expect_error(read_cytome_matrix(x, "tiles_counts"), "records 392 nonzeros but has no column chunks")
  expect_error(cytome_stream(x, "tiles_counts", identity, along = "features"), "no column chunks")
})
