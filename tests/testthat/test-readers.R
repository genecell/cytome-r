# What read_cytome() and cytome_obs() bring across from a store written by Python's
# from_anndata: column types, layers of the right shape only, cell embeddings only, obs
# graphs only (audit part 4, major 21).

.check_column <- function(got, e, info) {
  if (e$r == "factor") {
    expect_s3_class(got, "factor")
    expect_equal(is.ordered(got), isTRUE(e$ordered), info = info)
    expect_equal(levels(got), as.character(e$categories), info = info)
    want <- e$values
    if (e$kind == "numeric") expect_equal(as.numeric(as.character(got)), as.numeric(want), info = info)
    else if (e$kind == "logical") expect_equal(as.logical(as.character(got)), as.logical(want), info = info)
    else expect_equal(as.character(got), as.character(want), info = info)
  } else if (e$r == "logical") {
    expect_type(got, "logical"); expect_equal(got, as.logical(e$values), info = info)
  } else if (e$r == "integer") {
    expect_type(got, "integer"); expect_equal(got, as.integer(e$values), info = info)
  } else if (e$r == "double") {
    expect_type(got, "double"); expect_equal(got, as.numeric(e$values), tolerance = 1e-6, info = info)
  } else {
    expect_type(got, "character"); expect_equal(got, as.character(e$values), info = info)
  }
}

test_that("obs and var columns come back in the types their writer recorded", {
  p <- fixture_store("py_types"); e <- fixture_expected("py_types")
  x <- cytome_open(p); on.exit(cytome_close(x))
  obs <- cytome_obs(x)
  expect_equal(obs$barcode, e$barcodes)
  for (col in names(e$obs)) {
    expect_true(col %in% names(obs), info = col)               # "col/with slash" keeps its name
    .check_column(obs[[col]], e$obs[[col]], col)
  }
  var <- cytome_var(x, "RNA_counts")
  for (col in names(e$var)) .check_column(var[[col]], e$var[[col]], col)
})

test_that("a value outside the recorded categories is kept as an extra level", {
  p <- fixture_copy("py_types")
  con <- DBI::dbConnect(RSQLite::SQLite(), p)
  DBI::dbExecute(con, "UPDATE cells SET cluster = 'z' WHERE cell_idx = 0")      # _column_meta not updated
  DBI::dbExecute(con, "UPDATE cells SET flag = 7 WHERE cell_idx = 0")
  DBI::dbDisconnect(con)
  x <- cytome_open(p); on.exit(cytome_close(x))
  expect_warning(obs <- cytome_obs(x), "cells.cluster are not among")
  expect_equal(as.character(obs$cluster[1]), "z")
  expect_equal(utils::tail(levels(obs$cluster), 1), "z")
  expect_true(is.numeric(obs$flag))                               # left as stored
})

test_that("only matrices of the main assay's shape become layers; the rest are named", {
  skip_if_not_installed("SingleCellExperiment")
  p <- fixture_store("py_types"); e <- fixture_expected("py_types")
  expect_message(sce <- read_cytome(p), "RNA_wide \\(23 x 600\\).*RNA_raw_X \\(23 x 10\\)|RNA_raw_X \\(23 x 10\\).*RNA_wide \\(23 x 600\\)")
  expect_setequal(SummarizedExperiment::assayNames(sce), c("counts", "log1p", "counts_int", "zero"))
  expect_equal(unname(as.matrix(SummarizedExperiment::assay(sce, "counts_int"))),
               t(as.matrix(expected_matrix(e$matrices$RNA_counts_int))))
  expect_equal(rownames(sce), e$var$gene_id$values)
  rd <- SingleCellExperiment::reducedDimNames(sce)
  expect_setequal(rd, names(e$embeddings))                       # no var embedding
  for (nm in names(e$embeddings))
    expect_equal(unname(SingleCellExperiment::reducedDim(sce, nm)), expected_dense(e$embeddings[[nm]]),
                 tolerance = 1e-6, info = nm)
  obs_graphs <- names(e$graphs)[vapply(e$graphs, function(g) g$axis == "obs", logical(1))]
  expect_setequal(SingleCellExperiment::colPairNames(sce), obs_graphs)
  expect_s3_class(SummarizedExperiment::colData(sce)$cluster, "factor")
})

test_that("graphs over features are not listed as cell graphs", {
  p <- fixture_copy("py_types")
  con <- DBI::dbConnect(RSQLite::SQLite(), p)
  DBI::dbExecute(con, "INSERT INTO graph_edges (graph_name, axis, entity_table, row_idx, col_idx, value)
                       VALUES ('gene_net', 'var', 'genes', 0, 1, 0.5)")
  DBI::dbDisconnect(con)
  x <- cytome_open(p); on.exit(cytome_close(x))
  expect_false(any(c("gene_net", "RNA_varp_corr") %in% cytome_graphs(x)))
  expect_error(cytome_graph(x, "gene_net"), "not found")
})

test_that("Seurat gets the layers under the names Seurat gives the features", {
  skip_if_not_installed("SeuratObject")
  p <- fixture_store("py_types"); e <- fixture_expected("py_types")
  so <- suppressMessages(read_cytome(p, as = "Seurat"))
  expect_true(all(c("counts", "log1p", "counts_int", "zero") %in% SeuratObject::Layers(so[["RNA"]])))
  expect_equal(rownames(so), gsub("_", "-", e$var$gene_id$values))
  expect_true(all(names(e$embeddings) %in% SeuratObject::Reductions(so)))
  expect_s3_class(so$cluster, "factor")
  lg <- SeuratObject::LayerData(so, assay = "RNA", layer = "log1p")
  expect_equal(unname(as.matrix(lg)), t(as.matrix(expected_matrix(e$matrices$RNA_log1p))), tolerance = 1e-6)
})

test_that("delayed = TRUE backs every assay with the file, layers and altExps too", {
  skip_if_not_installed("DelayedArray")
  skip_if_not_installed("SingleCellExperiment")
  sce <- read_cytome(fixture_store("py_10x"), delayed = TRUE)
  expect_s4_class(SummarizedExperiment::assay(sce), "DelayedArray")
  expect_s4_class(SummarizedExperiment::assay(SingleCellExperiment::altExp(sce, "ATAC")), "DelayedArray")
  st <- suppressMessages(read_cytome(fixture_store("py_types"), delayed = TRUE))
  expect_s4_class(SummarizedExperiment::assay(st, "log1p"), "DelayedArray")
})

# cytome 0.3.9 stores a categorical's values in the type of its categories (INTEGER, REAL or
# TEXT; never codes) and records nullable booleans as "boolean"; earlier versions stored every
# categorical as text, a float category as the 15 significant digits SQLite gives it. Both
# forms must read alike. The cell names follow `_anndata_obs_index` when Python recorded one
# (obs with its own barcode column, as after merging libraries whose barcodes repeat).
test_that("categoricals read alike in either storage, and cell names follow the recorded index", {
  p <- tempfile(fileext = ".cytome")
  con <- DBI::dbConnect(RSQLite::SQLite(), p)
  DBI::dbExecute(con, paste("CREATE TABLE cells (cell_idx INTEGER PRIMARY KEY, barcode TEXT NOT NULL,",
                            "_obs_name TEXT, i_new INTEGER, f_new REAL, b_new INTEGER,",
                            "f_old TEXT, b_old TEXT, nb INTEGER)"))
  DBI::dbExecute(con, "CREATE TABLE _column_meta (table_name TEXT, column_name TEXT, dtype TEXT, categories TEXT)")
  DBI::dbExecute(con, "CREATE TABLE _metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL, value_type TEXT NOT NULL)")
  DBI::dbExecute(con, "INSERT INTO cells VALUES (0, 'BC0', 'a_BC0', 1, 0.30000000000000004, 1, '0.3', 'True', 1),
                                                 (1, 'BC0', 'b_BC0', 4, NULL, 0, '1.0', '0', NULL)")
  meta <- data.frame(table_name = "cells",
                     column_name = c("i_new", "f_new", "b_new", "f_old", "b_old", "nb"),
                     dtype = c(rep("categorical", 5), "boolean"),
                     categories = c("[1, 2, 3, 4]", "[0.30000000000000004, 1.0]", "[false, true]",
                                    "[0.30000000000000004, 1.0]", "[false, true]", NA))
  DBI::dbWriteTable(con, "_column_meta", meta, append = TRUE)
  DBI::dbExecute(con, "INSERT INTO _metadata VALUES ('_anndata_obs_index', '\"_obs_name\"', 'scalar')")
  DBI::dbDisconnect(con)
  x <- cytome_open(p); on.exit(cytome_close(x))
  o <- cytome_obs(x)
  expect_equal(levels(o$i_new), c("1", "2", "3", "4")); expect_equal(as.integer(o$i_new), c(1L, 4L))
  expect_equal(as.integer(o$f_new), c(1L, NA)); expect_equal(as.integer(o$f_old), c(1L, 2L))
  expect_equal(levels(o$b_new), c("FALSE", "TRUE"))
  expect_equal(as.character(o$b_new), c("TRUE", "FALSE")); expect_equal(as.character(o$b_old), c("TRUE", "FALSE"))
  expect_identical(o$nb, c(TRUE, NA))
  expect_equal(cytome:::.cell_names(x, o), c("a_BC0", "b_BC0"))
})

test_that("repeated barcodes become unique names in Seurat, which requires them", {
  skip_if_not_installed("SeuratObject")
  p <- fixture_copy("py_10x")
  con <- DBI::dbConnect(RSQLite::SQLite(), p)
  DBI::dbExecute(con, "UPDATE cells SET barcode = 'AAAC-1' WHERE cell_idx < 3")
  DBI::dbDisconnect(con)
  expect_warning(so <- suppressMessages(read_cytome(p, as = "Seurat")), "repeat")
  expect_false(anyDuplicated(colnames(so)) > 0)
})
