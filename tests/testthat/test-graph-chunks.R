# Graphs written by the Python package from cytome 0.3.6 are stored as
# compressed row chunks (graph_meta / graph_chunks), not one row per edge.
# 0.1.0 read only graph_edges, so these graphs were silently absent: an empty
# cytome_graphs() and a Seurat / SCE object without its neighbour graph.
# The fixture is written by tools/make_graph_fixture.py with Python cytome
# 0.3.6: 60 cells, two graphs in 7 chunks each, the last cell without edges.

.fixture <- function() system.file("extdata", "graph_chunks.cytome", package = "cytome")

.expected <- function(name, n = 60L) {
  # values are hex floats, which R parses exactly on every platform
  e <- utils::read.csv(system.file("extdata", paste0("expected_graph_", name, ".csv"),
                                   package = "cytome"),
                       colClasses = c("integer", "integer", "character"))
  Matrix::sparseMatrix(i = e$i, j = e$j, x = as.numeric(e$x), dims = c(n, n))
}

test_that("chunked graphs are listed", {
  x <- cytome_open(.fixture()); on.exit(cytome_close(x))
  expect_setequal(cytome_graphs(x), c("connectivities", "distances"))
})

test_that("a chunked graph reads back entry for entry, float32 and float64", {
  x <- cytome_open(.fixture()); on.exit(cytome_close(x))
  g <- cytome_graph(x, "connectivities")
  expect_s4_class(g, "dgCMatrix")
  expect_equal(dim(g), c(60L, 60L))
  expect_equal(g, .expected("connectivities"), tolerance = 1e-6)   # float32 on disk
  expect_equal(cytome_graph(x, "distances"), .expected("distances"), tolerance = 0)
})

test_that("the shape is kept, including a last cell without edges", {
  x <- cytome_open(.fixture()); on.exit(cytome_close(x))
  g <- cytome_graph(x, "connectivities")
  expect_equal(nrow(g), 60L)
  expect_equal(sum(g[60, ] != 0), 0)
})

test_that("n pads a chunked graph and refuses to cut it", {
  x <- cytome_open(.fixture()); on.exit(cytome_close(x))
  expect_equal(dim(cytome_graph(x, "connectivities", n = 70L)), c(70L, 70L))
  expect_error(cytome_graph(x, "connectivities", n = 50L), "larger than")
})

test_that("an unknown graph says which graphs there are", {
  x <- cytome_open(.fixture()); on.exit(cytome_close(x))
  expect_error(cytome_graph(x, "knn"), "connectivities")
})

test_that("chunked graphs reach Seurat and SingleCellExperiment", {
  skip_if_not_installed("SeuratObject")
  so <- read_cytome(.fixture(), as = "Seurat")
  expect_true(all(c("connectivities", "distances") %in% names(so@graphs)))
  skip_if_not_installed("SingleCellExperiment")
  sce <- read_cytome(.fixture())
  expect_true(all(c("connectivities", "distances") %in%
                  SingleCellExperiment::colPairNames(sce)))
})
