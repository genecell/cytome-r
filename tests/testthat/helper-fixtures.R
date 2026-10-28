# Stores written by the Python package and by PIASO's Rust importer, with the values R must
# read back (tools/make_python_fixtures.py writes both). The stores are kept gzip-compressed:
# each carries every table of the schema, close to a megabyte unpacked and a few tens of
# kilobytes packed.

.fixture_dir <- function() test_path("fixtures")

# The store `name`, unpacked once per session into a temporary folder.
fixture_store <- function(name, dir = file.path(tempdir(), "cytome_fixtures")) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  out <- file.path(dir, paste0(name, ".cytome"))
  if (!file.exists(out)) {
    src <- file.path(.fixture_dir(), paste0(name, ".cytome.gz"))
    if (!file.exists(src)) skip(paste("fixture", name, "not present"))
    con <- gzfile(src, "rb")
    chunks <- list()
    repeat {
      b <- readBin(con, "raw", 1048576L)
      if (!length(b)) break
      chunks[[length(chunks) + 1L]] <- b
    }
    close(con)
    writeBin(unlist(chunks, use.names = FALSE), out)
  }
  out
}

# A fresh copy of the store, for tests that change its folder or its files.
fixture_copy <- function(name, dir = tempfile("fixture_")) {
  dir.create(dir, recursive = TRUE)
  out <- file.path(dir, paste0(name, ".cytome"))
  file.copy(fixture_store(name), out)
  out
}

fixture_expected <- function(name) {
  jsonlite::fromJSON(file.path(.fixture_dir(), paste0(name, ".json")), simplifyVector = TRUE)
}

# An expected sparse matrix (cells x features, as stored), from its 0-based triplets.
expected_matrix <- function(e) {
  Matrix::sparseMatrix(i = as.integer(e$i) + 1L, j = as.integer(e$j) + 1L,
                       x = as.numeric(e$x), dims = as.integer(e$shape))
}

# An expected dense matrix, from its row-major hex values.
expected_dense <- function(e) {
  matrix(as.numeric(e$x), nrow = e$shape[1], ncol = e$shape[2], byrow = TRUE)
}
