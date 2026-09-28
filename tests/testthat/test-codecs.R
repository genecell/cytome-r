# Without the optional zstandard package, the Python writer compresses with
# zlib and still labels the blob "zstd" (cytome/io/compression.py falls back
# silently). The Python reader decides by the bytes; the R reader trusted the
# label and failed with "zstd frame content size unknown". It now decides by
# the bytes too, exactly when the label is "zstd".

.bytes <- function() as.raw(rep(0:255, 20))

test_that("zlib bytes labelled zstd are read", {
  z <- memCompress(.bytes(), type = "gzip")              # a zlib stream, 78 9c
  expect_identical(z[1], as.raw(0x78))
  expect_identical(cytome:::cytome_decompress(z, "zstd"), .bytes())
})

test_that("lz4 bytes labelled zstd are read", {
  l <- cytome:::cytome_compress(.bytes(), "lz4")
  expect_identical(cytome:::cytome_decompress(l, "zstd"), .bytes())
})

test_that("real zstd frames, and the other labels, are unchanged", {
  for (m in c("zstd", "lz4", "zlib")) {
    b <- cytome:::cytome_compress(.bytes(), m)
    expect_identical(cytome:::cytome_decompress(b, m), .bytes(), info = m)
  }
})

test_that("an uncompressed blob passes through", {
  expect_identical(cytome:::cytome_decompress(.bytes(), "none"), .bytes())
})

test_that("a whole store written without zstandard reads", {
  p <- system.file("extdata", "zlib_as_zstd.cytome", package = "cytome")
  skip_if(!nzchar(p))
  x <- cytome_open(p); on.exit(cytome_close(x))
  m <- read_cytome_matrix(x, "RNA_counts")
  e <- utils::read.csv(system.file("extdata", "expected_zlib_as_zstd.csv", package = "cytome"))
  expect_equal(dim(m), c(50L, 120L))
  expect_equal(unname(as.matrix(m)),
               unname(as.matrix(Matrix::sparseMatrix(i = e$i, j = e$j, x = e$x, dims = c(50L, 120L)))),
               tolerance = 1e-6)
  expect_equal(dim(cytome_embedding(x, "RNA_umap")), c(120L, 2L))
})
