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

# The lz4 size header is little-endian on every platform (python lz4.block writes it so),
# and the format caps a block at LZ4_MAX_INPUT_SIZE (about 2.1 GB). The header used to be
# read through the host's integer layout, and a size past the cap was first allocated, then
# failed inside the decoder.

test_that("the lz4 size header is read as little-endian", {
  payload <- as.raw(rep(c(1, 2, 3), 100))                        # 300 bytes
  l <- cytome:::cytome_compress(payload, "lz4")
  expect_identical(l[1:4], as.raw(c(0x2c, 0x01, 0x00, 0x00)))    # 300, low byte first
  expect_identical(cytome:::cytome_decompress(l, "lz4"), payload)
})

test_that("an lz4 header past the format's limit is refused before anything is allocated", {
  l <- cytome:::cytome_compress(as.raw(1:100), "lz4")
  l[1:4] <- as.raw(c(0xf0, 0xff, 0xff, 0xff))                    # says 4 GB
  expect_error(cytome:::cytome_decompress(l, "lz4"), "past the format's limit")
})

test_that("zlib streams of any compression ratio decode, and a truncated one is refused", {
  zeros <- raw(5e6)                                              # compresses about 1000-fold
  z <- cytome:::cytome_compress(zeros, "zlib")
  expect_lt(length(z), 1e5)
  expect_identical(cytome:::cytome_decompress(z, "zlib"), zeros)
  expect_identical(cytome:::cytome_decompress(c(z, as.raw(1:3)), "zlib"), zeros)   # trailing bytes ignored, as in Python
  expect_error(cytome:::cytome_decompress(z[seq_len(length(z) - 4)], "zlib"), "truncated")
  expect_error(cytome:::cytome_decompress(raw(0), "zlib"), "truncated")
})

test_that("empty payloads round-trip under every codec", {
  for (m in c("zstd", "lz4", "zlib"))
    expect_identical(cytome:::cytome_decompress(cytome:::cytome_compress(raw(0), m), m), raw(0), info = m)
})
