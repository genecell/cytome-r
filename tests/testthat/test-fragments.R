# Fragments in compressed chunks (fragment_chunks) -- what PIASO's importer and Python's
# streaming import write, and the only form subset() keeps -- were never read: a query
# returned 0 of 800 fragments. A query without an interval failed ("Query does not require
# parameters") on every store (audit part 4, major 21).

.frag_df <- function(e) {
  data.frame(start = as.integer(e$start), end = as.integer(e$end), cell_idx = as.integer(e$cell_idx))
}
.sorted <- function(df) {
  df <- df[order(df$start, df$end, df$cell_idx), c("start", "end", "cell_idx")]
  rownames(df) <- NULL
  df
}

test_that("fragments in compressed chunks are read, whole and by interval", {
  p <- fixture_store("piaso_atac"); e <- fixture_expected("piaso_atac")
  x <- cytome_open(p); on.exit(cytome_close(x))
  n_chunks <- DBI::dbGetQuery(x$con, "SELECT COUNT(*) n FROM fragment_chunks WHERE chrom = 'chr1'")$n
  expect_gt(n_chunks, 1)
  for (chrom in c("chr1", "chr2")) {
    want <- .sorted(.frag_df(e$fragments[[chrom]]))
    got <- cytome_fragments(x, chrom)
    expect_equal(nrow(got), nrow(want))
    expect_equal(.sorted(got), want, info = chrom)
    expect_true(all(is.na(got$dup_count)))
    expect_equal(got$barcode, e$barcodes[got$cell_idx + 1L])
    sub <- cytome_fragments(x, chrom, 3000, 6000)
    keep <- want$start < 6000 & want$end > 3000
    expect_gt(sum(keep), 0)
    expect_equal(.sorted(sub), .sorted(want[keep, ]), info = chrom)
  }
})

test_that("fragments in the per-row tables are read without an interval too", {
  p <- fixture_store("py_frag_rows"); e <- fixture_expected("py_frag_rows")
  x <- cytome_open(p); on.exit(cytome_close(x))
  want <- .frag_df(e$fragments$chr1)
  got <- cytome_fragments(x, "chr1")
  expect_equal(.sorted(got), .sorted(want))
  expect_equal(got$dup_count[order(got$start, got$end, got$cell_idx)],
               as.integer(e$fragments$chr1$dup_count)[order(want$start, want$end, want$cell_idx)])
  sub <- cytome_fragments(x, "chr1", 2000, 4000)
  expect_equal(.sorted(sub), .sorted(want[want$start < 4000 & want$end > 2000, ]))
})

test_that("chunks stored with encoding 0 (absolute starts and ends) read the same", {
  p <- fixture_copy("piaso_atac"); e <- fixture_expected("piaso_atac")
  con <- DBI::dbConnect(RSQLite::SQLite(), p)
  ch <- DBI::dbGetQuery(con, "SELECT id, starts_blob, ends_blob, cell_idx_blob, compression, encoding
                              FROM fragment_chunks WHERE chrom = 'chr2'")
  for (k in seq_len(nrow(ch))) {
    f <- cytome:::cytome_decode_fragments(ch$starts_blob[[k]], ch$ends_blob[[k]], ch$cell_idx_blob[[k]],
                                          ch$compression[k], ch$encoding[k])
    blob <- function(v) cytome:::cytome_compress(writeBin(v, raw(), size = 4, endian = "little"), ch$compression[k])
    DBI::dbExecute(con, "UPDATE fragment_chunks SET starts_blob = ?, ends_blob = ?, encoding = 0 WHERE id = ?",
                   params = list(list(blob(f$start)), list(blob(f$end)), ch$id[k]))
  }
  DBI::dbDisconnect(con)
  x <- cytome_open(p); on.exit(cytome_close(x))
  expect_equal(.sorted(cytome_fragments(x, "chr2")), .sorted(.frag_df(e$fragments$chr2)))
})

test_that("an unknown chromosome names the ones there are; an empty one is empty", {
  p <- fixture_store("piaso_atac")
  x <- cytome_open(p); on.exit(cytome_close(x))
  expect_error(cytome_fragments(x, "chrUn_x"), "chr1, chr2")
  none <- cytome_fragments(x, "chr5")                          # a table, no fragments
  expect_equal(nrow(none), 0L)
  expect_named(none, c("chrom", "start", "end", "cell_idx", "dup_count", "barcode"))
  expect_error(cytome_fragments(x, "chr1", start = 10), "both start and end")
})
