# Opening for reading writes nothing. A store in WAL mode in a folder that is not writable
# could not be opened ("attempt to write a readonly database"), and the package's own
# reference.cytome was in WAL mode, so its examples wrote -wal and -shm files into the
# installed library (audit part 4, major 27).

.sidecars <- function(path) file.exists(paste0(path, c("-wal", "-shm", "-journal")))

test_that("the bundled reference store is not in WAL mode, and reading it leaves nothing", {
  ref <- system.file("extdata", "reference.cytome", package = "cytome")
  skip_if(!nzchar(ref))
  h <- readBin(ref, "raw", 20L)
  expect_equal(as.integer(h[19:20]), c(1L, 1L))                 # rollback journal, not WAL
  before <- .sidecars(ref)
  x <- cytome_open(ref)
  invisible(read_cytome_matrix(x, "RNA_counts"))
  cytome_close(x)
  expect_equal(.sidecars(ref), before)
  expect_false(any(.sidecars(ref)))
})

.read_only_dir <- function(path) {
  Sys.chmod(dirname(path), "0555")
  # root, or a filesystem that ignores modes, can still write: nothing to test then
  if (file.access(dirname(path), 2L) == 0L) {
    Sys.chmod(dirname(path), "0755")
    skip("this folder stays writable (root, or modes ignored)")
  }
}

test_that("a WAL-mode store in a read-only folder opens immutable, and writes nothing", {
  skip_on_os("windows")
  p <- fixture_copy("piaso_atac"); e <- fixture_expected("piaso_atac")
  expect_equal(as.integer(readBin(p, "raw", 20L)[19:20]), c(2L, 2L))   # WAL mode, as PIASO writes it
  .read_only_dir(p)
  on.exit(Sys.chmod(dirname(p), "0755"), add = TRUE)
  x <- cytome_open(p)
  M <- read_cytome_matrix(x, "ATAC_counts", features_x_cells = FALSE)
  cytome_close(x)
  expect_equal(unname(as.matrix(M)), as.matrix(expected_matrix(e$matrices$ATAC_counts)))
  expect_false(any(.sidecars(p)))
})

test_that("committed changes left in the WAL of a read-only folder are refused, not read around", {
  skip_on_os("windows")
  src <- fixture_copy("piaso_atac")
  w <- DBI::dbConnect(RSQLite::SQLite(), src)
  DBI::dbExecute(w, "PRAGMA wal_autocheckpoint = 0")
  DBI::dbExecute(w, "UPDATE cells SET sample_id = 'S1'")
  dst <- file.path(tempfile("walcopy_"), "s.cytome")
  dir.create(dirname(dst))
  file.copy(src, dst); file.copy(paste0(src, "-wal"), paste0(dst, "-wal"))   # no -shm
  DBI::dbDisconnect(w)
  expect_gt(file.size(paste0(dst, "-wal")), 0)
  .read_only_dir(dst)
  on.exit(Sys.chmod(dirname(dst), "0755"), add = TRUE)
  expect_error(cytome_open(dst), "committed changes are still in s.cytome-wal")
})

test_that("readers wait for a writer's lock, and paths need no escaping", {
  d <- file.path(tempdir(), if (.Platform$OS.type == "windows") "a b%20c#d" else "a b%20c#d?e")   # no '?' on Windows
  dir.create(d, showWarnings = FALSE)
  p <- file.path(d, "x.cytome")
  file.copy(fixture_store("py_10x"), p, overwrite = TRUE)
  x <- cytome_open(p); on.exit(cytome_close(x))
  expect_equal(DBI::dbGetQuery(x$con, "PRAGMA busy_timeout")[[1]], 60000)
  expect_equal(nrow(cytome_obs(x)), fixture_expected("py_10x")$n_cells)
  expect_equal(x$path, normalizePath(p))
})
