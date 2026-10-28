#' cytome: native R reader for the Cytome single-cell store
#'
#' A `.cytome` file is a plain SQLite database. `cytome` opens it read-only with DBI,
#' reads the relational tables (cells / genes / peaks, embeddings, graphs, fragments)
#' directly, and decodes the chunked compressed CSR count matrices via a small Rcpp
#' layer (lz4-block / zstd / zlib) — with zero Python dependency.
#'
#' @name cytome
#' @useDynLib cytome, .registration = TRUE
#' @importFrom Rcpp evalCpp
NULL

#' Open a .cytome file (read-only)
#'
#' Nothing is written, to the file or to its folder beyond the `-wal` and `-shm` files SQLite
#' keeps beside a store in WAL mode where the folder lets it. A store in WAL mode in a folder
#' that is not writable (an installed package, a read-only share) is opened immutable when its
#' WAL holds nothing, which needs no `-shm`; when its WAL still holds committed changes it
#' cannot be read there, and the error says so. As in the Python package, readers wait up to
#' 60 s for a writer's lock instead of failing at once.
#'
#' @param path Path to a `.cytome` SQLite file.
#' @return A `cytome` handle (a list with the open DBI connection). Close with [cytome_close()].
#' @export
cytome_open <- function(path) {
  if (!file.exists(path)) stop("cytome: file not found: ", path)
  path <- normalizePath(path, mustWork = TRUE)
  ## synchronous = NULL: RSQLite would otherwise set it, a write a read-only connection
  ## refuses with a warning. bigint = "numeric": counts past 2^31 come back as doubles.
  con <- DBI::dbConnect(RSQLite::SQLite(), .read_only_uri(path), flags = RSQLite::SQLITE_RO,
                        synchronous = NULL, bigint = "numeric")
  DBI::dbExecute(con, "PRAGMA busy_timeout = 60000")
  obj <- list(con = con, path = path)
  class(obj) <- "cytome"
  obj
}

## The URI that reads `path` without writing to it or its folder, chosen as the Python
## package chooses it (sqlite_engine._read_only_uri): `mode=ro` when the file is not in WAL
## mode, or SQLite can use its -wal and -shm (they exist, or the folder lets SQLite create
## them); otherwise `immutable=1`, which takes no locks and ignores the WAL, and is therefore
## used only when the WAL is absent or empty.
.read_only_uri <- function(path) {
  wal <- paste0(path, "-wal"); shm <- paste0(path, "-shm")
  if (!.in_wal_mode(path) || file.access(dirname(path), 2L) == 0L ||
      (file.exists(wal) && file.exists(shm)))
    return(.sqlite_uri(path, "mode=ro"))
  wal_bytes <- if (file.exists(wal)) file.size(wal) else 0
  if (wal_bytes > 0)
    stop("cytome: cannot read ", path, ": ", format(wal_bytes, big.mark = ","),
         " bytes of committed changes are still in ", basename(wal), ", and SQLite cannot ",
         "open it here because the folder ", dirname(path), " is not writable. Checkpoint or ",
         "close the store where it was written, or copy all of ", basename(path),
         "* to a writable folder and open the copy.")
  .sqlite_uri(path, "mode=ro&immutable=1")
}

## Bytes 18 and 19 of an SQLite header are its write and read versions: 2 means WAL.
.in_wal_mode <- function(path) {
  con <- file(path, "rb"); on.exit(close(con))
  h <- readBin(con, "raw", 20L)
  length(h) == 20L && as.integer(h[19L]) == 2L
}

## A file: URI for SQLite, the path percent-encoded (a '?', '#' or '%' in a folder name
## would otherwise end the path) and a Windows drive letter placed after the authority.
.sqlite_uri <- function(path, query) {
  p <- gsub("\\\\", "/", path)
  if (grepl("^[A-Za-z]:", p)) p <- paste0("/", p)
  b <- as.integer(charToRaw(enc2utf8(p)))
  keep <- (b >= 48L & b <= 57L) | (b >= 65L & b <= 90L) | (b >= 97L & b <= 122L) |
    b %in% as.integer(charToRaw("._~/:-"))
  ch <- character(length(b))
  ch[keep] <- vapply(b[keep], function(v) rawToChar(as.raw(v)), "")
  ch[!keep] <- sprintf("%%%02X", b[!keep])
  paste0("file://", paste(ch, collapse = ""), "?", query)
}

#' Close a cytome handle
#' @param x A `cytome` handle from [cytome_open()].
#' @return `NULL`, invisibly.
#' @export
cytome_close <- function(x) {
  if (inherits(x, "cytome") && DBI::dbIsValid(x$con)) DBI::dbDisconnect(x$con)
  invisible(NULL)
}

.tables <- function(x) DBI::dbListTables(x$con)
.has_table <- function(x, t) t %in% .tables(x)

#' @export
print.cytome <- function(x, ...) {
  m <- cytome_matrices(x)
  cat("<cytome>", x$path, "\n")
  cat("  cells:", tryCatch(DBI::dbGetQuery(x$con, "SELECT count(*) n FROM cells")$n, error = function(e) NA), "\n")
  if (nrow(m)) {
    cat("  matrices:\n")
    for (i in seq_len(nrow(m)))
      cat(sprintf("    %-16s %d x %d (%s, %s)\n", m$matrix_name[i], m$n_rows[i], m$n_cols[i],
                  m$dtype[i], paste0(m$row_entity[i], " x ", m$col_entity[i])))
  }
  emb <- cytome_embeddings(x)
  if (length(emb)) cat("  embeddings:", paste(emb, collapse = ", "), "\n")
  invisible(x)
}

#' Matrix catalogue (from `matrix_meta`)
#' @param x A `cytome` handle.
#' @return A data.frame: matrix_name, n_rows, n_cols, dtype, row_entity, col_entity.
#' @export
cytome_matrices <- function(x) {
  if (!.has_table(x, "matrix_meta")) return(data.frame())
  DBI::dbGetQuery(x$con,
    "SELECT matrix_name, n_rows, n_cols, dtype, row_entity, col_entity FROM matrix_meta")
}

#' High-level summary of a cytome
#' @param x A `cytome` handle.
#' @return A list: `path`, `n_cells`, `matrices` (as [cytome_matrices()]), `embeddings`
#'   (names) and `tables` (the SQLite tables of the store).
#' @export
cytome_info <- function(x) {
  list(path = x$path,
       n_cells = DBI::dbGetQuery(x$con, "SELECT count(*) n FROM cells")$n,
       matrices = cytome_matrices(x),
       embeddings = cytome_embeddings(x),
       tables = .tables(x))
}

#' Read an entity table (cells / genes / peaks / tiles / samples) as a data.frame
#'
#' Rows come in their stored order, column names as stored, and each column in the type its
#' writer recorded (`_column_meta`): a categorical as a factor with its levels in order, a
#' boolean as logical, an integer as integer.
#'
#' @param x A `cytome` handle.
#' @param entity Entity table name. Default `"cells"` (obs). Use `"genes"`/`"peaks"` for var.
#' @return A data.frame.
#' @export
cytome_obs <- function(x, entity = "cells") {
  if (!.has_table(x, entity)) stop("cytome: no entity table '", entity, "'")
  qt <- DBI::dbQuoteIdentifier(x$con, entity)
  info <- DBI::dbGetQuery(x$con, paste0("PRAGMA table_info(", qt, ")"))
  pk <- info$name[info$pk == 1]
  q <- paste("SELECT * FROM", qt)
  if (length(pk) == 1L) q <- paste(q, "ORDER BY", DBI::dbQuoteIdentifier(x$con, pk))
  .restore_column_types(x, entity, DBI::dbGetQuery(x$con, q))
}

#' Read the feature (var) table for a modality's matrix
#'
#' Resolves the matrix's `col_entity` (e.g. `RNA_counts` -> `genes`, `ATAC_counts` -> `peaks`).
#' @param x A `cytome` handle.
#' @param matrix_name Matrix name, e.g. `"RNA_counts"`. If `NULL`, returns the `genes` table.
#' @return A data.frame of features.
#' @export
cytome_var <- function(x, matrix_name = NULL) {
  if (is.null(matrix_name)) return(cytome_obs(x, "genes"))
  ce <- DBI::dbGetQuery(x$con,
    "SELECT col_entity FROM matrix_meta WHERE matrix_name = ?", params = list(matrix_name))
  if (!nrow(ce)) stop("cytome: matrix '", matrix_name, "' not found")
  cytome_obs(x, ce$col_entity[1])
}

# Column names of a table, or character(0) when it does not exist.
.table_columns <- function(x, table) {
  if (!.has_table(x, table)) return(character(0))
  DBI::dbGetQuery(x$con, paste0("PRAGMA table_info(", DBI::dbQuoteIdentifier(x$con, table), ")"))$name
}

# The identifiers of an entity's rows, in row order: barcodes for cells, the id column of a
# feature table.
.entity_ids <- function(x, entity) {
  if (entity == "cells") {
    return(as.character(DBI::dbGetQuery(x$con,
      "SELECT barcode FROM cells ORDER BY cell_idx")$barcode))
  }
  .feature_ids(cytome_obs(x, entity))
}
