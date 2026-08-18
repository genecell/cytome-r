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
#' @param path Path to a `.cytome` SQLite file.
#' @return A `cytome` handle (a list with the open DBI connection). Close with [cytome_close()].
#' @export
cytome_open <- function(path) {
  if (!file.exists(path)) stop("cytome: file not found: ", path)
  con <- DBI::dbConnect(RSQLite::SQLite(), path, flags = RSQLite::SQLITE_RO)
  obj <- list(con = con, path = normalizePath(path))
  class(obj) <- "cytome"
  obj
}

#' Close a cytome handle
#' @param x A `cytome` handle from [cytome_open()].
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
#' @param x A `cytome` handle.
#' @param entity Entity table name. Default `"cells"` (obs). Use `"genes"`/`"peaks"` for var.
#' @return A data.frame.
#' @export
cytome_obs <- function(x, entity = "cells") {
  if (!.has_table(x, entity)) stop("cytome: no entity table '", entity, "'")
  DBI::dbReadTable(x$con, entity)
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
