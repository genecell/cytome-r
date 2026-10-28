# Native R writer for the .cytome SQLite format (no Python). Mirrors
# cytome/io/sqlite_engine.py (_create_schema), cytome/io/chunked_io.py
# (write_sparse_chunked), cytome/core/graph.py (GraphStore.write_sparse) and
# cytome/io/convert_anndata.py (_store_column_meta), so a store written here is
# read by the Python package as it reads its own.
#
# A store is built in a temporary file beside the destination, inside one
# transaction, and renamed over the destination only once it is complete: a
# write that fails leaves the old file as it was and no partial store. Values
# are placed by name, never by position: a matrix whose features or cells come
# in another order is written under its own features and cells.
#
# cytome never stores values that were not measured. A layer over only some of
# the features (a Seurat scale.data over the variable genes) is not written,
# with a warning that says how to recompute it after reading: stored over all
# the features, the zeros of the ones it lacks would read as values (0 in
# scale.data is the mean). An assay, layer or graph over only some of the cells
# (an assay measured in some batches only) is not written either, with a warning
# that says to subset the object to those cells and write that as its own store:
# its zeros would read as cells without counts, and a graph's as cells without
# neighbours. An embedding keeps NaN for the cells it lacks, which reads as
# missing.

## ---- modalities ------------------------------------------------------------

## Modality -> feature table, mirroring MODALITY_REGISTRY in
## cytome/utils/modality.py. These are the modalities the Python package
## knows; anything else (SCT, ADT, ERCC, ...) has no feature table there.
.MODALITY_TABLE <- c(RNA = "genes", GA = "GA_genes", ATAC = "peaks", tiles = "tiles")

## The canonical name of a modality (RNA, GA, ATAC, tiles), or NA.
.canonical_modality <- function(m) {
  hit <- names(.MODALITY_TABLE)[toupper(names(.MODALITY_TABLE)) == toupper(m)]
  if (length(hit)) hit else NA_character_
}

.modality_table <- function(modality) {
  m <- .canonical_modality(modality)
  if (is.na(m))
    stop("cytome: cannot write modality '", modality, "' from R. ",
         "Supported: ", paste(names(.MODALITY_TABLE), collapse = ", "),
         ". Build the cytome on the Python side for other modalities.")
  .MODALITY_TABLE[[m]]
}

## Which of `available` alternative modalities to write: all the writable ones
## when `requested` is NULL (the others named in a warning), else exactly those
## requested (an error for one that does not exist or cannot be written).
.resolve_alt <- function(requested, available, kind) {
  if (is.null(requested)) {
    bad <- available[is.na(vapply(available, .canonical_modality, ""))]
    if (length(bad))
      warning("cytome: not written, as cytome has no feature table for them: ",
              paste(bad, collapse = ", "), " (", kind, "). Writable modalities: ",
              paste(names(.MODALITY_TABLE), collapse = ", "), ".", call. = FALSE)
    return(setdiff(available, bad))
  }
  miss <- setdiff(requested, available)
  if (length(miss))
    stop("cytome: no such ", kind, ": ", paste(miss, collapse = ", "),
         ". Available: ", if (length(available)) paste(available, collapse = ", ") else "none")
  for (m in requested) .modality_table(m)
  requested
}

## `layers`: FALSE (none, the default), TRUE (all available), or a character
## vector naming the ones to carry. Returns the names to actually write.
.resolve_layers <- function(layers, available) {
  if (isFALSE(layers) || is.null(layers)) return(character(0))
  if (isTRUE(layers)) return(available)
  miss <- setdiff(layers, available)
  if (length(miss))
    stop("cytome: no such layer(s): ", paste(miss, collapse = ", "),
         ". Available: ", if (length(available)) paste(available, collapse = ", ")
                          else "none")
  intersect(layers, available)
}

## Detect a gene-symbol column in a per-feature annotation frame (SCE rowData /
## Seurat meta.features) so the cytome `symbol` column is populated with real
## symbols rather than a copy of the row ids — lets Python name resolution
## (COSG/dotplot/inferGRN) show symbols. Returns a character vector or NULL.
.detect_feature_symbol <- function(meta, n) {
  if (is.null(meta) || ncol(meta) == 0L || nrow(meta) != n) return(NULL)
  aliases <- c("symbol", "Symbol", "gene_symbols", "gene_symbol",
               "gene_name", "feature_name", "symbols")
  hit <- intersect(aliases, colnames(meta))
  if (length(hit) == 0L) return(NULL)
  as.character(meta[[hit[1]]])
}

## Parse "chr1:100-200" / "chr1-100-200" into a data frame. Python's peaks
## table declares chr/start/end_ NOT NULL, so an unparseable id has to fail
## here rather than produce a file that Python refuses.
.parse_peak_ids <- function(ids) {
  m <- regmatches(ids, regexec("^(.+?)[:_-]([0-9]+)[-_]([0-9]+)$", ids))
  ok <- vapply(m, length, integer(1)) == 4L
  if (!all(ok))
    stop("cytome: ATAC feature ids must look like 'chr1:100-200'; ",
         sum(!ok), " did not, e.g. ", ids[which(!ok)[1]])
  data.frame(chr = vapply(m, `[`, character(1), 2),
             start = as.numeric(vapply(m, `[`, character(1), 3)),
             end_ = as.numeric(vapply(m, `[`, character(1), 4)),
             stringsAsFactors = FALSE)
}

## ---- aligning values by name ---------------------------------------------

## Any matrix-like object as a dgCMatrix (logical, pattern, symmetric, dense
## and DelayedArray inputs included).
.as_dgC <- function(M, what) {
  tryCatch({
    if (!methods::is(M, "CsparseMatrix")) M <- methods::as(M, "CsparseMatrix")
    M <- methods::as(methods::as(M, "generalMatrix"), "dMatrix")
    methods::as(M, "dgCMatrix")
  }, error = function(e) stop("cytome: cannot read ", what, " as a sparse matrix: ",
                              conditionMessage(e), call. = FALSE))
}

## Place the rows (features) of `M` (features x cells) at their names' positions
## in `ids`, so that a matrix whose features come in another order is stored by
## name. A matrix over only some of the features is not placed: NULL is
## returned, and the caller leaves it out.
.align_features <- function(M, ids, what) {
  rn <- rownames(M)
  if (is.null(rn) || identical(as.character(rn), as.character(ids))) {
    if (nrow(M) != length(ids))
      stop("cytome: ", what, " has ", nrow(M), " features and no names; the assay has ", length(ids))
    return(M)
  }
  if (anyDuplicated(rn)) stop("cytome: ", what, " has duplicated feature names")
  pos <- match(rn, ids)
  if (anyNA(pos))
    stop("cytome: ", what, " has features that are not the assay's, e.g. ", rn[which(is.na(pos))[1]])
  if (length(rn) < length(ids)) return(NULL)
  M <- methods::as(.as_dgC(M, what), "TsparseMatrix")
  Matrix::sparseMatrix(i = pos[M@i + 1L], j = M@j + 1L, x = M@x,
                       dims = c(length(ids), ncol(M)), dimnames = list(NULL, colnames(M)))
}

## Place the columns (cells) of `M` (features x cells) at their names' positions
## among `cells`, so that a matrix whose cells come in another order is stored
## by name. A matrix over only some of the cells is not placed: NULL is
## returned, and the caller leaves it out.
.align_cells <- function(M, cells, what) {
  cn <- colnames(M)
  if (is.null(cn) || identical(as.character(cn), cells)) {
    if (ncol(M) != length(cells))
      stop("cytome: ", what, " has ", ncol(M), " cells and no names; the object has ", length(cells))
    return(M)
  }
  pos <- match(cn, cells)
  if (anyNA(pos) || anyDuplicated(cn))
    stop("cytome: ", what, " has cells that are not the object's, or repeats one")
  if (length(cn) < length(cells)) return(NULL)
  M <- methods::as(.as_dgC(M, what), "TsparseMatrix")
  Matrix::sparseMatrix(i = M@i + 1L, j = pos[M@j + 1L], x = M@x, dims = c(nrow(M), length(cells)))
}

## The warning for a matrix or graph over only some of the cells, left out.
.some_cells_warning <- function(what, n, total)
  warning("cytome: ", what, " is not written: it covers only ", n, " of ", total, " cells. ",
          "To keep it, subset the object to those cells and write that as a separate store.",
          call. = FALSE)

## An embedding's rows in the order of `cells` (rows it lacks are NaN).
.align_embedding <- function(E, cells, what) {
  E <- as.matrix(E)
  storage.mode(E) <- "double"
  rn <- rownames(E)
  if (is.null(rn) || identical(as.character(rn), cells)) {
    if (nrow(E) != length(cells))
      stop("cytome: embedding ", what, " has ", nrow(E), " rows; the object has ", length(cells), " cells")
    return(E)
  }
  pos <- match(cells, rn)
  if (anyNA(pos))
    warning("cytome: embedding ", what, " has no row for ", sum(is.na(pos)),
            " cells; they are stored as NaN", call. = FALSE)
  out <- E[pos, , drop = FALSE]
  out[is.na(pos), ] <- NaN
  out
}

## A cells x cells graph in the order of `cells`, or NULL, with a warning, for a
## graph over only some of the cells: a store records a graph's shape, not the
## cells it covers, so every reader would take a cell it lacks for a cell
## without neighbours.
.align_graph <- function(G, cells, what) {
  G <- .as_dgC(G, paste("graph", what))
  dn <- dimnames(G)
  if (nrow(G) != ncol(G)) stop("cytome: graph ", what, " is not square")
  if (is.null(dn[[1]]) || identical(as.character(dn[[1]]), cells)) {
    if (nrow(G) != length(cells))
      stop("cytome: graph ", what, " is ", nrow(G), " x ", ncol(G), "; the object has ",
           length(cells), " cells")
    return(G)
  }
  pr <- match(dn[[1]], cells); pc <- match(dn[[2]], cells)
  if (anyNA(pr) || anyNA(pc)) stop("cytome: graph ", what, " has cells that are not the object's")
  if (nrow(G) < length(cells)) {
    .some_cells_warning(paste0("graph '", what, "'"), nrow(G), length(cells))
    return(NULL)
  }
  T <- methods::as(G, "TsparseMatrix")
  Matrix::sparseMatrix(i = pr[T@i + 1L], j = pc[T@j + 1L], x = T@x,
                       dims = c(length(cells), length(cells)))
}

## ---- the store ---------------------------------------------------------

## The tables of Python's schema that its readers and validate() expect,
## created empty where this writer leaves them so.
.cytome_write_schema <- function(con) {
  stmts <- c(
    "CREATE TABLE _manifest (key TEXT PRIMARY KEY, value TEXT);",
    paste0("CREATE TABLE _provenance (id INTEGER PRIMARY KEY AUTOINCREMENT, timestamp TEXT NOT NULL, ",
           "operation TEXT NOT NULL, package_name TEXT NOT NULL, package_version TEXT NOT NULL, python_version TEXT, ",
           "function_name TEXT NOT NULL, parameters TEXT NOT NULL, dependency_versions TEXT, input_objects TEXT, ",
           "output_objects TEXT, random_seed INTEGER, hardware_info TEXT, duration_seconds REAL, peak_memory_mb REAL, notes TEXT);"),
    "CREATE TABLE _schema_migrations (migration_id INTEGER PRIMARY KEY, description TEXT, applied_at TEXT, from_version TEXT, to_version TEXT);",
    "CREATE TABLE _metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL, value_type TEXT NOT NULL);",
    "CREATE TABLE _column_meta (table_name TEXT NOT NULL, column_name TEXT NOT NULL, dtype TEXT NOT NULL, categories TEXT, PRIMARY KEY (table_name, column_name));",
    "CREATE TABLE genes (gene_idx INTEGER PRIMARY KEY, gene_id TEXT NOT NULL UNIQUE, symbol TEXT, chr TEXT, start INTEGER, end_ INTEGER, biotype TEXT);",
    "CREATE INDEX idx_genes_symbol ON genes(symbol);",
    "CREATE TABLE GA_genes (gene_idx INTEGER PRIMARY KEY, gene_id TEXT NOT NULL UNIQUE, symbol TEXT, chr TEXT, start INTEGER, end_ INTEGER, biotype TEXT);",
    ## ATAC features. Python's schema declares chr/start/end_ NOT NULL, so a
    ## peak_id that does not parse as chr:start-end would be rejected there --
    ## .write_feature_table() parses it and errors here instead, where the user
    ## can still see which id was the problem.
    "CREATE TABLE peaks (peak_idx INTEGER PRIMARY KEY, peak_id TEXT NOT NULL, chr TEXT NOT NULL, start INTEGER NOT NULL, end_ INTEGER NOT NULL, annotation TEXT, nearest_gene TEXT, distance_to_tss INTEGER);",
    "CREATE TABLE tiles (tile_idx INTEGER PRIMARY KEY, tile_id TEXT NOT NULL, chr TEXT NOT NULL, start INTEGER NOT NULL, end_ INTEGER NOT NULL);",
    "CREATE TABLE samples (sample_idx INTEGER PRIMARY KEY, sample_id TEXT NOT NULL UNIQUE, batch TEXT, condition TEXT, donor TEXT, chemistry TEXT, n_cells INTEGER);",
    "CREATE TABLE proteins (protein_idx INTEGER PRIMARY KEY, protein_id TEXT NOT NULL UNIQUE, gene_mapping TEXT);",
    paste0("CREATE TABLE matrix_chunks (id INTEGER PRIMARY KEY, matrix_name TEXT NOT NULL, ",
           "chunk_idx INTEGER NOT NULL, row_start INTEGER NOT NULL, row_end INTEGER NOT NULL, n_nonzero INTEGER NOT NULL, ",
           "data_blob BLOB NOT NULL, indices_blob BLOB NOT NULL, indptr_blob BLOB NOT NULL, dtype TEXT NOT NULL, ",
           "compression TEXT DEFAULT 'zstd', UNIQUE(matrix_name, chunk_idx));"),
    "CREATE INDEX idx_matrix_chunks ON matrix_chunks(matrix_name, chunk_idx);",
    paste0("CREATE TABLE matrix_csc_chunks (id INTEGER PRIMARY KEY, matrix_name TEXT NOT NULL, ",
           "chunk_idx INTEGER NOT NULL, col_start INTEGER NOT NULL, col_end INTEGER NOT NULL, n_nonzero INTEGER NOT NULL, ",
           "data_blob BLOB NOT NULL, indices_blob BLOB NOT NULL, indptr_blob BLOB NOT NULL, dtype TEXT NOT NULL, ",
           "compression TEXT DEFAULT 'zstd', UNIQUE(matrix_name, chunk_idx));"),
    paste0("CREATE TABLE matrix_meta (matrix_name TEXT PRIMARY KEY, n_rows INTEGER NOT NULL, ",
           "n_cols INTEGER NOT NULL, n_nonzero INTEGER NOT NULL, dtype TEXT NOT NULL, row_entity TEXT NOT NULL, ",
           "col_entity TEXT NOT NULL, chunk_size INTEGER NOT NULL, n_chunks INTEGER NOT NULL, has_csc INTEGER DEFAULT 0, ",
           "csc_chunk_size INTEGER, csc_n_chunks INTEGER, created_at TEXT, provenance_id INTEGER, is_integer INTEGER);"),
    paste0("CREATE TABLE dense_chunks (id INTEGER PRIMARY KEY, array_name TEXT NOT NULL, ",
           "chunk_idx INTEGER NOT NULL, row_start INTEGER NOT NULL, row_end INTEGER NOT NULL, n_cols INTEGER NOT NULL, ",
           "data_blob BLOB NOT NULL, dtype TEXT NOT NULL, compression TEXT DEFAULT 'zstd', UNIQUE(array_name, chunk_idx));"),
    "CREATE INDEX idx_dense_chunks ON dense_chunks(array_name, chunk_idx);",
    paste0("CREATE TABLE embedding_meta (array_name TEXT PRIMARY KEY, n_rows INTEGER NOT NULL, ",
           "n_cols INTEGER NOT NULL, dtype TEXT NOT NULL, entity TEXT NOT NULL, chunk_size INTEGER NOT NULL, ",
           "n_chunks INTEGER NOT NULL, created_at TEXT, provenance_id INTEGER);"),
    "CREATE TABLE graph_knn (cell_i INTEGER NOT NULL, cell_j INTEGER NOT NULL, weight REAL NOT NULL, PRIMARY KEY (cell_i, cell_j));",
    paste0("CREATE TABLE graph_peak_gene (peak_idx INTEGER NOT NULL, gene_idx INTEGER NOT NULL, distance INTEGER NOT NULL, ",
           "correlation REAL, p_value REAL, p_adj REAL, provenance_id INTEGER, PRIMARY KEY (peak_idx, gene_idx));"),
    paste0("CREATE TABLE graph_edges (graph_name TEXT NOT NULL, axis TEXT NOT NULL DEFAULT 'obs', ",
           "entity_table TEXT NOT NULL DEFAULT 'cells', row_idx INTEGER NOT NULL, col_idx INTEGER NOT NULL, ",
           "value REAL NOT NULL, PRIMARY KEY (graph_name, axis, row_idx, col_idx));"),
    paste0("CREATE TABLE graph_meta (graph_name TEXT NOT NULL, axis TEXT NOT NULL DEFAULT 'obs', ",
           "entity_table TEXT NOT NULL DEFAULT 'cells', n_rows INTEGER NOT NULL, n_cols INTEGER NOT NULL, ",
           "n_nonzero INTEGER NOT NULL, dtype TEXT NOT NULL, n_chunks INTEGER NOT NULL, PRIMARY KEY (graph_name, axis));"),
    paste0("CREATE TABLE graph_chunks (graph_name TEXT NOT NULL, axis TEXT NOT NULL DEFAULT 'obs', ",
           "chunk_idx INTEGER NOT NULL, row_start INTEGER NOT NULL, row_end INTEGER NOT NULL, n_nonzero INTEGER NOT NULL, ",
           "data_blob BLOB NOT NULL, indices_blob BLOB NOT NULL, indptr_blob BLOB NOT NULL, compression TEXT NOT NULL, ",
           "PRIMARY KEY (graph_name, axis, chunk_idx));"))
  for (s in stmts) DBI::dbExecute(con, s)
}

.now <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
.json <- function(x, ...) as.character(jsonlite::toJSON(x, auto_unbox = TRUE, null = "null", digits = NA, ...))

## int32 values as little-endian bytes, whatever the host.
.int32_raw <- function(v) writeBin(as.integer(v), raw(), size = 4L, endian = "little")
## Non-negative values below 2^31 as little-endian int64 bytes.
.int64_raw <- function(v) writeBin(as.vector(rbind(as.integer(v), 0L)), raw(), size = 4L, endian = "little")

.value_raw <- function(x, dtype) {
  switch(dtype,
         float32 = writeBin(as.double(x), raw(), size = 4L, endian = "little"),
         float64 = writeBin(as.double(x), raw(), size = 8L, endian = "little"),
         int32 = .int32_raw(x))
}

.check_dtype <- function(dtype) {
  if (!(is.character(dtype) && length(dtype) == 1L && dtype %in% c("float32", "float64", "int32")))
    stop("cytome: dtype must be \"float32\", \"float64\" or \"int32\", not ", deparse(dtype))
  dtype
}

## The cells table: the barcodes, then every column of `obs` in its own SQL
## type, with `_column_meta` saying which R type it was. Names are kept as they
## are. A column named like one of the table's own (cell_idx, barcode,
## sample_id, n_fragments) fills that column when it means the same thing.
.write_cells <- function(con, barcodes, obs) {
  n <- length(barcodes)
  cols <- list(cell_idx = seq_len(n) - 1L, barcode = as.character(barcodes),
               sample_id = rep(NA_character_, n), n_fragments = rep(0L, n))
  sql_types <- c(cell_idx = "INTEGER PRIMARY KEY", barcode = "TEXT NOT NULL",
                 sample_id = "TEXT", n_fragments = "INTEGER DEFAULT 0")
  meta <- list(list("barcode", "string", NA_character_))
  skipped <- character(0)
  if (!is.null(obs) && ncol(obs)) {
    if (nrow(obs) != n) stop("cytome: the cell table has ", nrow(obs), " rows for ", n, " cells")
    for (k in seq_along(obs)) {
      v <- obs[[k]]
      cn <- names(obs)[k]
      if (is.na(cn) || !nzchar(cn)) cn <- paste0("V", k)
      tv <- .sql_column(v)
      if (is.null(tv)) { skipped <- c(skipped, cn); next }
      name <- cn
      if (tolower(cn) == "cell_idx") next                              # positional; not data
      if (tolower(cn) == "barcode") {
        if (identical(as.character(v), cols$barcode)) next
        name <- "barcode_obs"
      } else if (tolower(cn) == "sample_id") {
        cols$sample_id <- as.character(v)
        meta[[length(meta) + 1L]] <- list("sample_id", tv$dtype, tv$categories)
        next
      } else if (tolower(cn) == "n_fragments") {
        if (is.numeric(v) && all(is.na(v) | v == round(v))) { cols$n_fragments <- as.integer(v); next }
        name <- "n_fragments_obs"
      }
      ## SQLite column names are case-insensitive
      while (tolower(name) %in% tolower(names(cols))) name <- paste0(name, "_")
      if (name != cn)
        warning("cytome: column '", cn, "' is stored as '", name, "'", call. = FALSE)
      cols[[name]] <- tv$values
      sql_types[[name]] <- tv$sql
      meta[[length(meta) + 1L]] <- list(name, tv$dtype, tv$categories)
    }
  }
  if (length(skipped))
    warning("cytome: not written, as cytome cannot store their type: ",
            paste(skipped, collapse = ", "), call. = FALSE)
  q <- function(s) DBI::dbQuoteIdentifier(con, s)
  DBI::dbExecute(con, sprintf("CREATE TABLE cells (%s);",
                              paste(q(names(sql_types)), sql_types, collapse = ", ")))
  DBI::dbExecute(con, "CREATE INDEX idx_cells_sample ON cells(sample_id);")
  if (n)
    DBI::dbExecute(con, sprintf("INSERT INTO cells (%s) VALUES (%s)",
                                paste(q(names(cols)), collapse = ", "),
                                paste(rep("?", length(cols)), collapse = ", ")),
                   params = unname(cols))
  .write_column_meta(con, "cells", meta)
  ## One row per sample, as the format's samples table holds them.
  sid <- cols$sample_id[!is.na(cols$sample_id)]
  if (length(sid)) {
    tab <- table(factor(sid, levels = unique(sid)))
    DBI::dbExecute(con, "INSERT INTO samples (sample_idx, sample_id, n_cells) VALUES (?, ?, ?)",
                   params = list(seq_along(tab) - 1L, names(tab), as.integer(tab)))
  }
}

## One column as SQLite stores it and as `_column_meta` describes it, or NULL
## for a type cytome cannot store.
.sql_column <- function(v) {
  if (!is.null(dim(v)) || is.list(v) || is.complex(v) || is.raw(v)) return(NULL)
  if (is.factor(v)) {
    lv <- levels(v)
    return(list(values = as.character(v), sql = "TEXT", categories = .json(I(lv[!is.na(lv)])),   # addNA(): missing, not a category
                dtype = if (is.ordered(v)) "ordered_categorical" else "categorical"))
  }
  if (is.logical(v))
    return(list(values = as.integer(v), sql = "INTEGER", dtype = "bool", categories = NA_character_))
  if (inherits(v, "integer64"))
    return(list(values = as.numeric(v), sql = "INTEGER", dtype = "int64", categories = NA_character_))
  if (is.integer(v))
    return(list(values = as.integer(v), sql = "INTEGER", dtype = "int32", categories = NA_character_))
  if (is.numeric(v) && !is.object(v))
    return(list(values = as.double(v), sql = "REAL", dtype = "float64", categories = NA_character_))
  if (is.atomic(v))      # character, and dates and times as their text
    return(list(values = as.character(v), sql = "TEXT", dtype = "string", categories = NA_character_))
  NULL
}

.write_column_meta <- function(con, table, meta) {
  if (!length(meta)) return(invisible())
  DBI::dbExecute(con,
    "INSERT OR REPLACE INTO _column_meta (table_name, column_name, dtype, categories) VALUES (?, ?, ?, ?)",
    params = list(rep(table, length(meta)),
                  vapply(meta, `[[`, "", 1L), vapply(meta, `[[`, "", 2L),
                  vapply(meta, function(m) as.character(m[[3L]]), "")))
}

## One feature table. Ids must be unique (the id column is UNIQUE in Python's
## schema): repeats get -1, -2 suffixes, as from_anndata gives them, and the
## names as given stay in `symbol`.
.write_feature_table <- function(con, modality, ids, symbols = NULL) {
  tbl <- .modality_table(modality)
  ids <- as.character(ids)
  if (is.null(symbols)) symbols <- ids
  if (anyDuplicated(ids)) {
    warning("cytome: ", sum(duplicated(ids)), " repeated ", modality, " feature ids were made ",
            "unique with -1, -2 suffixes; the names as given are kept as symbols", call. = FALSE)
    ids <- make.unique(ids, sep = "-")
  }
  n <- length(ids)
  if (tbl %in% c("genes", "GA_genes")) {
    DBI::dbExecute(con, sprintf(
      "INSERT INTO %s (gene_idx, gene_id, symbol) VALUES (?, ?, ?)", tbl),
      params = list(seq_len(n) - 1L, ids, as.character(symbols)))
    .write_column_meta(con, tbl, list(list("gene_id", "string", NA_character_),
                                      list("symbol", "string", NA_character_)))
  } else {
    coords <- .parse_peak_ids(ids)
    id_col <- if (tbl == "peaks") "peak" else "tile"
    DBI::dbExecute(con, sprintf(
      "INSERT INTO %s (%s_idx, %s_id, chr, start, end_) VALUES (?, ?, ?, ?, ?)", tbl, id_col, id_col),
      params = list(seq_len(n) - 1L, ids, coords$chr, coords$start, coords$end_))
  }
  invisible(tbl)
}

## One features x cells matrix, stored cells x features as row chunks of
## `chunk_size` cells. A dgCMatrix's columns are the cells, so its column
## pointers and row indices are already the CSR arrays of the stored matrix:
## each chunk is a slice of them, with no transpose.
.cytome_write_matrix <- function(con, name, M, row_entity, col_entity,
                                 dtype = "float32", chunk_size = 2000L,
                                 compression = "zstd") {
  M <- .as_dgC(M, name)
  x <- M@x
  whole <- all(x == round(x), na.rm = TRUE)
  if (dtype == "int32" && (!whole || anyNA(x) || any(abs(x) > .Machine$integer.max)))
    stop("cytome: ", name, " holds values that are not whole numbers in the int32 range; ",
         "write it with dtype = \"float32\" or \"float64\"")
  if (dtype == "float32" && any(is.finite(x) & abs(x) > 3.4028234663852886e38))
    stop("cytome: ", name, " holds values past the float32 range; write it with dtype = \"float64\"")
  n_cells <- ncol(M); p <- M@p; ci <- 0L
  stmt <- paste("INSERT INTO matrix_chunks (matrix_name, chunk_idx, row_start, row_end, n_nonzero,",
                "data_blob, indices_blob, indptr_blob, dtype, compression) VALUES (?,?,?,?,?,?,?,?,?,?)")
  for (lo in .chunk_starts(n_cells, chunk_size)) {
    hi <- min(lo + chunk_size, n_cells)
    a <- p[lo + 1L]; b <- p[hi + 1L]
    sel <- if (b > a) (a + 1L):b else integer(0)
    DBI::dbExecute(con, stmt, params = list(
      name, ci, lo, hi, b - a,
      list(cytome_compress(.value_raw(x[sel], dtype), compression)),
      list(cytome_compress(.int32_raw(M@i[sel]), compression)),
      list(cytome_compress(.int32_raw(p[(lo + 1L):(hi + 1L)] - a), compression)),
      dtype, compression))
    ci <- ci + 1L
  }
  DBI::dbExecute(con, paste(
    "INSERT INTO matrix_meta (matrix_name, n_rows, n_cols, n_nonzero, dtype, row_entity, col_entity,",
    "chunk_size, n_chunks, has_csc, created_at, is_integer) VALUES (?,?,?,?,?,?,?,?,?,0,?,?)"),
    params = list(name, n_cells, nrow(M), length(x), dtype, row_entity, col_entity,
                  as.integer(chunk_size), ci, .now(), as.integer(whole && !anyNA(x))))
}

## Embeddings: dense, row-major float32 chunks of `chunk_size` cells. Names
## become X_<name> unless they already start with X_.
.write_embeddings <- function(con, emb_list, chunk_size, compression) {
  keys <- vapply(names(emb_list), function(nm) if (startsWith(nm, "X_")) nm else paste0("X_", tolower(nm)), "")
  if (anyDuplicated(tolower(keys))) keys <- make.unique(keys, sep = "_")
  for (k in seq_along(emb_list)) {
    E <- emb_list[[k]]; key <- keys[[k]]
    n <- nrow(E); ci <- 0L
    for (lo in .chunk_starts(n, chunk_size)) {
      hi <- min(lo + chunk_size, n)
      blk <- E[(lo + 1L):hi, , drop = FALSE]
      DBI::dbExecute(con, paste(
        "INSERT INTO dense_chunks (array_name, chunk_idx, row_start, row_end, n_cols, data_blob, dtype, compression)",
        "VALUES (?,?,?,?,?,?,?,?)"),
        params = list(key, ci, lo, hi, ncol(E),
                      list(cytome_compress(writeBin(as.double(t(blk)), raw(), size = 4L, endian = "little"),
                                           compression)),
                      "float32", compression))
      ci <- ci + 1L
    }
    DBI::dbExecute(con, paste(
      "INSERT INTO embedding_meta (array_name, n_rows, n_cols, dtype, entity, chunk_size, n_chunks, created_at)",
      "VALUES (?,?,?,?,?,?,?,?)"),
      params = list(key, n, ncol(E), "float32", "cells", as.integer(chunk_size), ci, .now()))
  }
}

## Graphs as Python's GraphStore writes them (cytome 0.3.6): rows in chunks of
## about 2^20 entries, values float64, column indices int32 and a chunk-local
## int64 index pointer, each lz4-compressed; graph_meta records the shape. Each
## is also named in `_anndata_obsp_map`, so Python's to_anndata() puts it in obsp.
.GRAPH_CHUNK_NNZ <- 2^20

.write_graphs <- function(con, graph_list) {
  for (nm in names(graph_list)) {
    G <- graph_list[[nm]]
    R <- Matrix::t(G)                         # columns of t(G) are the rows of G: CSR
    p <- R@p; n <- ncol(R); r0 <- 0L; ci <- 0L
    repeat {
      r1 <- findInterval(p[r0 + 1L] + .GRAPH_CHUNK_NNZ, p) - 1L
      r1 <- min(max(r1, r0 + 1L), n)
      if (n == 0L) r1 <- 0L
      a <- p[r0 + 1L]; b <- p[r1 + 1L]
      sel <- if (b > a) (a + 1L):b else integer(0)
      DBI::dbExecute(con, paste(
        "INSERT INTO graph_chunks (graph_name, axis, chunk_idx, row_start, row_end, n_nonzero,",
        "data_blob, indices_blob, indptr_blob, compression) VALUES (?, 'obs', ?, ?, ?, ?, ?, ?, ?, 'lz4')"),
        params = list(nm, ci, r0, r1, b - a,
                      list(cytome_compress(writeBin(as.double(R@x[sel]), raw(), size = 8L, endian = "little"), "lz4")),
                      list(cytome_compress(.int32_raw(R@i[sel]), "lz4")),
                      list(cytome_compress(.int64_raw(p[(r0 + 1L):(r1 + 1L)] - a), "lz4"))))
      ci <- ci + 1L
      r0 <- r1
      if (r0 >= n) break
    }
    DBI::dbExecute(con, paste(
      "INSERT INTO graph_meta (graph_name, axis, entity_table, n_rows, n_cols, n_nonzero, dtype, n_chunks)",
      "VALUES (?, 'obs', 'cells', ?, ?, ?, 'float64', ?)"),
      params = list(nm, nrow(G), ncol(G), length(R@x), ci))
  }
  if (length(graph_list))
    .write_metadata(con, "_anndata_obsp_map",
                    structure(as.list(names(graph_list)), names = names(graph_list)), "dict")
}

## The 0-based first rows of chunks of `size` over `n` rows (none when n is 0).
.chunk_starts <- function(n, size) if (n > 0L) seq.int(0L, n - 1L, by = size) else integer(0)

.write_metadata <- function(con, key, value, value_type) {
  DBI::dbExecute(con, "INSERT OR REPLACE INTO _metadata (key, value, value_type) VALUES (?, ?, ?)",
                 params = list(key, .json(value), value_type))
}

.write_manifest <- function(con, path, n_cells, modalities) {
  ver <- tryCatch(as.character(utils::packageVersion("cytome")), error = function(e) "unknown")
  kv <- list(format_version = .json("1.0.0"), min_reader_version = .json("1.0.0"),
             writer_version = .json(paste("cytome-r", ver)), created_at = .json(.now()),
             dataset_name = .json(sub("\\.cytome$", "", basename(path))), genome = "null",
             modalities = .json(I(modalities)), n_cells = .json(n_cells))
  DBI::dbExecute(con, "INSERT INTO _manifest (key, value) VALUES (?, ?)",
                 params = list(names(kv), unlist(kv, use.names = FALSE)))
}

## Write `spec` to `path`: a temporary file beside it, one transaction, then a
## rename over the destination. spec: barcodes, obs, modalities (each: name, ids,
## symbols, matrices -- a named list of features x cells over all of `ids`, the
## counts first), embeddings, graphs.
.write_store <- function(path, spec, dtype, chunk_size, compression) {
  dtype <- .check_dtype(dtype)
  if (!(is.character(compression) && length(compression) == 1L && compression %in% c("zstd", "lz4", "zlib")))
    stop("cytome: compression must be \"zstd\", \"lz4\" or \"zlib\"")
  chunk_size <- as.integer(chunk_size)
  if (length(chunk_size) != 1L || is.na(chunk_size) || chunk_size < 1L)
    stop("cytome: chunk_size must be a positive number of cells")
  path <- normalizePath(path, mustWork = FALSE)
  if (!dir.exists(dirname(path))) stop("cytome: folder ", dirname(path), " does not exist")
  tmp <- tempfile(pattern = paste0(".", basename(path), ".partial-"), tmpdir = dirname(path))
  on.exit(unlink(paste0(tmp, c("", "-journal", "-wal", "-shm"))), add = TRUE)

  con <- DBI::dbConnect(RSQLite::SQLite(), tmp, synchronous = NULL)
  tryCatch({
    ## A new file in one transaction: no journal is needed, since a failure
    ## deletes the file, and the commit is synced before the rename.
    DBI::dbExecute(con, "PRAGMA journal_mode = OFF")
    DBI::dbExecute(con, "PRAGMA synchronous = FULL")
    DBI::dbBegin(con)
    .cytome_write_schema(con)
    .write_cells(con, spec$barcodes, spec$obs)
    for (m in spec$modalities) {
      tbl <- .write_feature_table(con, m$name, m$ids, m$symbols)
      for (lyr in names(m$matrices))
        .cytome_write_matrix(con, paste0(m$name, "_", lyr), m$matrices[[lyr]], "cells", tbl,
                             dtype, chunk_size, compression)
    }
    if (length(spec$embeddings)) .write_embeddings(con, spec$embeddings, chunk_size, compression)
    if (length(spec$graphs)) .write_graphs(con, spec$graphs)
    .write_manifest(con, path, length(spec$barcodes),
                    vapply(spec$modalities, `[[`, "", "name"))
    DBI::dbCommit(con)
  }, finally = DBI::dbDisconnect(con))
  ## The old store's -wal would be replayed over the new file: remove it first.
  unlink(paste0(path, c("-wal", "-shm", "-journal")))
  if (!file.rename(tmp, path))
    stop("cytome: could not move the new store to ", path, " (is it open elsewhere?)")
  invisible(path)
}

## How a layer left out of the store is recomputed after reading. `sct`: the
## layer is from a Seurat SCTransform assay, whose model is not stored.
.recompute_hint <- function(layer, sct = FALSE) {
  if (sct) return("with SCTransform(); the SCT model is not stored")
  switch(layer, scale.data = "with ScaleData()", data = "with NormalizeData()", "from the counts")
}

## One modality's entry of a write spec: its counts and layers aligned to its
## features and to `cells`, or NULL when its counts cover only some of the cells
## (refused for the main modality, which the store is written around). A layer
## over only some of the features or cells is left out. Each is named in a
## warning that says how to get it back.
.spec_modality <- function(name, ids, symbols, counts, layers, cells, what, sct = FALSE,
                           main = FALSE) {
  mod <- .canonical_modality(name)
  M <- .align_features(counts, ids, what)
  if (is.null(M)) stop("cytome: ", what, " covers only ", nrow(counts), " of ", length(ids), " features")
  M <- .align_cells(M, cells, what)
  if (is.null(M)) {
    if (main)
      stop("cytome: ", what, " covers only ", ncol(counts), " of the object's ", length(cells),
           " cells. Subset the object to those cells and write that, or choose another main ",
           "assay.", call. = FALSE)
    .some_cells_warning(what, ncol(counts), length(cells))
    return(NULL)
  }
  mats <- list(counts = M)
  for (lyr in names(layers)) {
    w <- paste0(what, " layer '", lyr, "'")
    L <- .align_features(layers[[lyr]], ids, w)
    if (is.null(L)) {
      warning("cytome: ", mod, " layer '", lyr, "' is not written: it covers only ",
              nrow(layers[[lyr]]), " of ", length(ids), " features. Recompute it after reading, ",
              .recompute_hint(lyr, sct), ".", call. = FALSE)
      next
    }
    L <- .align_cells(L, cells, w)
    if (is.null(L)) {
      .some_cells_warning(paste0(mod, " layer '", lyr, "'"), ncol(layers[[lyr]]), length(cells))
      next
    }
    mats[[lyr]] <- L
  }
  list(name = mod, ids = ids, symbols = symbols, matrices = mats)
}

#' Write an object to a `.cytome` file
#'
#' Generic. The method is chosen from the class of `x`, so the same call works
#' whether you are holding a Seurat object or a SingleCellExperiment -- the
#' thing you would otherwise have to look up is the thing R already knows.
#' Third-party packages can add methods for their own classes.
#'
#' The store is built in a temporary file beside `path`, in one transaction,
#' and moves to `path` only once it is complete: a write that fails leaves an
#' existing file as it was. Values are placed by feature and cell names, so a
#' matrix whose genes or cells come in another order is stored under its own.
#' Column types travel: factors with their levels, logicals and integers.
#'
#' cytome never stores values that were not measured. A layer over only some
#' of the genes (a Seurat `scale.data` over the variable genes) is not written,
#' with a warning that names it: stored over all the genes, the zeros of the
#' genes it lacks would read as values, and in `scale.data` 0 is the mean.
#' Recompute it after reading, with `ScaleData()`. An assay, layer or graph
#' over only some of the cells (an assay measured in some batches only) is not
#' written either, with a warning that names it: to keep it, subset the object
#' to those cells and write that as a separate store. An embedding keeps `NaN`
#' for the cells it lacks.
#'
#' @param x A `Seurat` or `SingleCellExperiment` object.
#' @param path Output `.cytome` path.
#' @param ... Passed to the method.
#' @return `path`, invisibly. The file is readable by Python `cytome.open(path)`.
#' @seealso [read_cytome()]
#' @examples
#' \dontrun{
#' write_cytome(seurat_obj, "out.cytome")
#' write_cytome(sce, "out.cytome")
#' }
#' @export
write_cytome <- function(x, path, ...) UseMethod("write_cytome")

#' @rdname write_cytome
#' @export
write_cytome.default <- function(x, path, ...)
  stop("cytome: no write_cytome() method for class '", paste(class(x), collapse = "/"),
       "'. Supported: Seurat, SingleCellExperiment.")

#' Write a SingleCellExperiment to a `.cytome` file
#'
#' @param x A [SingleCellExperiment::SingleCellExperiment].
#' @param path Output `.cytome` path.
#' @param assay Assay to store as the main counts matrix (default `"counts"`).
#' @param main_modality Modality name for the main assay (default `"RNA"`; also
#'   `"ATAC"`, `"GA"` or `"tiles"`).
#' @param alt_modalities `altExp` names to write as extra modalities. `NULL`
#'   (default) writes every one cytome has a feature table for (RNA, GA, ATAC,
#'   tiles) and names the others (ERCC, ADT, ...) in a warning; naming one that
#'   cannot be written is an error. `character(0)` writes none.
#' @param layers Additional assays of the main experiment to carry.
#'   `FALSE` (default) writes none, `TRUE` writes all, or name them. Off by
#'   default because normalized values are recomputable from the counts and
#'   cost more space than the counts do; turn it on when the point is to hand
#'   somebody your object rather than to archive the measurements.
#' @param graphs Write `colPairs` as the Python package stores graphs
#'   (default `TRUE`). A neighbour graph is the expensive thing to recompute,
#'   and carrying it means the recipient reproduces the same clusters.
#' @param dtype `"float32"` (default), `"float64"`, or `"int32"` (whole numbers
#'   only; an error otherwise).
#' @param chunk_size Rows (cells) per storage chunk.
#' @param compression `"zstd"` (default), `"lz4"`, or `"zlib"`.
#' @param embeddings Store `reducedDims` (default `TRUE`).
#' @param ... Unused.
#' @return `path`, invisibly.
#' @export
write_cytome.SingleCellExperiment <- function(x, path, assay = "counts",
                                              main_modality = "RNA",
                                              alt_modalities = NULL,
                                              layers = FALSE, graphs = TRUE,
                                              dtype = "float32", chunk_size = 2000L,
                                              compression = "zstd", embeddings = TRUE, ...) {
  stopifnot(requireNamespace("SingleCellExperiment", quietly = TRUE))
  .modality_table(main_modality)
  cells <- colnames(x); if (is.null(cells)) cells <- paste0("cell", seq_len(ncol(x)))
  fid <- rownames(x); if (is.null(fid)) fid <- paste0("gene", seq_len(nrow(x)))
  if (anyDuplicated(cells)) stop("cytome: the cell names (colnames) repeat; cells must be unique")
  sym <- .detect_feature_symbol(as.data.frame(SummarizedExperiment::rowData(x)), nrow(x))

  ## Additional assays of the MAIN experiment -- logcounts and friends. Off by
  ## default: for an archive they are recomputable and cost more space than the
  ## counts (a scaled float matrix compresses far worse than integers). For a
  ## handoff they are the difference between sending your object and sending
  ## counts plus an instruction to guess the parameters.
  ## Assays are taken with the object's own names: withDimnames = FALSE returns
  ## the names they had when they were set, which renaming the genes or cells
  ## (Ensembl ids to symbols, say) leaves behind.
  lyr <- .resolve_layers(layers, setdiff(SummarizedExperiment::assayNames(x), assay))
  lmats <- lapply(lyr, function(a) SummarizedExperiment::assay(x, a))
  names(lmats) <- lyr
  main <- .spec_modality(main_modality, fid, sym, SummarizedExperiment::assay(x, assay),
                         lmats, cells, paste0("assay '", assay, "'"), main = TRUE)
  mods <- list(main)

  ## altExps as additional modalities. Previously they were silently dropped,
  ## so a multimodal object round-tripped as RNA only and nothing said so.
  for (m in .resolve_alt(alt_modalities, SingleCellExperiment::altExpNames(x), "altExp")) {
    ae <- SingleCellExperiment::altExp(x, m)
    if (.canonical_modality(m) %in% vapply(mods, `[[`, "", "name"))
      stop("cytome: two modalities would both be ", .canonical_modality(m))
    aid <- rownames(ae); if (is.null(aid)) aid <- paste0(m, seq_len(nrow(ae)))
    spec <- .spec_modality(
      m, aid, .detect_feature_symbol(as.data.frame(SummarizedExperiment::rowData(ae)), nrow(ae)),
      SummarizedExperiment::assay(ae, 1L), list(), cells, paste0("altExp '", m, "'"))
    if (!is.null(spec)) mods[[length(mods) + 1L]] <- spec
  }

  emb <- if (embeddings) as.list(SingleCellExperiment::reducedDims(x)) else list()
  emb <- lapply(names(emb), function(e) .align_embedding(emb[[e]], cells, e))
  names(emb) <- if (embeddings) SingleCellExperiment::reducedDimNames(x) else NULL
  gr <- if (graphs) SingleCellExperiment::colPairNames(x) else character(0)
  gl <- lapply(gr, function(g) .align_graph(.colpair_matrix(x, g), cells, g))
  names(gl) <- gr
  gl <- gl[!vapply(gl, is.null, TRUE)]

  .write_store(path, list(barcodes = cells,
                          obs = as.data.frame(SingleCellExperiment::colData(x), optional = TRUE),
                          modalities = mods, embeddings = emb, graphs = gl),
               dtype, chunk_size, compression)
}

## A colPair as a sparse matrix, as colPair(x, g, asSparse = TRUE) makes it
## (the first mcols column as values, the last value of a repeated pair), with
## the repeats found from one numeric key per pair instead of a duplicated()
## over a two-column matrix, which took 1.4 s of a 4.7 s write of 300,000 edges.
.colpair_matrix <- function(x, g) {
  h <- SingleCellExperiment::colPair(x, g)
  n <- S4Vectors::nnode(h)
  i <- S4Vectors::queryHits(h); j <- S4Vectors::subjectHits(h)
  m <- S4Vectors::mcols(h)
  v <- if (ncol(m)) m[[1L]] else rep(1, length(i))
  if (!is.numeric(v) && !is.logical(v))
    stop("cytome: colPair '", g, "' holds values of type '", typeof(v), "'")
  keep <- !duplicated((as.numeric(i) - 1) * n + j, fromLast = TRUE)
  Matrix::sparseMatrix(i = i[keep], j = j[keep], x = as.numeric(v[keep]), dims = c(n, n))
}

## The base name of a Seurat layer: "counts" for the split layers counts.A,
## counts.B of a Seurat v5 object.
.layer_base <- function(l) {
  for (b in c("scale.data", "counts", "data"))
    if (l == b || startsWith(l, paste0(b, "."))) return(b)
  l
}

## Seurat layer names other than the one already written as counts.
.seurat_layers <- function(obj, assay, written) {
  nms <- tryCatch(SeuratObject::Layers(obj[[assay]]), error = function(e) NULL)
  if (is.null(nms))
    nms <- intersect(c("counts", "data", "scale.data"),
                     methods::slotNames(obj@assays[[assay]]))
  setdiff(unique(vapply(nms, .layer_base, "")), written)
}

## Read one assay's matrix out of a Seurat object without requiring the full
## Seurat package. A Seurat v5 assay split by sample (counts.A, counts.B, ...)
## is joined first, so every cell is written. GetAssayData() dispatches
## through Seurat's S4 method, which calls UpdateSlots() and therefore needs
## Seurat itself -- deposited v4 objects then fail on a machine that only has
## SeuratObject. Try the supported call, fall back to the slot.
.seurat_layer <- function(obj, assay, layer) {
  a <- obj[[assay]]
  if (methods::is(a, "Assay5")) {
    found <- suppressWarnings(SeuratObject::Layers(a, search = layer))
    if (!length(found)) stop(sprintf("assay '%s' has no layer '%s'", assay, layer))
    if (!identical(found, layer)) a <- SeuratObject::JoinLayers(a, layers = layer)
    return(SeuratObject::LayerData(a, layer = layer))
  }
  tryCatch(
    SeuratObject::GetAssayData(obj, assay = assay, layer = layer),
    error = function(e) {
      as_ <- obj@assays[[assay]]
      slot_name <- if (layer %in% methods::slotNames(as_)) layer else "counts"
      m <- methods::slot(as_, slot_name)
      if (is.null(m) || !length(m))
        stop(sprintf("assay '%s' has no data in slot '%s'", assay, slot_name))
      m
    })
}

#' Write a Seurat object to a `.cytome` file
#'
#' @param x A `Seurat` object.
#' @param path Output `.cytome` path.
#' @param assay Assay for the main modality (default the object's default assay).
#'   It must be one cytome has a feature table for (RNA, GA, ATAC, tiles).
#' @param layer Layer/slot to store (default `"counts"`). The split layers of a
#'   Seurat v5 assay (`counts.A`, `counts.B`, ...) are joined.
#' @param alt_assays Other assays to write as extra modalities. `NULL` (default)
#'   writes every one cytome has a feature table for and names the others (SCT,
#'   ADT, ...) in a warning; naming one that cannot be written is an error.
#'   `character(0)` writes none. An assay over only some of the object's cells
#'   is not written, with a warning; the main assay must cover them all.
#' @param layers Other layers of the main assay (`data`, `scale.data`) to
#'   carry. `FALSE` (default), `TRUE`, or name them. See the
#'   `SingleCellExperiment` method for why this is off by default. A layer over
#'   only some of the features, such as a `scale.data` over the variable genes,
#'   is not written: a warning names it, and it is recomputed after reading, with
#'   `ScaleData()` (`SCTransform()` for an SCT assay). Nor is a layer over only
#'   some of the cells.
#' @param graphs Write the `graphs` slot as the Python package stores graphs
#'   (default `TRUE`). A graph over only some of the cells is not written, with
#'   a warning: a reader would take a cell it lacks for a cell without
#'   neighbours.
#' @param dtype,chunk_size,compression,embeddings As for the
#'   `SingleCellExperiment` method.
#' @param ... Unused.
#' @return `path`, invisibly.
#' @export
write_cytome.Seurat <- function(x, path, assay = NULL, layer = "counts",
                                alt_assays = NULL, layers = FALSE, graphs = TRUE,
                                dtype = "float32",
                                chunk_size = 2000L, compression = "zstd",
                                embeddings = TRUE, ...) {
  stopifnot(requireNamespace("SeuratObject", quietly = TRUE))
  if (is.null(assay)) assay <- SeuratObject::DefaultAssay(x)
  if (is.na(.canonical_modality(assay)))
    stop("cytome: cannot write assay '", assay, "' as the main modality: cytome has no feature ",
         "table for it (writable: ", paste(names(.MODALITY_TABLE), collapse = ", "), "). ",
         "Choose one with assay = , e.g. assay = \"RNA\".")
  cells <- colnames(x)
  A <- .seurat_layer(x, assay, layer)                      # features x cells
  fid <- rownames(A)
  rmeta <- tryCatch(as.data.frame(x[[assay]][[]]), error = function(e) NULL)
  if (!is.null(rmeta) && !is.null(rownames(rmeta)) && !identical(rownames(rmeta), fid))
    rmeta <- rmeta[match(fid, rownames(rmeta)), , drop = FALSE]

  ## Other layers of the MAIN assay (data, scale.data). See the SCE method for
  ## why this is off by default.
  lmats <- list()
  for (lyr in .resolve_layers(layers, .seurat_layers(x, assay, layer))) {
    lm <- tryCatch(.seurat_layer(x, assay, lyr), error = function(e) NULL)
    if (is.null(lm) || !length(lm)) next
    lmats[[lyr]] <- lm
  }
  mods <- list(.spec_modality(assay, fid, .detect_feature_symbol(rmeta, length(fid)), A, lmats,
                              cells, paste0("assay '", assay, "'"),
                              sct = inherits(x[[assay]], "SCTAssay"), main = TRUE))

  ## Other assays. They used to be dropped silently: a Seurat object with RNA
  ## and ATAC round-tripped as RNA only, and the file gave no sign that half of
  ## it was missing.
  for (m in .resolve_alt(alt_assays, setdiff(names(x@assays), assay), "assay")) {
    am <- .seurat_layer(x, m, layer)
    if (.canonical_modality(m) %in% vapply(mods, `[[`, "", "name"))
      stop("cytome: two modalities would both be ", .canonical_modality(m))
    am_meta <- tryCatch(as.data.frame(x[[m]][[]]), error = function(e) NULL)
    spec <- .spec_modality(
      m, rownames(am), .detect_feature_symbol(am_meta, nrow(am)), am, list(), cells,
      paste0("assay '", m, "'"))
    if (!is.null(spec)) mods[[length(mods) + 1L]] <- spec
  }

  reds <- if (embeddings) SeuratObject::Reductions(x) else character(0)
  emb <- lapply(reds, function(r) .align_embedding(SeuratObject::Embeddings(x, r), cells, r))
  names(emb) <- reds
  grs <- if (graphs) tryCatch(names(x@graphs), error = function(e) NULL) else NULL
  gl <- lapply(grs, function(g) .align_graph(x@graphs[[g]], cells, g))
  names(gl) <- grs
  gl <- gl[!vapply(gl, is.null, TRUE)]

  .write_store(path, list(barcodes = cells, obs = x[[]], modalities = mods,
                          embeddings = emb, graphs = gl),
               dtype, chunk_size, compression)
}
