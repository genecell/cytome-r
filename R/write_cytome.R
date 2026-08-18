# Native R writer for the .cytome SQLite format (no Python). Mirrors
# cytome/io/sqlite_engine.py (_create_schema) + cytome/io/chunked_io.py
# (write_sparse_chunked) so a cytome written here is read bit-for-bit by the
# Python `cytome` package (cross-language write conformance).

.cytome_write_schema <- function(con) {
  DBI::dbExecute(con, "PRAGMA journal_mode=WAL;")
  stmts <- c(
    "CREATE TABLE IF NOT EXISTS _manifest (key TEXT PRIMARY KEY, value TEXT);",
    "CREATE TABLE IF NOT EXISTS _metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL, value_type TEXT NOT NULL);",
    "CREATE TABLE IF NOT EXISTS _column_meta (table_name TEXT NOT NULL, column_name TEXT NOT NULL, dtype TEXT NOT NULL, categories TEXT, PRIMARY KEY (table_name, column_name));",
    "CREATE TABLE IF NOT EXISTS cells (cell_idx INTEGER PRIMARY KEY, barcode TEXT NOT NULL, sample_id TEXT, n_fragments INTEGER DEFAULT 0);",
    "CREATE TABLE IF NOT EXISTS genes (gene_idx INTEGER PRIMARY KEY, gene_id TEXT NOT NULL UNIQUE, symbol TEXT, chr TEXT, start INTEGER, end_ INTEGER, biotype TEXT);",
    ## ATAC features. Python's schema declares chr/start/end_ NOT NULL, so a
    ## peak_id that does not parse as chr:start-end would be rejected there --
    ## .write_feature_table() parses it and errors here instead, where the user
    ## can still see which id was the problem.
    "CREATE TABLE IF NOT EXISTS peaks (peak_idx INTEGER PRIMARY KEY, peak_id TEXT NOT NULL, chr TEXT NOT NULL, start INTEGER NOT NULL, end_ INTEGER NOT NULL, annotation TEXT, nearest_gene TEXT, distance_to_tss INTEGER);",
    paste0("CREATE TABLE IF NOT EXISTS matrix_chunks (id INTEGER PRIMARY KEY, matrix_name TEXT NOT NULL, ",
           "chunk_idx INTEGER NOT NULL, row_start INTEGER NOT NULL, row_end INTEGER NOT NULL, n_nonzero INTEGER NOT NULL, ",
           "data_blob BLOB NOT NULL, indices_blob BLOB NOT NULL, indptr_blob BLOB NOT NULL, dtype TEXT NOT NULL, ",
           "compression TEXT DEFAULT 'zstd', UNIQUE(matrix_name, chunk_idx));"),
    "CREATE INDEX IF NOT EXISTS idx_matrix_chunks ON matrix_chunks(matrix_name, chunk_idx);",
    paste0("CREATE TABLE IF NOT EXISTS matrix_meta (matrix_name TEXT PRIMARY KEY, n_rows INTEGER NOT NULL, ",
           "n_cols INTEGER NOT NULL, n_nonzero INTEGER NOT NULL, dtype TEXT NOT NULL, row_entity TEXT NOT NULL, ",
           "col_entity TEXT NOT NULL, chunk_size INTEGER NOT NULL, n_chunks INTEGER NOT NULL, has_csc INTEGER DEFAULT 0, ",
           "csc_chunk_size INTEGER, csc_n_chunks INTEGER, created_at TEXT, provenance_id INTEGER);"),
    paste0("CREATE TABLE IF NOT EXISTS dense_chunks (id INTEGER PRIMARY KEY, array_name TEXT NOT NULL, ",
           "chunk_idx INTEGER NOT NULL, row_start INTEGER NOT NULL, row_end INTEGER NOT NULL, n_cols INTEGER NOT NULL, ",
           "data_blob BLOB NOT NULL, dtype TEXT NOT NULL, compression TEXT DEFAULT 'zstd', UNIQUE(array_name, chunk_idx));"),
    paste0("CREATE TABLE IF NOT EXISTS embedding_meta (array_name TEXT PRIMARY KEY, n_rows INTEGER NOT NULL, ",
           "n_cols INTEGER NOT NULL, dtype TEXT NOT NULL, entity TEXT NOT NULL, chunk_size INTEGER NOT NULL, ",
           "n_chunks INTEGER NOT NULL, created_at TEXT, provenance_id INTEGER);"))
  for (s in stmts) DBI::dbExecute(con, s)
}

.float_raw <- function(x, dtype) writeBin(as.double(x), raw(),
                                          size = if (dtype == "float32") 4L else 8L, endian = "little")
.int_raw   <- function(x) writeBin(as.integer(x), raw(), size = 4L, endian = "little")

# Write one cells x features matrix into chunked, compressed CSR storage.
.cytome_write_matrix <- function(con, name, M, row_entity, col_entity,
                                 dtype = "float32", chunk_size = 2000L,
                                 compression = "zstd") {
  # M: cells x features; store row-major (CSR). Convert to RsparseMatrix (CSR).
  Mr <- methods::as(M, "RsparseMatrix")
  n_rows <- nrow(Mr); n_cols <- ncol(Mr); total_nnz <- 0L; ci <- 0L
  for (lo in seq.int(1L, max(n_rows, 1L), by = chunk_size)) {
    hi <- min(lo + chunk_size - 1L, n_rows)
    sub <- if (n_rows == 0L) Mr else methods::as(Mr[lo:hi, , drop = FALSE], "RsparseMatrix")
    p <- sub@p; j <- sub@j; x <- sub@x                    # CSR: p (nrow+1), j (0-based), x
    DBI::dbExecute(con,
      "INSERT INTO matrix_chunks(matrix_name, chunk_idx, row_start, row_end, n_nonzero, data_blob, indices_blob, indptr_blob, dtype, compression) VALUES (?,?,?,?,?,?,?,?,?,?)",
      params = list(name, ci, lo - 1L, hi, length(x),
                    list(cytome_compress(.float_raw(x, dtype), compression)),
                    list(cytome_compress(.int_raw(j), compression)),
                    list(cytome_compress(.int_raw(p), compression)),
                    dtype, compression))
    total_nnz <- total_nnz + length(x); ci <- ci + 1L
  }
  DBI::dbExecute(con,
    "INSERT INTO matrix_meta(matrix_name, n_rows, n_cols, n_nonzero, dtype, row_entity, col_entity, chunk_size, n_chunks, has_csc, created_at) VALUES (?,?,?,?,?,?,?,?,?,0,?)",
    params = list(name, n_rows, n_cols, total_nnz, dtype, row_entity, col_entity,
                  chunk_size, ci, format(Sys.time(), "%Y-%m-%dT%H:%M:%S")))
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

## Modality -> feature table, mirroring MODALITY_REGISTRY in
## cytome/utils/modality.py. Kept deliberately short: R writes the two
## modalities it can construct feature rows for. Anything else is an error
## naming what is supported, rather than a cytome Python cannot open.
.MODALITY_TABLE <- c(RNA = "genes", ATAC = "peaks")

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

## Seurat layer names other than the one already written as counts.
.seurat_layers <- function(obj, assay, written) {
  nms <- tryCatch(SeuratObject::Layers(obj[[assay]]), error = function(e) NULL)
  if (is.null(nms))
    nms <- intersect(c("counts", "data", "scale.data"),
                     methods::slotNames(obj@assays[[assay]]))
  setdiff(nms, written)
}

.modality_table <- function(modality) {
  key <- toupper(modality)
  if (!key %in% names(.MODALITY_TABLE))
    stop("cytome: cannot write modality '", modality, "' from R. ",
         "Supported: ", paste(names(.MODALITY_TABLE), collapse = ", "),
         ". Build the cytome on the Python side for other modalities.")
  .MODALITY_TABLE[[key]]
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
             start = as.integer(vapply(m, `[`, character(1), 3)),
             end_ = as.integer(vapply(m, `[`, character(1), 4)),
             stringsAsFactors = FALSE)
}

## One feature table for one modality.
.write_feature_table <- function(con, modality, feature_ids, feature_symbol = NULL) {
  tbl <- .modality_table(modality)
  ids <- as.character(feature_ids)
  if (tbl == "genes") {
    DBI::dbWriteTable(con, "genes",
      data.frame(gene_idx = seq_along(ids) - 1L, gene_id = ids,
                 symbol = if (is.null(feature_symbol)) ids else as.character(feature_symbol),
                 chr = NA_character_, start = NA_integer_, end_ = NA_integer_,
                 biotype = NA_character_, stringsAsFactors = FALSE),
      append = TRUE, row.names = FALSE)
  } else {
    coords <- .parse_peak_ids(ids)
    DBI::dbWriteTable(con, "peaks",
      data.frame(peak_idx = seq_along(ids) - 1L, peak_id = ids,
                 chr = coords$chr, start = coords$start, end_ = coords$end_,
                 annotation = NA_character_, nearest_gene = NA_character_,
                 distance_to_tss = NA_integer_, stringsAsFactors = FALSE),
      append = TRUE, row.names = FALSE)
  }
  invisible(tbl)
}

.write_entities <- function(con, barcodes, feature_ids, feature_symbol = NULL) {
  DBI::dbWriteTable(con, "cells",
    data.frame(cell_idx = seq_along(barcodes) - 1L, barcode = as.character(barcodes),
               sample_id = NA_character_, n_fragments = 0L, stringsAsFactors = FALSE),
    append = TRUE, row.names = FALSE)
  .write_feature_table(con, "RNA", feature_ids, feature_symbol)
}

.write_obs_cols <- function(con, obs) {
  # add each obs column to the cells table (TEXT/REAL), so Python cytome_obs() sees them
  if (is.null(obs) || ncol(obs) == 0L) return(invisible())
  for (cn in colnames(obs)) {
    v <- obs[[cn]]; sqlcol <- gsub("[^A-Za-z0-9_]", "_", cn)
    if (sqlcol %in% c("cell_idx", "barcode", "sample_id", "n_fragments")) sqlcol <- paste0(sqlcol, "_obs")
    sqltype <- if (is.numeric(v)) "REAL" else "TEXT"
    DBI::dbExecute(con, sprintf("ALTER TABLE cells ADD COLUMN \"%s\" %s;", sqlcol, sqltype))
    vv <- if (is.numeric(v)) as.numeric(v) else as.character(v)
    DBI::dbExecute(con, sprintf("UPDATE cells SET \"%s\" = ? WHERE cell_idx = ?;", sqlcol),
                   params = list(vv, seq_along(vv) - 1L))
  }
}

## Graphs go into the SAME table Python writes, so there is no second
## convention to drift from:
##   graph_edges(graph_name, axis, entity_table, row_idx, col_idx, value)
.write_graphs <- function(con, graph_list, axis = "obs", entity_table = "cells") {
  if (!length(graph_list)) return(invisible())
  DBI::dbExecute(con, paste0(
    "CREATE TABLE IF NOT EXISTS graph_edges (graph_name TEXT NOT NULL, ",
    "axis TEXT NOT NULL DEFAULT 'obs', entity_table TEXT NOT NULL DEFAULT 'cells', ",
    "row_idx INTEGER NOT NULL, col_idx INTEGER NOT NULL, value REAL NOT NULL);"))
  DBI::dbExecute(con,
    "CREATE INDEX IF NOT EXISTS idx_graph_edges ON graph_edges(graph_name, row_idx);")
  for (nm in names(graph_list)) {
    G <- methods::as(graph_list[[nm]], "TsparseMatrix")     # i, j, x triplets
    if (!length(G@x)) next
    n <- length(G@x)
    ## Every bound parameter must have the same length: RSQLite recycles
    ## nothing, and a scalar next to a vector is "Parameter 4 does not have
    ## length 1", which reads as a type error rather than a length one.
    DBI::dbExecute(con,
      "INSERT INTO graph_edges(graph_name, axis, entity_table, row_idx, col_idx, value) VALUES (?,?,?,?,?,?)",
      params = list(rep(nm, n), rep(axis, n), rep(entity_table, n),
                    as.integer(G@i), as.integer(G@j), as.numeric(G@x)))
  }
}

.write_embeddings <- function(con, emb_list) {
  # emb_list: named list of (cells x k) matrices → dense_chunks + embedding_meta
  for (nm in names(emb_list)) {
    E <- as.matrix(emb_list[[nm]]); key <- if (startsWith(nm, "X_")) nm else paste0("X_", tolower(nm))
    storage.mode(E) <- "double"
    raw <- cytome_compress(writeBin(as.double(t(E)), raw(), size = 4L, endian = "little"), "zstd")  # row-major float32
    DBI::dbExecute(con,
      "INSERT INTO dense_chunks(array_name, chunk_idx, row_start, row_end, n_cols, data_blob, dtype, compression) VALUES (?,0,0,?,?,?,?,?)",
      params = list(key, nrow(E), ncol(E), list(raw), "float32", "zstd"))
    DBI::dbExecute(con,
      "INSERT INTO embedding_meta(array_name, n_rows, n_cols, dtype, entity, chunk_size, n_chunks, created_at) VALUES (?,?,?,?,?,?,1,?)",
      params = list(key, nrow(E), ncol(E), "float32", "cells", nrow(E),
                    format(Sys.time(), "%Y-%m-%dT%H:%M:%S")))
  }
}

#' Write an object to a `.cytome` file
#'
#' Generic. The method is chosen from the class of `x`, so the same call works
#' whether you are holding a Seurat object or a SingleCellExperiment -- the
#' thing you would otherwise have to look up is the thing R already knows.
#' Third-party packages can add methods for their own classes.
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

#' @export
write_cytome.default <- function(x, path, ...)
  stop("cytome: no write_cytome() method for class '", paste(class(x), collapse = "/"),
       "'. Supported: Seurat, SingleCellExperiment.")

#' Write a SingleCellExperiment to a `.cytome` file
#'
#' @param x A [SingleCellExperiment::SingleCellExperiment].
#' @param path Output `.cytome` path.
#' @param assay Assay to store as the main counts matrix (default `"counts"`).
#' @param main_modality Modality name for the main assay (default `"RNA"`).
#' @param alt_modalities `altExp` names to write as extra modalities. `NULL`
#'   (default) writes all of them; `character(0)` writes none.
#' @param layers Additional assays of the main experiment to carry.
#'   `FALSE` (default) writes none, `TRUE` writes all, or name them. Off by
#'   default because normalized values are recomputable from the counts and
#'   cost more space than the counts do; turn it on when the point is to hand
#'   somebody your object rather than to archive the measurements.
#' @param graphs Write `colPairs` to the shared `graph_edges` table
#'   (default `TRUE`). A neighbour graph is the expensive thing to recompute,
#'   and carrying it means the recipient reproduces the same clusters.
#' @param dtype `"float32"` (default) or `"float64"`.
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
  if (file.exists(path)) file.remove(path)
  con <- DBI::dbConnect(RSQLite::SQLite(), path)
  on.exit(DBI::dbDisconnect(con))
  .cytome_write_schema(con)

  A <- SummarizedExperiment::assay(x, assay)               # features x cells
  bc <- colnames(x); if (is.null(bc)) bc <- paste0("cell", seq_len(ncol(x)))
  fid <- rownames(x); if (is.null(fid)) fid <- paste0("gene", seq_len(nrow(x)))
  sym <- .detect_feature_symbol(as.data.frame(SummarizedExperiment::rowData(x)), nrow(x))

  DBI::dbWriteTable(con, "cells",
    data.frame(cell_idx = seq_along(bc) - 1L, barcode = as.character(bc),
               sample_id = NA_character_, n_fragments = 0L, stringsAsFactors = FALSE),
    append = TRUE, row.names = FALSE)
  main_tbl <- .write_feature_table(con, main_modality, fid, sym)
  .cytome_write_matrix(con, paste0(main_modality, "_counts"),
                       Matrix::t(methods::as(A, "CsparseMatrix")),
                       "cells", main_tbl, dtype, chunk_size, compression)

  ## Additional assays of the MAIN experiment -- logcounts and friends. Off by
  ## default: for an archive they are recomputable and cost more space than the
  ## counts (a scaled float matrix compresses far worse than integers). For a
  ## handoff they are the difference between sending your object and sending
  ## counts plus an instruction to guess the parameters.
  for (a in .resolve_layers(layers, setdiff(SummarizedExperiment::assayNames(x), assay))) {
    .cytome_write_matrix(con, paste0(main_modality, "_", a),
                         Matrix::t(methods::as(SummarizedExperiment::assay(x, a),
                                               "CsparseMatrix")),
                         "cells", main_tbl, dtype, chunk_size, compression)
  }

  ## altExps as additional modalities. Previously they were silently dropped,
  ## so a multimodal object round-tripped as RNA only and nothing said so.
  alts <- SingleCellExperiment::altExpNames(x)
  if (is.null(alt_modalities)) alt_modalities <- alts
  for (m in intersect(alt_modalities, alts)) {
    ae <- SingleCellExperiment::altExp(x, m)
    aid <- rownames(ae); if (is.null(aid)) aid <- paste0(m, seq_len(nrow(ae)))
    tbl <- .write_feature_table(con, m, aid,
                                .detect_feature_symbol(
                                  as.data.frame(SummarizedExperiment::rowData(ae)), nrow(ae)))
    am <- SummarizedExperiment::assay(ae, 1L)
    .cytome_write_matrix(con, paste0(toupper(m), "_counts"),
                         Matrix::t(methods::as(am, "CsparseMatrix")),
                         "cells", tbl, dtype, chunk_size, compression)
  }

  .write_obs_cols(con, as.data.frame(SingleCellExperiment::colData(x)))
  if (embeddings && length(SingleCellExperiment::reducedDims(x)))
    .write_embeddings(con, as.list(SingleCellExperiment::reducedDims(x)))
  if (graphs && length(SingleCellExperiment::colPairNames(x)))
    .write_graphs(con, stats::setNames(
      lapply(SingleCellExperiment::colPairNames(x),
             function(g) SingleCellExperiment::colPair(x, g, asSparse = TRUE)),
      SingleCellExperiment::colPairNames(x)))
  DBI::dbExecute(con, "INSERT OR REPLACE INTO _manifest(key,value) VALUES('n_cells',?)",
                 params = list(as.character(ncol(x))))
  DBI::dbExecute(con, "PRAGMA wal_checkpoint(TRUNCATE);")
  invisible(path)
}

## Read one assay's matrix out of a Seurat object without requiring the full
## Seurat package. GetAssayData() dispatches through Seurat's S4 method, which
## calls UpdateSlots() and therefore needs Seurat itself -- deposited v4 objects
## then fail on a machine that only has SeuratObject. Try the supported call,
## fall back to the slot.
.seurat_layer <- function(obj, assay, layer) {
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
#' @param layer Layer/slot to store (default `"counts"`).
#' @param alt_assays Other assays to write as extra modalities. `NULL` (default)
#'   writes all of them; `character(0)` writes none.
#' @param layers Other layers of the main assay (`data`, `scale.data`) to
#'   carry. `FALSE` (default), `TRUE`, or name them. See the
#'   `SingleCellExperiment` method for why this is off by default.
#' @param graphs Write the `graphs` slot to the shared `graph_edges` table
#'   (default `TRUE`).
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
  A <- .seurat_layer(x, assay, layer)                      # features x cells
  if (file.exists(path)) file.remove(path)
  con <- DBI::dbConnect(RSQLite::SQLite(), path)
  on.exit(DBI::dbDisconnect(con))
  .cytome_write_schema(con)

  bc <- colnames(A); fid <- rownames(A)
  rmeta <- tryCatch(as.data.frame(x[[assay]][[]]), error = function(e) NULL)
  DBI::dbWriteTable(con, "cells",
    data.frame(cell_idx = seq_along(bc) - 1L, barcode = as.character(bc),
               sample_id = NA_character_, n_fragments = 0L, stringsAsFactors = FALSE),
    append = TRUE, row.names = FALSE)
  main_tbl <- .write_feature_table(con, assay, fid,
                                   .detect_feature_symbol(rmeta, length(fid)))
  .cytome_write_matrix(con, paste0(toupper(assay), "_counts"),
                       Matrix::t(methods::as(A, "CsparseMatrix")),
                       "cells", main_tbl, dtype, chunk_size, compression)

  ## Other layers of the MAIN assay (data, scale.data). See the SCE method for
  ## why this is off by default.
  for (lyr in .resolve_layers(layers, .seurat_layers(x, assay, layer))) {
    lm <- tryCatch(.seurat_layer(x, assay, lyr), error = function(e) NULL)
    if (is.null(lm) || !length(lm)) next
    .cytome_write_matrix(con, paste0(toupper(assay), "_", lyr),
                         Matrix::t(methods::as(lm, "CsparseMatrix")),
                         "cells", main_tbl, dtype, chunk_size, compression)
  }

  ## Other assays. They used to be dropped silently: a Seurat object with RNA
  ## and ATAC round-tripped as RNA only, and the file gave no sign that half of
  ## it was missing.
  others <- setdiff(names(x@assays), assay)
  if (is.null(alt_assays)) alt_assays <- others
  for (m in intersect(alt_assays, others)) {
    am <- tryCatch(.seurat_layer(x, m, layer), error = function(e) NULL)
    if (is.null(am)) next
    tbl <- .write_feature_table(con, m, rownames(am),
             .detect_feature_symbol(tryCatch(as.data.frame(x[[m]][[]]),
                                             error = function(e) NULL), nrow(am)))
    .cytome_write_matrix(con, paste0(toupper(m), "_counts"),
                         Matrix::t(methods::as(am, "CsparseMatrix")),
                         "cells", tbl, dtype, chunk_size, compression)
  }

  .write_obs_cols(con, x[[]])
  reds <- SeuratObject::Reductions(x)
  if (embeddings && length(reds)) .write_embeddings(con, stats::setNames(
    lapply(reds, function(r) SeuratObject::Embeddings(x, r)), reds))
  grs <- tryCatch(names(x@graphs), error = function(e) NULL)
  if (graphs && length(grs)) .write_graphs(con, stats::setNames(
    lapply(grs, function(g) x@graphs[[g]]), grs))
  DBI::dbExecute(con, "INSERT OR REPLACE INTO _manifest(key,value) VALUES('n_cells',?)",
                 params = list(as.character(ncol(A))))
  DBI::dbExecute(con, "PRAGMA wal_checkpoint(TRUNCATE);")
  invisible(path)
}
