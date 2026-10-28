# feature id column for dimnames (genes: gene_id/symbol; peaks: peak_id; tiles: tile_id)
.feature_ids <- function(var) {
  for (col in c("gene_id", "peak_id", "tile_id", "protein_id", "var_name", "symbol", "id", "name"))
    if (col %in% names(var)) return(as.character(var[[col]]))
  as.character(seq_len(nrow(var)))
}
.cell_ids <- function(obs) if ("barcode" %in% names(obs)) as.character(obs$barcode) else as.character(seq_len(nrow(obs)))

## The cells' names: the column the Python package recorded as the AnnData's obs_names
## (`_anndata_obs_index`, set when obs had its own barcode column that differs from them),
## else the barcodes.
.cell_names <- function(x, obs) {
  col <- .metadata_value(x, "_anndata_obs_index")
  if (is.character(col) && length(col) == 1L && col %in% names(obs)) as.character(obs[[col]])
  else .cell_ids(obs)
}

## A value of the store's `_metadata` (JSON), or NULL.
.metadata_value <- function(x, key) {
  if (!.has_table(x, "_metadata")) return(NULL)
  v <- DBI::dbGetQuery(x$con, "SELECT value FROM _metadata WHERE key = ?", params = list(key))$value
  if (!length(v)) return(NULL)
  tryCatch(jsonlite::fromJSON(v[1]), error = function(e) NULL)
}

.count_modalities <- function(mats) sub("_counts$", "", mats$matrix_name[grepl("_counts$", mats$matrix_name)])

## The main matrix must exist; when it does not, say which modalities the store has.
.check_main_matrix <- function(mats, main_name, main_modality) {
  if (main_name %in% mats$matrix_name) return(invisible())
  mods <- .count_modalities(mats)
  stop("cytome: no '", main_name, "' matrix",
       if (length(mods)) paste0("; this store has ", paste0(mods, "_counts", collapse = ", "),
                                ": choose one with main_modality = \"", mods[1], "\"")
       else paste0("; its matrices are ", paste(mats$matrix_name, collapse = ", ")))
}

## Matrices of `modality` that are not its counts and not another modality's:
## RNA_logcounts yes, RNA_counts no, ATAC_counts no. raw_X is the AnnData .raw,
## not a layer.
.extra_layer_names <- function(mats, modality) {
  all_mods <- .count_modalities(mats)
  nm <- mats$matrix_name[startsWith(mats$matrix_name, paste0(modality, "_"))]
  nm <- setdiff(nm, paste0(modality, "_counts"))
  nm[!sub(paste0("^", modality, "_"), "", nm) %in% c(all_mods, "raw_X")]
}

## Which of the main modality's other matrices are layers: the same cells and the same
## features as its counts. A matrix of another width (an AnnData .raw over more genes, a
## sparse obsm stored as a matrix) has other columns, and labelling them with the genes
## aborted read_cytome() on the dimnames. Returns list(layers, skipped).
.layer_plan <- function(x, mats, modality, main) {
  layers <- character(0); skipped <- character(0)
  for (nm in mats$matrix_name[startsWith(mats$matrix_name, paste0(modality, "_"))]) {
    if (nm == main$name || nm %in% paste0(.count_modalities(mats), "_counts")) next
    info <- .matrix_info(x, nm)
    same <- identical(info$n_rows, main$n_rows) && identical(info$n_cols, main$n_cols) &&
      identical(info$col_entity, main$col_entity)
    if (same && nm %in% .extra_layer_names(mats, modality)) layers <- c(layers, nm)
    else skipped <- c(skipped, sprintf("%s (%s x %s)", nm, info$n_rows, info$n_cols))
  }
  list(layers = layers, skipped = skipped)
}

## The embeddings of cells (a var embedding, e.g. varm PCs, has as many rows as features).
.cell_embeddings <- function(x) {
  if (!.has_table(x, "embedding_meta")) return(character(0))
  if (!"entity" %in% .table_columns(x, "embedding_meta")) return(cytome_embeddings(x))
  DBI::dbGetQuery(x$con,
    "SELECT array_name FROM embedding_meta WHERE entity = 'cells' OR entity IS NULL")$array_name
}

## What read_cytome() could not bring across, said once.
.report_skipped <- function(skipped, failed) {
  if (length(skipped))
    message("cytome: not read into the object, as their shape is not the main assay's: ",
            paste(skipped, collapse = ", "), ". Read them with read_cytome_matrix().")
  if (length(failed))
    warning("cytome: could not be read: ", paste(failed, collapse = "; "), call. = FALSE)
}

#' Read a .cytome into a SingleCellExperiment
#'
#' Maps the cytome's multimodal structure to an SCE: the main modality's counts to the main
#' assay, other modalities to `altExps`, embeddings to `reducedDims`, the cells table to
#' `colData`, and the feature table to `rowData`. Native R, zero Python.
#'
#' @param path Path to a `.cytome` file.
#' @param main_modality Modality for the main assay (default `"RNA"`).
#' @param alt_modalities Modalities to add as `altExps`. `NULL` (default) = all other `_counts`
#'   matrices; `character(0)` = none.
#' @param embeddings If `TRUE` (default) load embeddings into `reducedDims`.
#' @param delayed If `TRUE`, every assay (main, layers and altExps) is a DelayedArray.
#' @noRd
#' @return A [SingleCellExperiment::SingleCellExperiment].
.read_cytome_sce <- function(path, main_modality = "RNA", alt_modalities = NULL,
                             embeddings = TRUE, delayed = FALSE) {
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE))
    stop("cytome: install 'SingleCellExperiment' to use read_cytome()")
  x <- cytome_open(path); on.exit(cytome_close(x))
  mats <- cytome_matrices(x)
  obs <- cytome_obs(x, "cells")
  cells <- .cell_names(x, obs)
  failed <- character(0)

  main_name <- paste0(main_modality, "_counts")
  .check_main_matrix(mats, main_name, main_modality)
  main <- .matrix_info(x, main_name)
  var <- cytome_var(x, main_name)
  ## Backed by the file when delayed: nothing is read until something asks for values.
  assay_of <- function(nm, features) {
    M <- if (delayed) cytome_delayed(path, nm) else read_cytome_matrix(x, nm)   # features x cells
    dimnames(M) <- list(features, cells)
    M
  }
  M <- assay_of(main_name, .feature_ids(var))

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = M),
    colData = S4Vectors::DataFrame(obs, row.names = cells, check.names = FALSE),
    rowData = S4Vectors::DataFrame(var, row.names = .feature_ids(var), check.names = FALSE))

  if (is.null(alt_modalities)) alt_modalities <- setdiff(.count_modalities(mats), main_modality)
  skipped <- character(0)
  for (m in alt_modalities) {
    nm <- paste0(m, "_counts")
    if (!(nm %in% mats$matrix_name)) next
    info <- .matrix_info(x, nm); av <- cytome_var(x, nm)
    if (!identical(info$n_rows, main$n_rows) || nrow(av) != info$n_cols) {
      skipped <- c(skipped, sprintf("%s (%s x %s)", nm, info$n_rows, info$n_cols))
      next
    }
    SingleCellExperiment::altExp(sce, m) <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = assay_of(nm, .feature_ids(av))),
      rowData = S4Vectors::DataFrame(av, row.names = .feature_ids(av), check.names = FALSE))
  }

  ## Non-count matrices of the main modality with its shape come back as extra assays, so
  ## a file written with layers = TRUE round-trips rather than silently losing the
  ## normalization.
  plan <- .layer_plan(x, mats, main_modality, main)
  for (nm in plan$layers)
    SummarizedExperiment::assay(sce, sub(paste0("^", main_modality, "_"), "", nm)) <-
      assay_of(nm, .feature_ids(var))

  ## colPair<- goes through SelfHits, whose from/to must be plain integers.
  ## Handing it a dgCMatrix lets the coercion produce doubles and the object
  ## fails validity, so build the SelfHits explicitly.
  for (g in cytome_graphs(x)) {
    gm <- tryCatch(cytome_graph(x, g, n = length(cells)), error = function(e) {
      failed <<- c(failed, paste0("graph ", g, " (", conditionMessage(e), ")")); NULL })
    if (is.null(gm)) next
    tg <- methods::as(gm, "TsparseMatrix")
    if (!length(tg@x)) next
    hits <- S4Vectors::SelfHits(from = as.integer(tg@i) + 1L,
                                to = as.integer(tg@j) + 1L,
                                nnode = length(cells))
    S4Vectors::mcols(hits)$value <- as.numeric(tg@x)
    SingleCellExperiment::colPair(sce, g) <- hits
  }

  if (embeddings) {
    for (e in .cell_embeddings(x)) {
      em <- tryCatch(cytome_embedding(x, e), error = function(err) {
        failed <<- c(failed, paste0("embedding ", e, " (", conditionMessage(err), ")")); NULL })
      if (is.null(em)) next
      if (nrow(em) != ncol(sce)) {
        skipped <- c(skipped, sprintf("embedding %s (%d rows)", e, nrow(em)))
        next
      }
      rownames(em) <- cells
      SingleCellExperiment::reducedDim(sce, e) <- em
    }
  }
  .report_skipped(c(skipped, plan$skipped), failed)
  sce
}
