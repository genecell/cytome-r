#' Read a .cytome into a Seurat object
#'
#' Builds a Seurat object directly from the cytome (native R): the main modality becomes the
#' default assay, other modalities become additional assays, and embeddings become `DimReduc`s.
#' Requires the (suggested) `SeuratObject` package.
#'
#' @param path Path to a `.cytome` file.
#' @param main_modality Modality for the default assay (default `"RNA"`).
#' @param assay Assay name for the main modality (default = `main_modality`).
#' @param alt_modalities Modalities to add as extra assays. `NULL` (default) = all other
#'   `_counts` matrices; `character(0)` = none.
#' @param embeddings If `TRUE` (default) add embeddings as `DimReduc`s.
#' @noRd
#' @return A `Seurat` object.
.read_cytome_seurat <- function(path, main_modality = "RNA", assay = main_modality,
                               alt_modalities = NULL, embeddings = TRUE) {
  if (!requireNamespace("SeuratObject", quietly = TRUE))
    stop("cytome: install 'SeuratObject' (or 'Seurat') to use read_cytome()")
  x <- cytome_open(path); on.exit(cytome_close(x))
  mats <- cytome_matrices(x)
  obs <- cytome_obs(x, "cells"); cells <- .cell_names(x, obs)
  if (anyDuplicated(cells)) {
    warning("cytome: ", sum(duplicated(cells)), " cell names repeat (barcodes shared by merged ",
            "libraries?); Seurat needs unique names, so they are made unique", call. = FALSE)
    cells <- make.unique(cells)
  }
  failed <- character(0); skipped <- character(0)

  main_name <- paste0(main_modality, "_counts")
  .check_main_matrix(mats, main_name, main_modality)
  main <- .matrix_info(x, main_name)
  M <- read_cytome_matrix(x, main_name); var <- cytome_var(x, main_name)
  ## Seurat replaces "_" in feature names with "-"; doing it here, before the names are
  ## unique, keeps every layer's names equal to the assay's.
  dimnames(M) <- list(.seurat_features(.feature_ids(var)), cells)
  md <- as.data.frame(obs, optional = TRUE); rownames(md) <- cells

  so <- SeuratObject::CreateSeuratObject(counts = M, assay = assay, meta.data = md)

  if (is.null(alt_modalities)) alt_modalities <- setdiff(.count_modalities(mats), main_modality)
  for (m in alt_modalities) {
    nm <- paste0(m, "_counts")
    if (!(nm %in% mats$matrix_name)) next
    info <- .matrix_info(x, nm); av <- cytome_var(x, nm)
    if (!identical(info$n_rows, main$n_rows) || nrow(av) != info$n_cols) {
      skipped <- c(skipped, sprintf("%s (%s x %s)", nm, info$n_rows, info$n_cols))
      next
    }
    am <- read_cytome_matrix(x, nm)
    dimnames(am) <- list(.seurat_features(.feature_ids(av)), cells)
    so[[m]] <- SeuratObject::CreateAssayObject(counts = am)
  }

  plan <- .layer_plan(x, mats, main_modality, main)
  features <- rownames(so[[assay]])
  for (nm in plan$layers) {
    lm <- read_cytome_matrix(x, nm)
    dimnames(lm) <- list(features, cells)
    lyr <- sub(paste0("^", main_modality, "_"), "", nm)
    so <- tryCatch({ SeuratObject::LayerData(so, assay = assay, layer = lyr) <- lm; so },
                   error = function(e) {
                     failed <<- c(failed, paste0("layer ", lyr, " (", conditionMessage(e), ")")); so })
  }

  for (g in cytome_graphs(x)) {
    gm <- tryCatch(cytome_graph(x, g, n = length(cells)), error = function(e) {
      failed <<- c(failed, paste0("graph ", g, " (", conditionMessage(e), ")")); NULL })
    if (is.null(gm)) next
    dimnames(gm) <- list(cells, cells)
    so@graphs[[g]] <- SeuratObject::as.Graph(gm)
  }

  if (embeddings) {
    for (e in .cell_embeddings(x)) {
      em <- tryCatch(cytome_embedding(x, e), error = function(err) {
        failed <<- c(failed, paste0("embedding ", e, " (", conditionMessage(err), ")")); NULL })
      if (is.null(em)) next
      if (nrow(em) != length(cells)) {
        skipped <- c(skipped, sprintf("embedding %s (%d rows)", e, nrow(em)))
        next
      }
      rownames(em) <- cells
      key <- gsub("[^A-Za-z0-9]", "", e); if (key == "") key <- "dim"
      colnames(em) <- paste0(key, "_", seq_len(ncol(em)))
      so[[e]] <- SeuratObject::CreateDimReducObject(embeddings = em, key = paste0(key, "_"),
                                                    assay = assay)
    }
  }
  .report_skipped(c(skipped, plan$skipped), failed)
  so
}

## Feature names as Seurat keeps them: no underscores, unique.
.seurat_features <- function(ids) make.unique(gsub("_", "-", ids, fixed = TRUE))
