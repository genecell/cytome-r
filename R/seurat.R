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
  obs <- cytome_obs(x, "cells"); cells <- .cell_ids(obs)

  main_name <- paste0(main_modality, "_counts")
  if (!(main_name %in% mats$matrix_name)) stop("cytome: no '", main_name, "' matrix")
  M <- read_cytome_matrix(x, main_name); var <- cytome_var(x, main_name)
  dimnames(M) <- list(make.unique(.feature_ids(var)), cells)
  md <- as.data.frame(obs); rownames(md) <- cells

  so <- SeuratObject::CreateSeuratObject(counts = M, assay = assay, meta.data = md)

  if (is.null(alt_modalities)) alt_modalities <- setdiff(.count_modalities(mats), main_modality)
  for (m in alt_modalities) {
    nm <- paste0(m, "_counts")
    if (!(nm %in% mats$matrix_name)) next
    am <- read_cytome_matrix(x, nm); av <- cytome_var(x, nm)
    dimnames(am) <- list(make.unique(.feature_ids(av)), cells)
    so[[m]] <- SeuratObject::CreateAssayObject(counts = am)
  }

  for (nm in .extra_layer_names(mats, main_modality)) {
    lm <- read_cytome_matrix(x, nm)
    dimnames(lm) <- list(make.unique(.feature_ids(var)), cells)
    lyr <- sub(paste0("^", main_modality, "_"), "", nm)
    so <- tryCatch({ SeuratObject::LayerData(so, assay = assay, layer = lyr) <- lm; so },
                   error = function(e) so)
  }

  for (g in cytome_graphs(x)) {
    gm <- tryCatch(cytome_graph(x, g, n = length(cells)), error = function(e) NULL)
    if (is.null(gm)) next
    dimnames(gm) <- list(cells, cells)
    so@graphs[[g]] <- tryCatch(SeuratObject::as.Graph(gm), error = function(e) NULL)
  }

  if (embeddings) {
    for (e in cytome_embeddings(x)) {
      em <- tryCatch(cytome_embedding(x, e), error = function(err) NULL)
      if (is.null(em) || nrow(em) != length(cells)) next
      rownames(em) <- cells
      key <- gsub("[^A-Za-z0-9]", "", e); if (key == "") key <- "dim"
      colnames(em) <- paste0(key, "_", seq_len(ncol(em)))
      so[[e]] <- SeuratObject::CreateDimReducObject(embeddings = em, key = paste0(key, "_"),
                                                    assay = assay)
    }
  }
  so
}
