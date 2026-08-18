# feature id column for dimnames (genes: gene_id/symbol; peaks: peak_id)
.feature_ids <- function(var) {
  for (col in c("gene_id", "peak_id", "symbol", "id", "name")) if (col %in% names(var)) return(as.character(var[[col]]))
  as.character(seq_len(nrow(var)))
}
.cell_ids <- function(obs) if ("barcode" %in% names(obs)) as.character(obs$barcode) else as.character(seq_len(nrow(obs)))

.count_modalities <- function(mats) sub("_counts$", "", mats$matrix_name[grepl("_counts$", mats$matrix_name)])

## Matrices of `modality` that are not its counts and not another modality's:
## RNA_logcounts yes, RNA_counts no, ATAC_counts no.
.extra_layer_names <- function(mats, modality) {
  all_mods <- .count_modalities(mats)
  nm <- mats$matrix_name[startsWith(mats$matrix_name, paste0(modality, "_"))]
  nm <- setdiff(nm, paste0(modality, "_counts"))
  nm[!sub(paste0("^", modality, "_"), "", nm) %in% c(all_mods, "raw_X")]
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
#' @noRd
#' @return A [SingleCellExperiment::SingleCellExperiment].
.read_cytome_sce <- function(path, main_modality = "RNA", alt_modalities = NULL,
                             embeddings = TRUE, delayed = FALSE) {
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE))
    stop("cytome: install 'SingleCellExperiment' to use read_cytome()")
  x <- cytome_open(path); on.exit(cytome_close(x))
  mats <- cytome_matrices(x)
  obs <- cytome_obs(x, "cells")
  cells <- .cell_ids(obs)

  main_name <- paste0(main_modality, "_counts")
  if (!(main_name %in% mats$matrix_name)) stop("cytome: no '", main_name, "' matrix")
  var <- cytome_var(x, main_name)
  if (delayed) {
    ## Backed by the file: the assay is a DelayedArray and nothing is read
    ## until something asks for values.
    M <- cytome_delayed(path, main_name)
  } else {
    M <- read_cytome_matrix(x, main_name)                # features x cells
    dimnames(M) <- list(.feature_ids(var), cells)
  }

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = M),
    colData = S4Vectors::DataFrame(obs, row.names = cells),
    rowData = S4Vectors::DataFrame(var, row.names = .feature_ids(var)))

  if (is.null(alt_modalities)) alt_modalities <- setdiff(.count_modalities(mats), main_modality)
  for (m in alt_modalities) {
    nm <- paste0(m, "_counts")
    if (!(nm %in% mats$matrix_name)) next
    am <- read_cytome_matrix(x, nm); av <- cytome_var(x, nm)
    dimnames(am) <- list(.feature_ids(av), cells)
    SingleCellExperiment::altExp(sce, m) <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = am), rowData = S4Vectors::DataFrame(av, row.names = .feature_ids(av)))
  }

  ## Non-count matrices of the main modality come back as extra assays, so a
  ## file written with layers = TRUE round-trips rather than silently losing
  ## the normalization.
  for (nm in .extra_layer_names(mats, main_modality)) {
    lm <- read_cytome_matrix(x, nm)
    dimnames(lm) <- list(.feature_ids(var), cells)
    SummarizedExperiment::assay(sce, sub(paste0("^", main_modality, "_"), "", nm)) <- lm
  }

  ## colPair<- goes through SelfHits, whose from/to must be plain integers.
  ## Handing it a dgCMatrix lets the coercion produce doubles and the object
  ## fails validity, so build the SelfHits explicitly.
  for (g in cytome_graphs(x)) {
    gm <- tryCatch(cytome_graph(x, g, n = length(cells)), error = function(e) NULL)
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
    for (e in cytome_embeddings(x)) {
      em <- tryCatch(cytome_embedding(x, e), error = function(err) NULL)
      if (!is.null(em) && nrow(em) == ncol(sce)) {
        rownames(em) <- cells
        SingleCellExperiment::reducedDim(sce, e) <- em
      }
    }
  }
  sce
}
