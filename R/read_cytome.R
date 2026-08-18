#' Read a `.cytome` file
#'
#' One entry point. `as` chooses what you get back, mirroring
#' `anndataR::read_h5ad(path, as = ...)`:
#'
#' * `"SingleCellExperiment"` (default) -- main modality as the main assay,
#'   other modalities as `altExps`, embeddings as `reducedDims`.
#' * `"Seurat"` -- main modality as the default assay, other modalities as
#'   extra assays, embeddings as `DimReduc`s.
#' * `"cytome"` -- the open handle itself, for chunk-wise work with
#'   [cytome_stream()]. Close it with [cytome_close()] when done.
#'
#' @param path Path to a `.cytome` file.
#' @param as What to return. See above.
#' @param main_modality Modality for the main assay (default `"RNA"`).
#' @param alt_modalities Modalities to bring across as `altExps` / extra
#'   assays. `NULL` (default) = all other `_counts` matrices;
#'   `character(0)` = none.
#' @param embeddings If `TRUE` (default) load embeddings.
#' @param delayed If `TRUE`, back the main assay with a
#'   [DelayedArray::DelayedArray] that reads chunks from the file on demand,
#'   so the matrix is never held in memory. `SingleCellExperiment` only, and
#'   experimental -- see [cytome_delayed()].
#' @param ... Passed to the underlying reader.
#' @return A `SingleCellExperiment`, a `Seurat` object, or a `cytome` handle.
#' @seealso [write_cytome()], [cytome_stream()], [cytome_delayed()]
#' @examples
#' ref <- system.file("extdata", "reference.cytome", package = "cytome")
#' if (nzchar(ref)) {
#'   x <- read_cytome(ref, as = "cytome")
#'   print(x)
#'   cytome_close(x)
#' }
#' @export
read_cytome <- function(path,
                        as = c("SingleCellExperiment", "Seurat", "cytome"),
                        main_modality = "RNA", alt_modalities = NULL,
                        embeddings = TRUE, delayed = FALSE, ...) {
  as <- match.arg(as)
  if (as == "cytome") {
    if (delayed) warning("cytome: 'delayed' is ignored when as = \"cytome\"")
    return(cytome_open(path))
  }
  if (as == "Seurat") {
    if (delayed)
      stop("cytome: delayed = TRUE is only implemented for ",
           "as = \"SingleCellExperiment\". Seurat assays require an in-memory ",
           "matrix.")
    return(.read_cytome_seurat(path, main_modality = main_modality,
                               alt_modalities = alt_modalities,
                               embeddings = embeddings, ...))
  }
  .read_cytome_sce(path, main_modality = main_modality,
                   alt_modalities = alt_modalities,
                   embeddings = embeddings, delayed = delayed, ...)
}
