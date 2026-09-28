# cytome 0.1.1

* **Graphs written by the Python package from cytome 0.3.6 are read.** Python
  cytome 0.3.6 stores graphs as compressed row chunks (`graph_meta`,
  `graph_chunks`) instead of one row per edge. 0.1.0 read only the per-edge
  `graph_edges` table, so on such a store `cytome_graphs()` was empty and
  `read_cytome()` returned a `SingleCellExperiment` or `Seurat` object without
  its neighbour graphs, with no error. `cytome_graphs()` now lists both forms
  and `cytome_graph()` reads both, keeping the stored shape. Graphs written
  from R, or by Python before 0.3.6, read as before.

# cytome 0.1.0

First release.

* Reads `.cytome` files natively in R, with no Python dependency. A `.cytome`
  is a plain SQLite database, so cells / genes / peaks, embeddings and graphs
  are read through DBI; the chunked, compressed CSR matrices are decoded by a
  small Rcpp layer (lz4-block / zstd / zlib).
* `read_cytome_sce()` and `read_cytome_seurat()` build a `SingleCellExperiment`
  or a `Seurat` object, carrying alternative modalities and embeddings across.
* `write_cytome_sce()` and `write_cytome_seurat()` go the other way, so a
  `.cytome` written in R opens in the Python package and vice versa.
* `cytome_stream()` iterates the matrix in chunks for out-of-core work, so peak
  memory is set by the chunk size rather than by the number of cells.
* A `configure` script locates lz4 / zstd / zlib through pkg-config, with an
  environment-variable override and a compile test, rather than assuming a
  particular installation prefix.

## API, before the first publication

* `read_cytome(path, as = c("SingleCellExperiment", "Seurat", "cytome"))` and
  `write_cytome(x, path)` replace `read_cytome_sce()`, `read_cytome_seurat()`,
  `write_cytome_sce()` and `write_cytome_seurat()`. `write_cytome()` is an S3
  generic, so the method follows from the object you are holding and other
  packages can add their own; `as` is validated by `match.arg()`.
* **Alternative modalities no longer disappear on write.** A Seurat object with
  `RNA` and `ATAC` assays, or an SCE with an `ATAC` `altExp`, used to
  round-trip as RNA only, silently. Both now survive, and a reverse
  conformance check (R writes, Python reads) asserts it.
* `cytome_delayed()` and `read_cytome(delayed = TRUE)` back the assay with a
  `DelayedArray` that reads from the file on demand. Experimental.
* Writing a modality R cannot build a feature table for is an error naming the
  supported ones, instead of a file Python cannot open.

## Out-of-core, layers and graphs

* **The `DelayedArray` seed is chunk-aligned.** It reads only the storage
  chunks overlapping the requested cells, and `chunkdim()` reports the storage
  geometry so DelayedArray's block grid lands on chunk boundaries rather than
  straddling them. The previous seed read the whole matrix and then subset,
  which made peak memory the matrix plus the block — worse than loading it
  once, and not what the documentation claimed.
* **`layers =` on both write methods.** `FALSE` (default), `TRUE`, or a
  character vector. Off by default because an archive and a handoff want
  opposite things; naming a layer that is absent lists the ones that are not.
  They come back as extra assays / layers on read.
* **Graphs travel.** A Seurat `graphs` slot or an SCE `colPairs` is written to
  the same `graph_edges` table the Python package uses, so a KNN graph built
  on either side is readable from the other. On by default; `graphs = FALSE`
  opts out. `cytome_graphs()` and `cytome_graph()` read them directly.
