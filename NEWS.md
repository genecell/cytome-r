# cytome 0.1.4

Fixes from the code audit of 2026-10-27. Who was affected is said for each.

## Reading

* **Column-major matrices are read.** PIASO's importer and quantifier store
  ATAC peak and tile matrices by features (`matrix_csc_chunks`). They came
  back as an empty matrix with no error, streamed as no chunks and summed to
  zero through `cytome_delayed()`. Anyone who read a PIASO ATAC store in R
  analysed zeros. They are now read by `read_cytome_matrix()`,
  `read_cytome()` and `cytome_delayed()`, and `cytome_stream()` gains
  `along = "features"` to stream them by blocks of features. A matrix layout
  this version does not know, or a matrix recording values but holding no
  chunks, is now an error.
* **Every numeric dtype is read.** Only float32 and float64 were: a store with
  int32 counts (`from_10x_h5`, `from_cellranger`) could not be read at all,
  and an int64 embedding (a Visium `spatial`) was silently left out of
  `read_cytome()`. Each chunk is decoded in its own dtype (a `merge` of an
  int32 and a float32 store holds both), index arrays as int32 or int64 by
  their length, and every chunk is checked against what it records.
* **Fragments in compressed chunks are read.** `cytome_fragments()` read only
  the per-row tables, so stores from PIASO's importer, Python's streaming
  import or any `subset()` returned no fragments; and a query without an
  interval failed on every store. `dup_count` is `NA` for chunked fragments,
  which do not keep it, as the Python package reads them from cytome 0.3.9.
* **Column types come back.** Categoricals are factors with their levels in
  order, booleans logical, integers integer, as recorded in `_column_meta`;
  names such as `"col/with slash"` are no longer mangled. A value outside the
  recorded levels (the column was written since) is kept as an extra level,
  after the recorded ones, with a warning, as the Python package reads it
  from cytome 0.3.9.
* **`read_cytome()` takes as layers only matrices of the main assay's shape**
  (a `.raw` over more genes or a wide sparse obsm aborted it) and names what it
  leaves out. Feature embeddings and feature graphs are no longer offered as
  cell ones; a graph or embedding that fails to read is named in a warning.
  With `delayed = TRUE` every assay is backed by the file.
* **Stores in read-only folders open.** A store in WAL mode in a folder that is
  not writable failed with "attempt to write a readonly database"; it is now
  opened immutable when its WAL is empty, and refused with an explanation when
  committed changes are still in it. The bundled `reference.cytome` is no
  longer in WAL mode, so the examples no longer write `-wal` and `-shm` files
  into the installed library. Readers wait up to 60 s for a writer's lock.
* **`cytome_delayed()` is linear.** A block of 8,000 cells took 47 s, growing
  with the square of the cells; it takes 0.15 s.

## Writing

* **Values are placed by name, and cytome never stores values that were not
  measured.** A Seurat `scale.data` over the variable genes was written by
  position, so its values were stored under the wrong genes; an assay or a
  graph over only some of the cells (an assay measured in some batches only)
  was written by position too, so its values went to the wrong cells. Such a
  layer is now not written, with a warning that names it, says how many of
  the features it covers, and says to recompute it after reading
  (`ScaleData()`, or `SCTransform()` for an SCT assay): stored over all the
  features, the zeros of the features it lacks would read as values, and in
  `scale.data` 0 is the mean. An assay, layer or graph over only some of the
  cells is not written either, with a warning that names it, says how many of
  the cells it covers, and says to subset the object to those cells and write
  that as a separate store: its zeros would read as cells without counts, or
  without neighbours. The main assay must cover every cell. Matrices over all
  the features and cells are placed by name, whatever their order, and an
  embedding keeps `NaN` for the cells it lacks. Anyone who wrote
  `layers = TRUE` from a scaled Seurat object, or an object with an assay over
  some of its cells, should rewrite the file.
* **A failed write leaves the old file.** The old file was deleted first and a
  failure left a partial store. The store is now built beside the destination
  in one transaction and moved into place when complete.
* **Column types are kept**: factors with their levels, logicals and integers,
  recorded as the Python package reads them.
* **Graphs are written as the Python package stores them**, so Python's
  `to_anndata()` returns them in `obsp`; a colPair named `knn` read as an
  empty 0 x 0 graph there.
* **SCT, ADT and ERCC assays no longer fail the write**: they are skipped with
  a warning (cytome has no feature table for them). Split Seurat v5 layers are
  joined. GA and tiles modalities can be written.
* `dtype` is checked (`"int32"` now works, for whole numbers), and the writer
  is faster: a 20,000-cell object with five columns and a 15-NN graph took
  13.1 s to write those, 0.13 s now.

## Codecs

* The lz4 size header is read as little-endian on every platform, and blobs
  past the formats' 2 GB limits are refused by name instead of overflowing.

# cytome 0.1.3

* **The tests pass on arm64 macOS and Windows.** A test compared a float64
  graph with values read from a decimal CSV, exactly. Where R parses decimals
  without extended precision, some values came out one unit in the last place
  apart, so the check failed and no binary was built for those platforms. The
  expected values are now hex floats, which R parses exactly everywhere. The
  package code is unchanged from 0.1.2.

# cytome 0.1.2

* **Stores written without the Python `zstandard` package are read.** Without
  that optional package, the Python writer compresses with zlib but still labels
  each blob `zstd`. The Python reader decides by the bytes; this package trusted
  the label and failed with "zstd frame content size unknown". A blob labelled
  `zstd` is now decoded as what its bytes are (zstd, zlib or lz4), exactly as in
  Python; `lz4` and `zlib` labels are trusted, as there. A zstd frame that does
  not record its size is decoded in steps instead of refused, and a blob stored
  uncompressed (`none`) passes through.

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
