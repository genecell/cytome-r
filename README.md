# cytome

**Native R reader, writer and streaming interface for the `.cytome` single-cell
store.**

A `.cytome` file is a plain **SQLite** database, so R reads it directly with
**no Python dependency**: the relational tables (cells, genes, peaks,
embeddings, graphs) through `DBI`, and the chunked, compressed CSR count
matrices through a small Rcpp layer (lz4-block / zstd / zlib). It builds a
**SingleCellExperiment** or a **Seurat** object, writes either one back out,
and streams the matrix in chunks for work that must not load it whole.

```r
install.packages("cytome", repos = "https://genecell.r-universe.dev")
```

`cytome` is part of the **[PIASO](https://piaso.org)** ecosystem — *Precise
Integrative Analysis of Single-cell Omics* — a toolkit of companion packages
for single-cell and spatial omics covering normalization, dimensionality
reduction, marker gene identification, cell type annotation, gene-set scoring,
ligand–receptor analysis and more, together with the `.cytome` format that this
package reads and writes. Neither half requires the other: a `.cytome` file is
a plain SQLite database, and this reader has no Python dependency at all.

## What this package does, and what it does not

The R package is a **reader, writer and streaming interface for the format**.
It deliberately does not reimplement analysis: once your data is a Seurat
object or a SingleCellExperiment, use the tools you already use. Datasets are
*built* on the Python side and *analysed* on either.

| | Python `cytome` | R `cytome` |
|---|---|---|
| read a `.cytome` | yes | yes |
| write a `.cytome` | yes | yes — `write_cytome()`, RNA and ATAC modalities |
| normalized layers | yes | yes — `layers = TRUE` (off by default) |
| graphs / KNN | yes | yes — via the shared `graph_edges` table |
| convert to an in-memory object | AnnData | Seurat, SingleCellExperiment |
| chunked streaming | yes | yes — `cytome_stream()` |
| out-of-core matrix | yes | `cytome_delayed()` / `read_cytome(delayed = TRUE)`, chunk-aligned |
| build from Cell Ranger / 10x h5 | yes | no — use Python |
| merge, subset, filter cells | yes | no — use Python |
| end-to-end analysis — normalization, dimensionality reduction, clustering, annotation, gene-set scoring, cell–cell interaction | via [PIASO](https://piaso.org) | no — use PIASO, Seurat or Bioconductor |

The format itself is specified in the Python repository —
[`CYTOME_FORMAT_SPEC.md`](https://github.com/genecell/cytome/blob/main/CYTOME_FORMAT_SPEC.md),
which is the reference implementation. This package implements that spec and
is tested against a reference file written by Python, **in both directions**;
it deliberately does not carry its own copy of the document, because a second
copy is where drift starts.

## Reading and writing

One entry point each way. `read_cytome()` chooses what you get back;
`write_cytome()` is a generic that dispatches on what you hand it.

```r
library(cytome)

sce <- read_cytome("data.cytome")                          # SingleCellExperiment
so  <- read_cytome("data.cytome", as = "Seurat")           # Seurat
x   <- read_cytome("data.cytome", as = "cytome")           # the open handle

write_cytome(so,  "out.cytome")     # dispatches on the object -- no _seurat suffix
write_cytome(sce, "out.cytome")
```

Alternative modalities travel with the object: a Seurat with `RNA` and `ATAC`
assays, or an SCE with an `ATAC` `altExp`, round-trips with both.

## Out-of-core

```r
sce <- read_cytome("big.cytome", delayed = TRUE)
SummarizedExperiment::assay(sce)          # a DelayedArray, read on demand
DelayedArray::colSums(SummarizedExperiment::assay(sce))
```

Extraction is chunk-aligned — only the storage chunks overlapping the request
are read, and `chunkdim()` reports the storage geometry so DelayedArray's
blocks land on chunk boundaries. Peak memory is set by the block size, not the
matrix. For data that fits in memory, leave `delayed = FALSE`; the reason to
use this is data that does not.

For explicit chunk-wise work:

```r
x <- read_cytome("big.cytome", as = "cytome")
totals <- cytome_stream(x, "RNA_counts",
                        function(chunk, i0, i1, k) Matrix::rowSums(chunk))
per_cell <- Reduce(`+`, totals)          # stream returns one result per chunk
cytome_close(x)
```

## Install

```r
install.packages("cytome", repos = "https://genecell.r-universe.dev")   # binaries
# development version (needs a C++ toolchain and the headers below):
# remotes::install_github("genecell/cytome-r")
```

The codecs a `.cytome` is compressed with are required, not optional:

```
Debian/Ubuntu  apt-get install liblz4-dev libzstd-dev zlib1g-dev
Fedora/RHEL    dnf install lz4-devel libzstd-devel zlib-devel
macOS          brew install lz4 zstd
conda          conda install -c conda-forge lz4-c zstd zlib
```

`configure` finds them through `pkg-config`, accepts `CYTOME_CFLAGS` /
`CYTOME_LIBS` as an override, and stops with these instructions rather than
failing part-way through compiling. `SingleCellExperiment`, `SeuratObject` and
the `DelayedArray` stack are suggested, not required — install the ones whose
objects you actually want.

## Usage
```r
library(cytome)

x <- cytome_open("data.cytome")
x                                  # print summary: cells, matrices, embeddings
cytome_matrices(x)                 # catalogue
obs   <- cytome_obs(x, "cells")    # colData
M     <- read_cytome_matrix(x, "RNA_counts")   # features x cells dgCMatrix
umap  <- cytome_embedding(x, "X_umap")
frags <- cytome_fragments(x, "chr1", start = 1e6, end = 1.1e6)
cytome_close(x)

# whole objects — one entry point each way
sce <- read_cytome("data.cytome")               # SingleCellExperiment (+ altExps, reducedDims)
so  <- read_cytome("data.cytome", as = "Seurat")
write_cytome(so, "out.cytome")                  # dispatches on the object

# graphs stored in the file
cytome_graphs(x)
knn <- cytome_graph(x, "RNA_nn")
```

## Related

- **[PIASO](https://piaso.org)** — the ecosystem this format was built for:
  INFOG normalization, marker-gene-guided dimensionality reduction (GDR),
  gene-set scoring, cell type prediction, ligand–receptor analysis, and the
  companion packages around them (COSG, SCALAR, LARIS, PIASOmarkerDB).
- **[cytome (Python)](https://github.com/genecell/cytome)** — the reference
  implementation, and where `.cytome` files are usually built. It has the
  importers, the merge / subset / filter operations, and the
  [format specification](https://github.com/genecell/cytome/blob/main/CYTOME_FORMAT_SPEC.md)
  this package implements.
- **Documentation and tutorials** — <https://piaso.org>, including a page on
  converting between `.cytome`, AnnData, Seurat and SingleCellExperiment.

## Correctness

The package ships a reference `.cytome` **written by the Python
implementation** plus the expected values, and its tests assert it reads them
bit-for-bit — both codecs, zstd for RNA and lz4 for ATAC. CI runs the reverse
too: a file written by R, read back by Python, compared against the same
expectations. For two implementations of one format, testing only one
direction leaves the other free to rot.

## Development

`cytome` and the PIASO ecosystem are developed in the
**[Gord Fishell Laboratory](https://fishelllab.hms.harvard.edu)** at Harvard
Medical School and the Broad Institute.

Maintainer: Min Dai (<dai@broadinstitute.org>). Issues and feature requests:
<https://github.com/genecell/cytome-r/issues>.
