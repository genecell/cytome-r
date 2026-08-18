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
| normalization, dimensionality reduction, clustering | via PIASO | no — use Seurat / Bioconductor |

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
Needs `liblz4`, `libzstd`, `zlib` (e.g. from conda: `lz4-c zstd zlib`), plus `DBI`, `RSQLite`,
`Matrix`, `Rcpp`; `SingleCellExperiment` / `SeuratObject` are optional (for the object builders).

```r
Rcpp::compileAttributes(".")   # once, to generate RcppExports
# R CMD INSTALL .
```

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

# whole objects
sce <- read_cytome_sce("data.cytome")          # SingleCellExperiment (+ altExps, reducedDims)
so  <- read_cytome_seurat("data.cytome")        # Seurat

# out-of-core streaming, directly from the cytome (no full load):
x <- cytome_open("huge.cytome")
totals <- unlist(cytome_stream(x, "RNA_counts", function(m, s, e, i) Matrix::rowSums(m)))
cytome_close(x)
```

## Status (v1)
Read-only: matrices (lz4/zstd/zlib), obs/var, embeddings, fragments, KNN; SCE + Seurat builders;
chunk streaming. **Roadmap (v2):** writing (obs/embeddings/layers and beyond), a lazy/DelayedArray
backend for billion-nonzero stores, and an R port of PIASO streaming analysis on the cytome object.

The format contract is in [CYTOME_FORMAT_SPEC.md](CYTOME_FORMAT_SPEC.md); the cross-language
conformance test (`tests/testthat/test-conformance.R`) reads a cytome **written by Python** and
asserts a bit-for-bit match (both lz4 and zstd codecs).
