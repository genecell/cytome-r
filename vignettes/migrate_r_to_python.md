# Migrating SCE / Seurat → cytome → PIASO (Python)

The native R writer, `write_cytome()`, lets R users export a `SingleCellExperiment` or `Seurat`
object to a `.cytome` file **with zero Python dependency**, then analyse it in Python with the
PIASO + cytome ecosystem. A store written by R is read back by Python's `cytome.to_anndata()`
in the package's tests (`tests/testthat/test-python.R`, which runs wherever a Python with cytome
is found), and every store Python or PIASO's importer writes is read by R
(`tests/testthat/fixtures`).

## R: write the object

```r
library(cytome)

## --- A public Seurat object (10x pbmc3k) ---
## download once:  https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
##   tar xzf pbmc3k_filtered_gene_bc_matrices.tar.gz
counts <- Seurat::Read10X("filtered_gene_bc_matrices/hg19/")
pbmc   <- SeuratObject::CreateSeuratObject(counts, project = "pbmc3k")
write_cytome(pbmc, "pbmc3k.cytome")                       # counts, cell metadata, reductions, graphs

## --- A SingleCellExperiment ---
## sce <- scRNAseq::ZeiselBrainData()
write_cytome(sce, "brain.cytome",
             assay = "counts", compression = "zstd",      # or "lz4" / "zlib"
             layers = TRUE)                               # logcounts and other assays too
```

What is written:

- the counts matrix (`RNA_counts`, chunked and compressed), and with `layers = TRUE` the other
  assays or layers of the main modality, each placed by feature and cell name;
- the other modalities cytome has a feature table for (ATAC, GA, tiles) from `altExps` or extra
  assays; others (SCT, ADT, ERCC) are named in a warning and not written;
- every `colData` / `obj[[]]` column in the `cells` table, with its type: factors keep their
  levels, logicals and integers stay logical and integer;
- `reducedDims` / Seurat `Reductions` as embeddings (`X_pca`, `X_umap`, ...), and `colPairs` /
  Seurat graphs as Python stores graphs, which `to_anndata()` returns in `obsp`.

cytome never stores values that were not measured. A layer over only some of the features (a
Seurat `scale.data` over the variable genes) is not written: a warning names it, and
`ScaleData()` recomputes it after reading. An assay, layer or graph over only some of the cells
(an assay measured in some batches only) is not written either: a warning names it, and to keep
it you subset the object to those cells and write that as a separate store.

The file is built beside the destination and moved into place only when it is complete, so a
write that fails leaves an existing file as it was.

## Python: analyse with PIASO

```python
import cytome, piaso
ds = cytome.open("pbmc3k.cytome")          # the R-written file, no conversion step
piaso.pp.normalize_log1p(ds, key_added="log1p", save_layer=True)   # INFOG/log1p normalisation
# … piaso.tl.score / runGDR / cosg / leiden, piaso.pl.* — the full streaming pipeline
print(ds.n_cells, list(ds.cells.to_pandas().columns))   # R cell metadata is visible in Python
ds.close()
```

## Verified end-to-end

A 300-cell × 200-gene Seurat object (3 planted groups) written by the R writer opens in Python,
exposes its `true_group` metadata, and runs `piaso.pp.normalize_log1p` (writes `RNA_log1p`) —
i.e. **Seurat → native R writer → .cytome → PIASO (Python)** works end to end.
