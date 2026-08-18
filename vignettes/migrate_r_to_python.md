# Migrating SCE / Seurat → cytome → PIASO (Python)

The native R writer (`write_cytome_sce`, `write_cytome_seurat`) lets R users export a
`SingleCellExperiment` or `Seurat` object to a `.cytome` file **with zero Python dependency**, then
analyse it in Python with the PIASO + cytome ecosystem. The file is read bit-for-bit by Python
`cytome.open()` (cross-language write conformance is in `tests/testthat/test-write.R`).

## R: write the object

```r
library(cytome)

## --- A public Seurat object (10x pbmc3k) ---
## download once:  https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
##   tar xzf pbmc3k_filtered_gene_bc_matrices.tar.gz
counts <- Seurat::Read10X("filtered_gene_bc_matrices/hg19/")
pbmc   <- SeuratObject::CreateSeuratObject(counts, project = "pbmc3k")
write_cytome_seurat(pbmc, "pbmc3k.cytome")             # counts layer + cell metadata

## --- A SingleCellExperiment ---
## sce <- scRNAseq::ZeiselBrainData()
write_cytome_sce(sce, "brain.cytome",
                 assay = "counts", compression = "zstd",  # or "lz4" / "zlib"
                 embeddings = TRUE)                        # reducedDims -> X_<name>
```

What is written: the counts matrix (`RNA_counts`, chunked compressed CSR), every `colData` /
`obj[[]]` column (→ the `cells` table), and `reducedDims` / Seurat `Reductions` (→ embeddings
`X_pca`, `X_umap`, …). All read natively in Python.

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

A 300-cell × 200-gene Seurat object (3 planted groups) written by `write_cytome_seurat` opens in
Python, exposes its `true_group` metadata, and runs `piaso.pp.normalize_log1p` (writes `RNA_log1p`) —
i.e. **Seurat → native R writer → .cytome → PIASO (Python)** works end to end. Round-trip and
cross-language conformance (matrix + obs + embeddings, all three codecs) are covered by the package
tests.
