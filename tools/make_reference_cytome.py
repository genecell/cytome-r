"""Generate the cross-language conformance fixture: a small deterministic .cytome written by
Python cytome, plus the expected values as CSVs, into inst/extdata/. The R conformance test reads
the cytome with `cytomer` and asserts it matches these expectations — proving the format is
language-agnostic (Python writes, R reads, bit-for-bit).

Run once (in a python env with cytome):  python tools/make_reference_cytome.py
"""
import os, shutil, numpy as np, pandas as pd, scipy.sparse as sp
import cytome

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "inst", "extdata")
os.makedirs(OUT, exist_ok=True)
CYT = os.path.join(OUT, "reference.cytome")
if os.path.exists(CYT):
    os.remove(CYT)

rng = np.random.RandomState(0)
n_cells, n_genes, n_peaks = 12, 8, 5

# deterministic dense matrices (cells x features), some zeros
rna = (rng.poisson(1.5, size=(n_cells, n_genes)) * (rng.rand(n_cells, n_genes) > 0.3)).astype(np.float32)
atac = (rng.poisson(1.0, size=(n_cells, n_peaks)) * (rng.rand(n_cells, n_peaks) > 0.4)).astype(np.float32)
umap = rng.randn(n_cells, 2).astype(np.float32)

cells = pd.DataFrame({
    "barcode": [f"BC{i:02d}-1" for i in range(n_cells)],
    "sample_id": ["s1"] * (n_cells // 2) + ["s2"] * (n_cells - n_cells // 2),
    "n_counts": rna.sum(1).astype(int),
})
genes = pd.DataFrame({"gene_id": [f"ENSG{i}" for i in range(n_genes)],
                      "symbol": [f"Gene{i}" for i in range(n_genes)]})
peaks = pd.DataFrame({"peak_id": [f"chr1:{1000*i}-{1000*i+200}" for i in range(n_peaks)],
                      "chr": "chr1",
                      "start": [1000 * i for i in range(n_peaks)],
                      "end_": [1000 * i + 200 for i in range(n_peaks)]})

ds = cytome.create(CYT)
ds.set_entity("cells", cells)
ds.set_entity("genes", genes)
ds.set_entity("peaks", peaks)

# Write RNA with zstd (add_matrix default) and ATAC with lz4 (via the chunked writer) so the
# R conformance test exercises BOTH codecs.
ds.add_matrix("RNA_counts", sp.csr_matrix(rna))            # zstd
w = ds.create_layer_writer("ATAC_counts", n_rows=n_cells, n_cols=n_peaks, dtype="float32",
                           compression="lz4", col_entity="peaks")
w.write_chunk(sp.csr_matrix(atac), 0)
w.finalize()                                               # lz4
try:
    ds.add_embedding("X_umap", umap)
except Exception as e:
    print("embedding add failed (non-fatal):", e)
ds.flush()
# report compressions actually used
import sqlite3
rows = ds._conn.execute("SELECT matrix_name, compression FROM matrix_chunks GROUP BY matrix_name, compression").fetchall()
print("matrix compressions:", rows)
ds.close()

# expected CSVs (R reads these and compares)
pd.DataFrame(rna, index=cells.barcode, columns=genes.gene_id).to_csv(f"{OUT}/expected_rna.csv")
pd.DataFrame(atac, index=cells.barcode, columns=peaks.peak_id).to_csv(f"{OUT}/expected_atac.csv")
pd.DataFrame(umap, index=cells.barcode, columns=["umap1", "umap2"]).to_csv(f"{OUT}/expected_umap.csv")
cells.to_csv(f"{OUT}/expected_cells.csv", index=False)
genes.to_csv(f"{OUT}/expected_genes.csv", index=False)
print(f"wrote {CYT} + expected_*.csv  ({n_cells} cells, RNA {n_genes}g, ATAC {n_peaks}p)")
