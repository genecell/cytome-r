"""Write inst/extdata/zlib_as_zstd.cytome as the Python package writes it when
the optional zstandard package is missing: every blob is zlib, labelled "zstd".
Run from the package root, in an environment with cytome but WITHOUT zstandard:
    python tools/make_zlib_fixture.py"""
import os, sqlite3
import numpy as np, pandas as pd, scipy.sparse as sp
import cytome
from cytome.io import compression

assert not compression._ZSTD_AVAILABLE, "run this where zstandard is not installed"
OUT = "inst/extdata/zlib_as_zstd.cytome"
for sfx in ("", "-wal", "-shm"):
    if os.path.exists(OUT + sfx):
        os.remove(OUT + sfx)
rng = np.random.default_rng(0)          # the store the hub's CI writes
ds = cytome.create(OUT)
ds.set_entity("cells", pd.DataFrame({"barcode": [f"c{i}" for i in range(120)], "cell_type": rng.choice(["A", "B"], 120)}))
ds.set_entity("genes", pd.DataFrame({"gene_id": [f"G{i}" for i in range(50)]}))
X = sp.random(120, 50, density=0.2, format="csr", dtype=np.float32, random_state=rng)
ds.add_matrix("RNA_counts", X)
ds.add_embedding("RNA_umap", rng.normal(size=(120, 2)).astype(np.float32))
ds.flush(); ds.close()
con = sqlite3.connect(OUT)
labels = con.execute("SELECT DISTINCT compression FROM matrix_chunks").fetchall()
first = con.execute("SELECT data_blob FROM matrix_chunks LIMIT 1").fetchone()[0][:2].hex()
print("labels", labels, "first bytes", first)
assert labels == [("zstd",)] and first.startswith("78")
con.execute("PRAGMA wal_checkpoint(TRUNCATE)"); con.execute("PRAGMA journal_mode=DELETE"); con.execute("VACUUM"); con.close()
t = X.T.tocoo()                         # features x cells, as read_cytome_matrix returns it
pd.DataFrame({"i": t.row + 1, "j": t.col + 1, "x": t.data.astype(np.float64)}).to_csv(
    "inst/extdata/expected_zlib_as_zstd.csv", index=False, float_format="%.9g")
print("cytome", cytome.__version__, "->", OUT, os.path.getsize(OUT), "bytes")
