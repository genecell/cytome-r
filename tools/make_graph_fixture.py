"""Write inst/extdata/graph_chunks.cytome with the Python package (cytome >= 0.3.6),
whose graphs are stored as compressed row chunks, and the entries the R reader must
return. Run from the package root:  python tools/make_graph_fixture.py"""
import os, sqlite3
import numpy as np, pandas as pd, scipy.sparse as sp
import cytome
from cytome.core.graph import GraphStore

OUT = "inst/extdata/graph_chunks.cytome"
for sfx in ("", "-wal", "-shm"):
    if os.path.exists(OUT + sfx):
        os.remove(OUT + sfx)
GraphStore._CHUNK_NNZ = 150          # several chunks from a small graph
rng = np.random.default_rng(0)
n = 60
ds = cytome.create(OUT)
ds.set_entity("cells", pd.DataFrame({"barcode": [f"c{i}" for i in range(n)]}))
ds.set_entity("genes", pd.DataFrame({"gene_id": [f"G{i}" for i in range(10)]}))
ds.add_matrix("RNA_counts", sp.random(n, 10, density=0.3, format="csr", dtype=np.float32, random_state=rng))
G = sp.random(n, n, density=0.15, format="csr", dtype=np.float64, random_state=rng)
G = (G + G.T).tocsr()
G = G.tolil(); G[n - 1, :] = 0; G[:, n - 1] = 0; G = G.tocsr(); G.eliminate_zeros()   # a last cell without edges
ds.add_graph("connectivities", G.astype(np.float32))
ds.add_graph("distances", G)
ds.flush(); ds.close()
con = sqlite3.connect(OUT)
print("chunks per graph:", con.execute("SELECT graph_name, n_chunks, dtype FROM graph_meta").fetchall())
con.execute("PRAGMA wal_checkpoint(TRUNCATE)"); con.execute("PRAGMA journal_mode=DELETE"); con.execute("VACUUM"); con.close()
for name, M in (("connectivities", G.astype(np.float32)), ("distances", G)):
    t = M.tocoo()
    # hex floats: R parses them exactly on every platform; a 17-digit decimal
    # can land one ulp off where R has no extended precision (arm64 macOS, Windows)
    pd.DataFrame({"i": t.row + 1, "j": t.col + 1, "x": [float(v).hex() for v in t.data.astype(np.float64)]}).to_csv(
        f"inst/extdata/expected_graph_{name}.csv", index=False)
print("cytome", cytome.__version__, "->", OUT, os.path.getsize(OUT), "bytes")
