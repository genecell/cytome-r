"""Write the stores the R tests read, with the writers that make them in practice, and what
R must return for each. Run from the package root:

    python tools/make_python_fixtures.py [OUT]       # OUT defaults to tests/testthat/fixtures

Each store is written gzip-compressed (<name>.cytome.gz; the tests unpack it into a temporary
folder) next to <name>.json, the values R must read back. A Python store carries every table of
the schema, so an unpacked fixture is close to a megabyte, but compressed it is a few tens of
kilobytes.

    py_types       from_anndata: categoricals with string, integer, float and boolean categories,
                   booleans with and without missing values, nullable integers, strings with
                   non-ASCII text, a .raw wider than X, a sparse obsm of 600 columns (stored as a
                   matrix wider than the genes), int32 and float32 layers, float32/float64/int64
                   obsm, obs and var graphs, a var embedding, several row chunks
    py_10x         from_cellranger of a Cell Ranger multiome .h5: int32 RNA and ATAC counts, a
                   duplicated gene symbol, plus an int64 spatial embedding
    py_frag_rows   fragments in the per-row tables (the deprecated import_fragments)
    py_merge_mixed merge of an int32 store and a float32 store: the chunks differ in dtype
    piaso_atac     PIASO's Rust importer and quantifier (set PIASO_BIN to the folder holding
                   cytome-import-fragments and piaso-atac): column-major tiles and peak matrices,
                   fragments in compressed chunks (encoding 1); skipped when PIASO_BIN is unset

Expected values are the values written (for py_types and py_merge_mixed, where the Python reader
of today misreads some of them) or, for piaso_atac, Python's own reading of the store. Floats are
given as hex strings (float.hex), which R parses exactly on every platform.
"""
import gzip
import json
import os
import shutil
import sqlite3
import subprocess
import sys
import tempfile
from contextlib import closing
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

warnings.simplefilter("ignore")
import anndata as ad  # noqa: E402
import cytome  # noqa: E402
import cytome.core.dataset as dsmod  # noqa: E402
from cytome.core.graph import GraphStore  # noqa: E402
from cytome.core.measurement import MeasurementLayer  # noqa: E402

OUT = Path(sys.argv[1] if len(sys.argv) > 1 else "tests/testthat/fixtures")
OUT.mkdir(parents=True, exist_ok=True)
WORK = Path(tempfile.mkdtemp(prefix="cytome_fixtures_", dir=os.environ.get("TMPDIR")))


def hexs(a):
    return [float(v).hex() for v in np.asarray(a, dtype=np.float64).ravel()]


def triplets(M):
    t = sp.coo_matrix(M)
    t.sum_duplicates()
    o = np.lexsort((t.col, t.row))
    return {"shape": [int(M.shape[0]), int(M.shape[1])], "i": t.row[o].tolist(),
            "j": t.col[o].tolist(), "x": hexs(t.data[o])}


def dense(E):
    E = np.asarray(E)
    return {"shape": [int(E.shape[0]), int(E.shape[1])], "x": hexs(E)}  # row-major


def remove_store(p):
    for sfx in ("", "-wal", "-shm", "-journal"):
        if os.path.exists(str(p) + sfx):
            os.remove(str(p) + sfx)


def finish(name, store, expected, journal="DELETE"):
    """Checkpoint, optionally leave WAL mode, compress into OUT, write the expectations."""
    c = sqlite3.connect(str(store))
    c.execute("PRAGMA wal_checkpoint(TRUNCATE)")
    if journal:                     # a rollback-journal file, compacted
        c.execute(f"PRAGMA journal_mode={journal}")
        c.execute("VACUUM")
    c.close()
    assert not os.path.exists(str(store) + "-wal") or os.path.getsize(str(store) + "-wal") == 0
    with open(store, "rb") as fh, gzip.open(OUT / f"{name}.cytome.gz", "wb", compresslevel=9) as gz:
        shutil.copyfileobj(fh, gz)
    expected = {"writer": f"cytome {cytome.__version__}", **expected}
    (OUT / f"{name}.json").write_text(json.dumps(expected, indent=None, sort_keys=True))
    print(f"{name}: {os.path.getsize(store):,} bytes -> "
          f"{os.path.getsize(OUT / f'{name}.cytome.gz'):,} compressed")


def obs_expectation(df):
    """What R must return for each column: the type and the values, missing as None.

    A categorical's categories and values keep their own JSON type (string, number or
    boolean); R's levels are their as.character()."""
    out = {}
    for col in df.columns:
        s = df[col]
        if isinstance(s.dtype, pd.CategoricalDtype):
            cats = s.cat.categories
            kind = ("logical" if pd.api.types.is_bool_dtype(cats.dtype) else
                    "numeric" if pd.api.types.is_numeric_dtype(cats.dtype) else "character")
            conv = {"logical": bool, "numeric": float, "character": str}[kind]
            out[col] = {"r": "factor", "ordered": bool(s.cat.ordered), "kind": kind,
                        "categories": [conv(c) for c in cats.tolist()],
                        "values": [None if pd.isna(v) else conv(v) for v in s.astype(object)]}
        elif pd.api.types.is_bool_dtype(s.dtype) or str(s.dtype) == "boolean":
            out[col] = {"r": "logical", "values": [None if pd.isna(v) else bool(v) for v in s.astype(object)]}
        elif pd.api.types.is_integer_dtype(s.dtype):
            out[col] = {"r": "integer", "values": [None if pd.isna(v) else int(v) for v in s.astype(object)]}
        elif pd.api.types.is_float_dtype(s.dtype):
            out[col] = {"r": "double", "values": [None if pd.isna(v) else float(v).hex() for v in s.astype(object)]}
        else:
            out[col] = {"r": "character", "values": [None if v is None or (isinstance(v, float) and np.isnan(v))
                                                     else str(v) for v in s.astype(object)]}
    return out


# ---------------------------------------------------------------------------------------------
# py_types: from_anndata with every column type and the shapes that broke the R reader
# ---------------------------------------------------------------------------------------------
def make_types():
    rng = np.random.default_rng(0)
    n, m = 23, 7
    X = rng.poisson(1.0, size=(n, m)).astype(np.float32)
    X[0:4] = 0                                     # an all-empty first chunk
    X[7] = 0
    obs = pd.DataFrame({
        "cluster": pd.Categorical(rng.choice(["b", "a", "c"], n), categories=["b", "a", "c", "unused"]),
        "stage": pd.Categorical(rng.choice(["early", "mid", "late"], n),
                                categories=["early", "mid", "late"], ordered=True),
        "int_cat": pd.Categorical(rng.choice([1, 2, 3], n), categories=[1, 2, 3, 4]),
        "float_cat": pd.Categorical(rng.choice([0.5, 1.0, 2.25], n)),
        "bool_cat": pd.Categorical(rng.choice([True, False], n)),
        "zelltyp": pd.Categorical(rng.choice(["Neuron ä", "Glia ö"], n)),
        "flag": rng.random(n) > 0.5,
        "flag_na": pd.array([None if i % 5 == 0 else bool(i % 2) for i in range(n)], dtype="boolean"),
        "n_genes": rng.integers(0, 100, n).astype(np.int64),
        "n_umi_na": pd.array([None if i % 4 == 0 else 10 * i for i in range(n)], dtype="Int64"),
        "score": np.where(rng.random(n) > 0.7, np.nan, rng.random(n)),
        "ratio32": rng.random(n).astype(np.float32),
        "label": [None if i % 6 == 0 else f"cell µ {i}" for i in range(n)],
        "col/with slash": np.arange(n, dtype=np.float64),
    }, index=[f"AAAC-{i}" for i in range(n)])
    obs.loc[obs.index[3], "cluster"] = np.nan
    genes = ["G1", "Gene_A", "G3", "G1", "MT-CO1", "G6", "G7"]          # a duplicated symbol
    var = pd.DataFrame({"gene_ids": [f"ENSG{i:05d}" for i in range(m)],
                        "highly_variable": [True, False] * 3 + [True],
                        "means": rng.random(m)}, index=genes)
    a = ad.AnnData(X=sp.csr_matrix(X), obs=obs, var=var)
    a.layers["log1p"] = sp.csr_matrix(np.log1p(X).astype(np.float32))
    a.layers["counts_int"] = sp.csr_matrix(X.astype(np.int32))
    raw_X = np.hstack([X, rng.poisson(2.0, size=(n, 3)).astype(np.float32)])   # .raw wider than X
    a.raw = ad.AnnData(X=sp.csr_matrix(raw_X), obs=obs[[]],
                       var=pd.DataFrame(index=genes[:4] + ["G1b", "G6", "G7", "R1", "R2", "R3"]))
    a.obsm["X_umap"] = rng.random((n, 2)).astype(np.float32)
    a.obsm["X_pca"] = rng.normal(size=(n, 5))
    a.obsm["spatial"] = rng.integers(0, 5000, (n, 2)).astype(np.int64)
    a.obsm["X_wide"] = sp.random(n, 600, density=0.05, random_state=4, format="csr", dtype=np.float32)
    C = sp.random(n, n, density=0.2, random_state=1, format="csr", dtype=np.float32)
    a.obsp["connectivities"] = C
    a.obsp["distances"] = sp.csr_matrix(C.astype(np.float64) * 3.0)
    a.varm["PCs"] = rng.normal(size=(m, 3))
    a.varp["corr"] = sp.random(m, m, density=0.3, random_state=2, format="csr", dtype=np.float64)

    dsmod.compute_chunk_size = lambda *a_, **k_: 4          # several row chunks
    GraphStore._CHUNK_NNZ = 9                               # several graph chunks
    store = WORK / "py_types.cytome"
    remove_store(store)
    cytome.from_anndata(a, output=str(store), force=True).close()
    # an all-zero layer: chunks whose blobs are empty
    ds = cytome.open(str(store), mode="r+")
    w = ds.create_layer_writer("RNA_zero", n_rows=n, n_cols=m, dtype=np.float64, compression="lz4",
                               col_entity="genes", storage_chunk_size=8)
    w.write_chunk(sp.csr_matrix((n, m), dtype=np.float64), 0)
    w.finalize()
    # a row-major matrix that also carries column chunks (a feature index)
    MeasurementLayer(ds._conn, "RNA_log1p").build_feature_index(chunk_size=3, compression="zlib")
    ds.close()

    c = sqlite3.connect(str(store))
    mats = {r[0]: r for r in c.execute("SELECT matrix_name, n_rows, n_cols, col_entity FROM matrix_meta")}
    embs = {r[0]: r for r in c.execute("SELECT array_name, entity, dtype FROM embedding_meta")}
    c.close()
    print("py_types matrices", mats, "embeddings", embs)
    names = {"RNA_counts": X, "RNA_log1p": np.log1p(X), "RNA_counts_int": X,
             "RNA_zero": np.zeros((n, m))}
    raw_name = [k for k in mats if "raw" in k]
    wide_name = [k for k in mats if "wide" in k]
    assert len(raw_name) == 1 and len(wide_name) == 1, (raw_name, wide_name)
    names[raw_name[0]] = raw_X
    names[wide_name[0]] = a.obsm["X_wide"].toarray()
    expected = {
        "n_cells": n, "barcodes": list(obs.index),
        "matrices": {k: triplets(sp.csr_matrix(v)) for k, v in names.items()},
        "embeddings": {}, "graphs": {},
        "obs": obs_expectation(obs),
        "var": {"gene_id": {"r": "character", "values": list(var["gene_ids"])},
                "symbol": {"r": "character", "values": genes},
                "highly_variable": {"r": "logical", "values": [bool(v) for v in var["highly_variable"]]},
                "means": {"r": "double", "values": hexs(var["means"])}},
        "raw_matrix": raw_name[0], "wide_matrix": wide_name[0],
    }
    for key, arr in (("X_umap", a.obsm["X_umap"]), ("X_pca", a.obsm["X_pca"]),
                     ("spatial", a.obsm["spatial"])):
        hit = [e for e in embs if e.endswith(key.replace("X_", "")) and embs[e][1] == "cells"]
        assert len(hit) == 1, (key, hit)
        expected["embeddings"][hit[0]] = {**dense(arr), "dtype": embs[hit[0]][2]}
    varm = [e for e in embs if embs[e][1] != "cells"]
    expected["var_embeddings"] = varm
    with closing(sqlite3.connect(str(store))) as c:
        for g, axis in c.execute("SELECT graph_name, axis FROM graph_meta"):
            truth = C if g.endswith("connectivities") else (
                sp.csr_matrix(C.astype(np.float64) * 3.0) if g.endswith("distances") else a.varp["corr"])
            expected["graphs"][g] = {**triplets(truth), "axis": axis}
    finish("py_types", store, expected)


# ---------------------------------------------------------------------------------------------
# py_10x: Cell Ranger multiome .h5 through from_cellranger -> int32 matrices
# ---------------------------------------------------------------------------------------------
def make_10x():
    import h5py
    rng = np.random.default_rng(3)
    n_cells = 9
    genes = ["ENSG1", "ENSG2", "ENSG3", "ENSG4"]
    peaks = ["chr1:100-200", "chr1:300-400", "chr2:10-90"]
    ids = genes + peaks
    M = rng.poisson(1.5, size=(len(ids), n_cells)).astype(np.int32)     # features x cells
    M[:, 0] = 0
    Mc = sp.csc_matrix(M)
    d = WORK / "tenx_run"
    d.mkdir(exist_ok=True)
    bcs = [f"BC{i}-1" for i in range(n_cells)]
    with h5py.File(d / "filtered_feature_bc_matrix.h5", "w") as f:
        g = f.create_group("matrix")
        g["barcodes"] = np.array([b.encode() for b in bcs])
        g["data"] = Mc.data.astype(np.int32)
        g["indices"] = Mc.indices.astype(np.int64)
        g["indptr"] = Mc.indptr.astype(np.int64)
        g["shape"] = np.array(Mc.shape, dtype=np.int32)
        ft = g.create_group("features")
        ft["id"] = np.array([x.encode() for x in ids])
        ft["name"] = np.array([x.encode() for x in ["Sox2", "Pax6", "Sox2", "Gfap"] + peaks])
        ft["feature_type"] = np.array([b"Gene Expression"] * 4 + [b"Peaks"] * 3)
        ft["genome"] = np.array([b"GRCh38"] * 7)
    store = WORK / "py_10x.cytome"
    remove_store(store)
    ds = cytome.from_cellranger(str(d), str(store), force=True, import_fragments=False)
    spatial = rng.integers(0, 20000, (n_cells, 2)).astype(np.int64)
    ds.add_embedding("RNA_spatial", spatial)
    ds.close()
    with closing(sqlite3.connect(str(store))) as c:
        print("py_10x", c.execute("SELECT matrix_name, dtype, n_chunks, has_csc FROM matrix_meta").fetchall(),
              c.execute("SELECT array_name, dtype FROM embedding_meta").fetchall())
    expected = {
        "n_cells": n_cells, "barcodes": bcs,
        "matrices": {"RNA_counts": triplets(sp.csr_matrix(M[:4].T)),
                     "ATAC_counts": triplets(sp.csr_matrix(M[4:].T))},
        "embeddings": {"RNA_spatial": {**dense(spatial), "dtype": "int64"}},
        "var": {"gene_id": {"r": "character", "values": genes},
                "symbol": {"r": "character", "values": ["Sox2", "Pax6", "Sox2", "Gfap"]}},
        "peaks": peaks,
    }
    finish("py_10x", store, expected)


# ---------------------------------------------------------------------------------------------
# py_frag_rows: fragments in the per-row tables (fragments_<chrom> with an R*-tree)
# ---------------------------------------------------------------------------------------------
def write_fragments(path, rows):
    with gzip.open(path, "wt") as fh:
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")


def fragment_rows(bcs, rng, n=60):
    rows = []
    for chrom in ("chr1", "chr2"):
        starts = np.sort(rng.integers(100, 9000, n))
        for s in starts:
            rows.append((chrom, int(s), int(s) + int(rng.integers(30, 400)),
                         bcs[int(rng.integers(0, len(bcs)))], int(rng.integers(1, 4))))
    return rows


def make_frag_rows():
    rng = np.random.default_rng(7)
    n = 12
    bcs = [f"CELL{i:02d}-1" for i in range(n)]
    peaks = [f"chr1:{1000 * i}-{1000 * i + 300}" for i in range(1, 6)]
    a = ad.AnnData(X=sp.csr_matrix(rng.poisson(1.0, size=(n, len(peaks))).astype(np.float32)),
                   obs=pd.DataFrame(index=bcs), var=pd.DataFrame(index=peaks))
    store = WORK / "py_frag_rows.cytome"
    remove_store(store)
    ds = cytome.from_anndata(a, modality="ATAC", output=str(store), force=True)
    rows = fragment_rows(bcs, rng, 50)
    write_fragments(WORK / "frags_rows.tsv.gz", rows)
    ds.ATAC.import_fragments(str(WORK / "frags_rows.tsv.gz"))
    ds.close()
    with closing(sqlite3.connect(str(store))) as c:
        print("py_frag_rows", c.execute("SELECT COUNT(*) FROM fragments_chr1").fetchone(),
              c.execute("SELECT COUNT(*) FROM fragment_chunks").fetchone(),
              c.execute("SELECT COUNT(*) FROM fragments_chr1_rtree").fetchone())
    expected = {"n_cells": n, "barcodes": bcs, "fragments": {}}
    for chrom in ("chr1", "chr2"):
        rr = sorted([r for r in rows if r[0] == chrom], key=lambda r: (r[1], r[2], bcs.index(r[3])))
        expected["fragments"][chrom] = {"start": [r[1] for r in rr], "end": [r[2] for r in rr],
                                        "cell_idx": [bcs.index(r[3]) for r in rr],
                                        "dup_count": [r[4] for r in rr]}
    finish("py_frag_rows", store, expected)


# ---------------------------------------------------------------------------------------------
# py_merge_mixed: merge an int32 store with a float32 store
# ---------------------------------------------------------------------------------------------
def make_merge_mixed():
    rng = np.random.default_rng(11)
    m = 5
    genes = [f"g{i}" for i in range(m)]
    A = rng.poisson(2.0, size=(6, m)).astype(np.int32)
    B = rng.poisson(2.0, size=(7, m)).astype(np.float32)       # integral, so stored as counts
    pa, pb = WORK / "mix_a.cytome", WORK / "mix_b.cytome"
    for p, M, pre in ((pa, A, "a"), (pb, B, "b")):
        remove_store(p)
        aa = ad.AnnData(X=sp.csr_matrix(M), obs=pd.DataFrame(index=[f"{pre}{i}" for i in range(M.shape[0])]),
                        var=pd.DataFrame(index=genes))
        cytome.from_anndata(aa, output=str(p), force=True).close()
    store = WORK / "py_merge_mixed.cytome"
    remove_store(store)
    cytome.merge([str(pa), str(pb)], str(store), force=True).close()
    with closing(sqlite3.connect(str(store))) as c:
        meta = c.execute("SELECT matrix_name, dtype FROM matrix_meta").fetchall()
        chunks = c.execute("SELECT matrix_name, dtype, COUNT(*) FROM matrix_chunks GROUP BY 1, 2").fetchall()
        bcs = [r[0] for r in c.execute("SELECT barcode FROM cells ORDER BY cell_idx")]
    print("py_merge_mixed", meta, chunks)
    name = meta[0][0]
    truth = np.vstack([A.astype(np.float64), B.astype(np.float64)])
    expected = {"n_cells": truth.shape[0], "barcodes": bcs, "chunk_dtypes": sorted({r[1] for r in chunks}),
                "matrices": {name: triplets(sp.csr_matrix(truth))}}
    finish("py_merge_mixed", store, expected)


# ---------------------------------------------------------------------------------------------
# piaso_atac: PIASO's Rust importer and quantifier (column-major matrices, fragment chunks)
# ---------------------------------------------------------------------------------------------
def bgzip(src, dst):
    try:
        import pysam
        pysam.tabix_compress(str(src), str(dst), force=True)
    except ImportError:
        with open(dst, "wb") as out:
            subprocess.run([os.environ.get("BGZIP", "bgzip"), "-c", str(src)], stdout=out, check=True)


def make_piaso_atac():
    bindir = os.environ.get("PIASO_BIN")
    if not bindir:
        print("piaso_atac: PIASO_BIN not set, skipped")
        return
    rng = np.random.default_rng(5)
    n_cells = 30
    bcs = [f"AAACGA{i:04d}-1" for i in range(n_cells)]
    rows = []
    for chrom, length in (("chr1", 20000), ("chr2", 12000)):
        for s in np.sort(rng.integers(100, length - 600, 400)):
            rows.append((chrom, int(s), int(s) + int(rng.integers(40, 500)),
                         bcs[int(rng.integers(0, n_cells))], int(rng.integers(1, 4))))
    d = WORK / "piaso_run"
    d.mkdir(exist_ok=True)
    with open(d / "frags.tsv", "w") as fh:
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    bgzip(d / "frags.tsv", d / "frags.tsv.gz")
    (d / "bc.txt").write_text("\n".join(bcs) + "\n")
    (d / "genome.sizes").write_text("chr1\t20000\nchr2\t12000\n")
    peaks = [("chr1", 500, 1500), ("chr1", 3000, 3600), ("chr1", 8000, 9000), ("chr1", 15000, 16500),
             ("chr2", 200, 900), ("chr2", 4000, 5200), ("chr2", 9000, 9400)]
    (d / "peaks.bed").write_text("".join(f"{c}\t{s}\t{e}\n" for c, s, e in peaks))
    store = WORK / "piaso_atac.cytome"
    remove_store(store)
    subprocess.run([os.path.join(bindir, "cytome-import-fragments"), "--fragments", str(d / "frags.tsv.gz"),
                    "--barcodes", str(d / "bc.txt"), "--output", str(store), "--chunk-size", "60",
                    "--threads", "2", "--finalize-threads", "2", "--genome", str(d / "genome.sizes"),
                    "--tile-size", "2000", "--memory-budget-mb", "64"], check=True, capture_output=True)
    subprocess.run([os.path.join(bindir, "piaso-atac"), "quantify", "--cytome", str(store), "--peaks",
                    str(d / "peaks.bed"), "--threads", "2", "--memory-budget-mb", "64"],
                   check=True, capture_output=True)
    ds = cytome.open(str(store), mode="r")
    expected = {"n_cells": n_cells, "barcodes": bcs, "matrices": {}, "layouts": {}, "fragments": {}}
    for (name, n_chunks, has_csc) in ds._conn.execute(
            "SELECT matrix_name, n_chunks, has_csc FROM matrix_meta").fetchall():
        expected["matrices"][name] = triplets(MeasurementLayer(ds._conn, name).to_memory())
        expected["layouts"][name] = "csc" if (n_chunks == 0 and has_csc) else "csr"
    fs = ds.ATAC.fragments
    for chrom in ("chr1", "chr2"):
        s, e, cc = [np.concatenate(x) for x in zip(*fs.iter_chromosome_chunks(chrom))]
        o = np.lexsort((cc, e, s))
        expected["fragments"][chrom] = {"start": s[o].tolist(), "end": e[o].tolist(),
                                        "cell_idx": cc[o].tolist()}
    expected["peaks"] = [f"{c}:{s}-{e}" for c, s, e in peaks]
    ds.close()
    # The store as PIASO leaves it: WAL mode and 64 KB pages.
    finish("piaso_atac", store, expected, journal=None)


if __name__ == "__main__":
    make_types()
    make_10x()
    make_frag_rows()
    make_merge_mixed()
    make_piaso_atac()
    shutil.rmtree(WORK, ignore_errors=True)
