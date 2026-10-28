"""Read a store written by the R package with the Python package, and compare it with what R
wrote (a JSON the R test writes beside it). Prints one line per check and "ok" at the end;
exits 1 on the first difference.

    python check_r_store.py <store.cytome> <expected.json>
"""
import json
import sys

import numpy as np
import pandas as pd

import cytome


def fail(msg):
    print("MISMATCH:", msg)
    sys.exit(1)


def dense(m):
    return m.toarray() if hasattr(m, "toarray") else np.asarray(m)


def main(path, expected_path):
    e = json.load(open(expected_path))
    ds = cytome.open(path, mode="r")
    a = cytome.to_anndata(ds)
    if list(a.obs_names) != e["cells"]:
        fail("cell names")
    for mod, want in e["matrices"].items():
        ad = a if mod == "RNA" else cytome.to_anndata(ds, modality=mod)
        X = dense(ad.X)
        W = np.asarray(want["values"], dtype=float).reshape(want["shape"])   # cells x features, row-major
        if X.shape != W.shape or not np.allclose(X, W, atol=1e-5):
            fail(f"{mod} X")
        if list(ad.var[want["id_col"]]) != want["features"]:
            fail(f"{mod} features")
        print(f"{mod}: X {X.shape} matches")
    for lyr, want in e.get("layers", {}).items():
        L = dense(a.layers[lyr])
        W = np.asarray(want, dtype=float).reshape(L.shape)
        if not np.allclose(L, W, atol=1e-5):
            fail(f"layer {lyr}")
        print(f"layer {lyr} matches")
    for col, want in e["obs"].items():
        s = a.obs[col]
        kind = want["kind"]
        if kind == "factor":
            if not isinstance(s.dtype, pd.CategoricalDtype):
                fail(f"obs {col} is {s.dtype}, not categorical")
            if list(s.cat.categories) != want["levels"] or bool(s.cat.ordered) != want["ordered"]:
                fail(f"obs {col} categories {list(s.cat.categories)}")
            got = [None if pd.isna(v) else str(v) for v in s.astype(object)]
        elif kind == "logical":
            if s.dtype != bool and str(s.dtype) != "boolean":
                fail(f"obs {col} is {s.dtype}, not boolean")
            got = [None if pd.isna(v) else bool(v) for v in s.astype(object)]
        elif kind == "integer":
            if not pd.api.types.is_integer_dtype(s.dtype):
                fail(f"obs {col} is {s.dtype}, not integer")
            got = [None if pd.isna(v) else int(v) for v in s.astype(object)]
        elif kind == "double":
            got = [None if pd.isna(v) else float(v) for v in s.astype(object)]
            want["values"] = [None if v is None else float(v) for v in want["values"]]
        else:
            got = [None if pd.isna(v) else str(v) for v in s.astype(object)]
        if kind == "double":
            ok = all((g is None and w is None) or (g is not None and w is not None and abs(g - w) < 1e-9)
                     for g, w in zip(got, want["values"]))
        else:
            ok = got == want["values"]
        if not ok:
            fail(f"obs {col}: {got[:5]} != {want['values'][:5]}")
        print(f"obs {col} ({kind}) matches")
    for key, want in e.get("obsm", {}).items():
        if key not in a.obsm:
            fail(f"obsm {key} missing; have {list(a.obsm)}")
        M = np.asarray(a.obsm[key], dtype=float)
        if not np.allclose(M, np.asarray(want, dtype=float).reshape(M.shape), atol=1e-5):
            fail(f"obsm {key}")
        print(f"obsm {key} matches")
    for key, want in e.get("obsp", {}).items():
        if key not in a.obsp:
            fail(f"obsp {key} missing; have {list(a.obsp)}")
        G = dense(a.obsp[key])
        if not np.allclose(G, np.asarray(want, dtype=float).reshape(G.shape)):
            fail(f"obsp {key}")
        print(f"obsp {key} matches")
    from cytome.utils.validation import validate
    report = validate(ds)
    if not report.passed:
        fail(f"validate: {report}")
    print("validate passed")
    ds.close()
    print("ok")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
