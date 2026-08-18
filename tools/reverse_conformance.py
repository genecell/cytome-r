#!/usr/bin/env python3
"""Read a cytome WRITTEN BY R and check it against the same expectations.

The existing conformance test runs one way: Python writes `reference.cytome`,
R reads it, R compares to `expected_*.csv`. Nothing checked the other
direction -- which is why `write_cytome.Seurat()` could drop every assay but
the first and no test noticed.

    Rscript tools/write_conformance.R  /tmp/from_r.cytome
    python  tools/reverse_conformance.py /tmp/from_r.cytome inst/extdata

Exits non-zero on the first mismatch.
"""
import sys
import pathlib
import numpy as np
import pandas as pd
import cytome


def main(argv):
    if len(argv) != 3:
        print(__doc__)
        return 2
    path, extdata = argv[1], pathlib.Path(argv[2])
    ds = cytome.open(path)
    problems = []

    mats = sorted(r[0] for r in
                  ds._conn.execute("SELECT matrix_name FROM matrix_meta"))
    print(f"matrices: {mats}")
    for want in ("RNA_counts", "ATAC_counts"):
        if want not in mats:
            problems.append(f"{want} missing -- a modality was dropped on write")

    for name, csv, modality in (("RNA_counts", "expected_rna.csv", "RNA"),
                                ("ATAC_counts", "expected_atac.csv", "ATAC")):
        if name not in mats:
            continue
        exp = pd.read_csv(extdata / csv, index_col=0)      # cells x features
        got = ds.to_anndata(modality=modality)
        arr = got.X.toarray() if hasattr(got.X, "toarray") else np.asarray(got.X)
        if arr.shape != exp.shape:
            problems.append(f"{name}: shape {arr.shape} != expected {exp.shape}")
            continue
        if not np.allclose(arr, exp.values, atol=1e-6):
            worst = np.abs(arr - exp.values).max()
            problems.append(f"{name}: values differ, max |diff| = {worst:.3g}")
        else:
            print(f"{name}: {arr.shape} matches {csv}")

    ds.close()
    if problems:
        print("\nREVERSE CONFORMANCE FAILED")
        for p in problems:
            print("  -", p)
        return 1
    print("\nreverse conformance ok: a cytome written by R reads in Python")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
