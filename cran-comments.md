## Test environments

* local: Linux, R 4.3.3 (conda), `R CMD check --as-cran`
* GitHub Actions: ubuntu-latest (release, oldrel-1), macos-latest —
  `.github/workflows/R-CMD-check.yaml`

## R CMD check results

Tests pass and the package documents and installs cleanly. The remaining
findings on the local run are **all** either missing external tools on that
machine or a consequence of the repository not being published yet, not
package defects. Listing them so the next person does not re-diagnose them:

| finding | cause |
|---|---|
| `checking top-level files ... WARNING` | `pandoc` and `checkbashisms` are not installed locally, so `README.md` / `NEWS.md` cannot be checked |
| `checking PDF version of manual ... WARNING` / `... without index ... ERROR` | no LaTeX (`pdflatex`) on the machine |
| `checking HTML version of manual ... NOTE` | no `tidy` |
| `checking compilation flags used ... NOTE` — `-march=nocona` | comes from the **conda** R's `Makeconf`, not from this package's `Makevars`; absent on a standard R build |
| `checking for future file timestamps ... NOTE` | clock skew on the build machine |
| `checking CRAN incoming feasibility ... NOTE` — 404 on the GitHub URLs | the repository is not created yet; resolves on publication |
| `checking for non-standard things in the check directory ... NOTE` | `cytome-manual.tex`, a by-product of the failed PDF step above |

Re-run on a machine with pandoc + LaTeX before submitting, and paste the real
tail here.

## Notes for the reviewer

* **New submission.**
* `SystemRequirements: liblz4, libzstd, zlib`. These are the compression
  codecs a `.cytome` file is written with, so they are required rather than
  optional. `configure` locates them via `pkg-config`, accepts `CYTOME_CFLAGS`
  / `CYTOME_LIBS` as an override, and fails at configure time with
  distribution-specific install instructions rather than part-way through
  compilation. Vendoring the two libraries is the next step if CRAN would
  prefer no system dependency.
* The `DelayedArray` / `S4Arrays` / `SparseArray` stack is in `Suggests`. The
  array-seed class and its methods are therefore registered in `.onLoad()`
  only when those packages are present, rather than declared at build time.
