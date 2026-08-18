#!/usr/bin/env Rscript
## One command to take the package from "works on my machine" to "checked".
##
## Needs an R with a WORKING C/C++ TOOLCHAIN plus roxygen2 -- the conda env the
## package was developed in has neither (source installs of xfun, yaml,
## commonmark and xml2 all fail to compile there), which is why man/ has never
## been generated and R CMD check has never run.
##
##     Rscript tools/prepare_release.R
##
stopifnot(file.exists("DESCRIPTION"))

need <- c("roxygen2", "Rcpp", "testthat")
missing <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing))
  stop("install first: ", paste(missing, collapse = ", "))

## The Rcpp entry points are still named cytomer_* from before the package
## was renamed. Renaming them requires recompiling, so it happens here rather
## than in a commit that would leave RcppExports.R calling symbols the
## installed .so does not have.
message("== 0. rename the cytomer_* C++ entry points")
for (f in c(Sys.glob("R/*.R"), Sys.glob("src/*.cpp"))) {
  if (!file.exists(f)) next
  txt <- readLines(f, warn = FALSE)
  new_txt <- gsub("cytomer_", "cytome_", txt, fixed = TRUE)
  if (!identical(txt, new_txt)) { writeLines(new_txt, f); message("  renamed in ", f) }
}
unlink(c("src/cytomer.so", "src/cytome.so", Sys.glob("src/*.o")))

message("== 1. Rcpp::compileAttributes()")
Rcpp::compileAttributes(".")

message("== 2. roxygen2::roxygenise()  -> man/ and NAMESPACE")
roxygen2::roxygenise(".", roclets = c("rd", "namespace"))

message("== 3. R CMD build")
tarball <- system2("R", c("CMD", "build", "."), stdout = TRUE)
writeLines(tarball)
tgz <- list.files(".", pattern = "^cytome_.*\\.tar\\.gz$", full.names = TRUE)
tgz <- tgz[order(file.info(tgz)$mtime, decreasing = TRUE)][1]
message("built: ", tgz)

message("== 4. R CMD check --as-cran")
status <- system2("R", c("CMD", "check", "--as-cran", shQuote(tgz)))
message(if (status == 0) "check finished" else "check FAILED")
message("read cytome.Rcheck/00check.log, then paste its tail into cran-comments.md")
quit(status = status)
