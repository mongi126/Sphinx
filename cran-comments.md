## Test environments
* local, Windows 11, R 4.4.2
* GitHub Actions (Windows, macOS, Linux)

## R CMD check results
0 errors | 0 warnings | 0 notes

## Downstream dependencies
There are no downstream dependencies for this package at the time of submission.

## Comments for CRAN team
This is the first submission of Sphinx.

* Examples use synthetic data helpers defined in `R/example-utils.R`.
* Examples requiring external data files or optional Bioconductor packages
  are wrapped in `\dontrun{}`.
* Vignettes are built from `vignettes/` during `R CMD build`.
