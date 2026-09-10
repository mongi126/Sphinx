# Installation

## Requirements

  - R \>= 4.3
  - For enrichment: `enrichR` (suggested); optionally `clusterProfiler`
    / `org.Hs.eg.db`
  - Composition heatmaps use Bioconductor `ComplexHeatmap` / `circlize`
    (Imports)

## Install from GitHub

``` r
install.packages("devtools")
devtools::install_github("mongi126/Sphinx")
```

## Install from source tarball

``` r
install.packages("Sphinx_1.0.1.tar.gz", repos = NULL, type = "source")
```

## Load and check

``` r
library(Sphinx)
packageVersion("Sphinx")
help(package = "Sphinx")
```

## Example data in the docs

Tutorials use **reg055\_A** from Schurch *et al.* CODEX imaging
(`schurch2020coordinated`):

  - Technology: CODEX
  - \~3,887 cells, 58 protein channels
  - Coordinates in `obsm/spatial`
  - Published labels in `obs/ClusterName`

Any table with `Cell_ID`, `X`, `Y`, and a categorical annotation column
works the same way.
