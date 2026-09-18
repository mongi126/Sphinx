# Sphinx

[![R-CMD-check](https://github.com/mongi126/Sphinx/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/mongi126/Sphinx/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![GitHub last commit](https://img.shields.io/github/last-commit/mongi126/Sphinx)](https://github.com/mongi126/Sphinx/commits/main)
[![Docs](https://img.shields.io/badge/docs-pkgdown-lightgrey)](https://mongi126.github.io/Sphinx/)

Sphinx is an R toolkit for **spatial proteomics**: preprocessing,
annotation, spatial neighborhood graphs, and functional enrichment—with
publication-ready plots.

![Sphinx workflow](articles/workflow.png)

Documentation demos use **reg055\_A** from Schurch *et al.* CODEX CRC
(2020).

**Example data download:**
[SpatialOmics dataset 60](https://gene.ai.tencent.com/SpatialOmics/dataset?datasetID=60)
→ select **reg055\_A**.

## Installation

``` r
devtools::install_github("mongi126/Sphinx")
# or: install.packages("Sphinx_1.0.1.tar.gz", repos = NULL, type = "source")
```

## Quick start

``` r
library(Sphinx)

df <- prepare_data(meta, celltype_col = "celltype")
edges <- build_spatial_network(df, method = "auto")
feat <- calculate_neighborhood_features(df, edges)
clus <- cluster_neighborhoods(feat, edges, method = "kmeans", k = 10)

visualize_spatial_distribution(df)
visualize_spatial_network(clus, edges, edge_mode = "top", top_n = 1500, point_alpha = 1)
```

## Learn more

Use **Get started** and **Tutorials** for workflow, preprocessing,
annotation, spatial networks, and functional analysis.
