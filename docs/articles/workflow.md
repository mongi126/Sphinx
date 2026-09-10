# Workflow

## Overview

**Sphinx** is an R toolkit for **spatial proteomics**. It organizes
common analysis steps into four focused modules that you can use
together or independently:

  - Preprocessing and quality control
  - Marker discovery and cell-type annotation
  - Spatial neighborhood graphs and interaction analysis
  - Differential proteins, enrichment, and publication-ready figures

Tutorials use **reg055\_A** from the public CODEX colorectal carcinoma
cohort of Schurch *et al.*, *Cell* (2020)
([data](https://data.mendeley.com/datasets/mpjzbtfgfr/1)): \~3.9k cells
and 58 proteins, with published `ClusterName` labels that we map onto
reclustered identities for the worked examples.

![](workflow.png)

## Installation

``` r
# remotes::install_github("mongi126/Sphinx")
# or: install.packages("Sphinx_1.0.1.tar.gz", repos = NULL, type = "source")
library(Sphinx)
packageVersion("Sphinx")
```

See `vignette("installation", package = "Sphinx")` for details.

## Quick start

With a table of `Cell_ID`, `X`, `Y`, and `celltype`:

``` r
library(Sphinx)
library(data.table)

df <- prepare_data(meta, cell_id_col = "Cell_ID", celltype_col = "celltype")

dist_result <- calculate_celltype_distances(df)
edges <- build_spatial_network(df, method = "auto")
feat <- calculate_neighborhood_features(df, edges)
clus <- cluster_neighborhoods(feat, edges, method = "kmeans", k = 10)

visualize_spatial_distribution(df, point_size = 0.45)
visualize_distance_heatmap(dist_result, show_values = TRUE)
visualize_spatial_network(
  clus, edges, edge_mode = "top", top_n = 1500, point_alpha = 1
)
```

## Modules

### 1\. Preprocessing

Load counts and coordinates, filter low-quality cells, normalize, and
cluster.

**Tutorial:** `vignette("preprocessing", package = "Sphinx")`

### 2\. Annotation

Find markers, review violin / heatmap / UMAP views, and assign readable
cell-type labels.

**Tutorial:** `vignette("annotation", package = "Sphinx")`

### 3\. Spatial networks

Build adaptive graphs, derive neighborhood clusters, and plot distances,
composition, purity, and interactions.

**Tutorial:** `vignette("spatial-network", package = "Sphinx")`

### 4\. Functional analysis

Test neighborhood-associated proteins (spatial-block aware), draw
volcano plots, and summarize enrichment.

**Tutorial:** `vignette("functional", package = "Sphinx")`

## Visualization defaults

  - Stable candy qualitative colors via `assign_celltype_colors()`
  - Soft pink sequential heatmaps; fonts at least 8 pt
  - `visualize_spatial_network()` draws opaque points with a white halo;
    use `zoom_center` / `zoom_radius` for local views
  - Volcano y-axis caps `-log10(adj.P)` (default 50)

## Getting help

  - Website: <https://mongi126.github.io/Sphinx/>
  - Issues: <https://github.com/mongi126/Sphinx/issues>
