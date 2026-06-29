# Workflow

## Overview

**Sphinx** is an R package for **spatial proteomics analysis**. It
provides an end-to-end workflow built on
[Seurat](https://satijalab.org/seurat/) for:

- Data preprocessing and quality control
- Cell-type annotation and marker discovery
- Spatial neighborhood modeling and interaction analysis
- Functional enrichment and publication-ready visualization

Sphinx supports common spatial proteomics platforms (e.g., CODEX, CyCIF,
IMC) and accepts tabular count matrices with spatial coordinates.

## Installation

Install the development version from GitHub (see
[`vignette("installation", package = "Sphinx")`](https://mongi126.github.io/Sphinx/articles/installation.md)
for details):

``` r

devtools::install_github("mongi126/Sphinx")
library(Sphinx)
packageVersion("Sphinx")
```

## Quick start

The typical analysis pipeline follows four sequential modules. Each step
produces a Seurat object that feeds into the next:

``` r

library(Sphinx)
library(Seurat)

# 1. Preprocessing
obj <- load_spatial_data("your_data.csv")
obj <- filter_data(obj)
obj <- process_data(obj)

# 2. Cell annotation
markers <- find_top_markers(obj)
obj     <- annotate_celltypes(obj, markers)

# 3. Spatial network
obj <- prepare_data(obj)
net <- build_spatial_network(obj, method = "knn")
obj <- calculate_neighborhood_features(obj, net)

# 4. Functional analysis
de_results <- perform_differential_expression(obj)
enrichment <- perform_cluster_enrichment(de_results)
```

For platform-specific parameters and visualization options, see the
module tutorials linked below.

## Analysis workflow

The overall workflow is summarized below:

![](workflow.jpg)

### Module 1. Data preprocessing

- Import spatial data and standardize format
  ([`load_spatial_data()`](https://mongi126.github.io/Sphinx/reference/load_spatial_data.md))
- Filter low-quality cells and proteins
  ([`filter_data()`](https://mongi126.github.io/Sphinx/reference/filter_data.md))
- Normalize, scale, and cluster
  ([`process_data()`](https://mongi126.github.io/Sphinx/reference/process_data.md))
- Extract spatial coordinates
  ([`extract_spatial_coordinates()`](https://mongi126.github.io/Sphinx/reference/extract_spatial_coordinates.md))

**Tutorial:**
[`vignette("preprocessing", package = "Sphinx")`](https://mongi126.github.io/Sphinx/articles/preprocessing.md)

### Module 2. Cell annotation

- Identify cluster-specific marker proteins
  ([`find_top_markers()`](https://mongi126.github.io/Sphinx/reference/find_top_markers.md))
- Visualize markers on UMAP and spatial maps
- Annotate cell types automatically or manually
  ([`annotate_celltypes()`](https://mongi126.github.io/Sphinx/reference/annotate_celltypes.md))

**Tutorial:**
[`vignette("annotation", package = "Sphinx")`](https://mongi126.github.io/Sphinx/articles/annotation.md)

### Module 3. Spatial neighborhood and network

- Build spatial graphs (kNN, Delaunay, radius, or window methods)
- Compute neighborhood composition and spatial metrics
- Analyze cell–cell interactions and spatial organization

**Tutorial:**
[`vignette("spatial-network", package = "Sphinx")`](https://mongi126.github.io/Sphinx/articles/spatial-network.md)

### Module 4. Functional analysis

- Perform differential protein expression across clusters
- Run pathway enrichment (GO, KEGG, Reactome)
- Visualize enrichment and volcano plots

**Tutorial:**
[`vignette("functional", package = "Sphinx")`](https://mongi126.github.io/Sphinx/articles/functional.md)

## Getting help

``` r

# Package documentation
help(package = "Sphinx")

# Function reference
?load_spatial_data
?build_spatial_network
```

Report bugs and request features on [GitHub
Issues](https://github.com/mongi126/Sphinx/issues).
