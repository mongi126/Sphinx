# Sphinx

**Sphinx** is an R package for **spatial proteomics analysis**. It
provides an end-to-end workflow from raw counts to biological
interpretation — including preprocessing, cell annotation, spatial
neighborhood modeling, and functional enrichment.

Supported platforms include CODEX, CyCIF, IMC, and other tabular spatial
proteomics formats with x/y coordinates.

## Installation

``` r

# Install development version from GitHub
devtools::install_github("mongi126/Sphinx")
library(Sphinx)
packageVersion("Sphinx")
```

## Quick start

``` r

library(Sphinx)
library(Seurat)

# 1. Preprocessing
obj <- load_spatial_data("your_data.csv")
obj <- filter_data(obj)
obj <- process_data(obj)

# 2. Cell annotation
markers <- find_top_markers(obj)
obj     <- annotate_celltypes(obj, cluster, markers)

# 3. Spatial network
obj <- prepare_data(obj)
net <- build_spatial_network(obj, method = "knn")
obj <- calculate_neighborhood_features(obj, net)

# 4. Functional analysis
de_results <- perform_differential_expression(obj)
enrichment <- perform_cluster_enrichment(de_results)
```

See the
[Workflow](https://mongi126.github.io/Sphinx/articles/workflow.md)
article for a guided overview, or jump directly to module tutorials
below.

## Workflow

### Module 1. Data preprocessing

- Data import and format standardization
  ([`load_spatial_data()`](https://mongi126.github.io/Sphinx/reference/load_spatial_data.md))
- Filtering of low-quality cells and proteins
  ([`filter_data()`](https://mongi126.github.io/Sphinx/reference/filter_data.md))
- Normalization, scaling, and clustering
  ([`process_data()`](https://mongi126.github.io/Sphinx/reference/process_data.md))
- Extraction of spatial coordinates
  ([`extract_spatial_coordinates()`](https://mongi126.github.io/Sphinx/reference/extract_spatial_coordinates.md))

**Tutorial:** [Data
Preprocessing](https://mongi126.github.io/Sphinx/articles/preprocessing.md)

### Module 2. Cell annotation

- Discovery of cell-type-specific marker proteins
  ([`find_top_markers()`](https://mongi126.github.io/Sphinx/reference/find_top_markers.md))
- UMAP and spatial visualization
- Automatic or semi-automatic cell type annotation
  ([`annotate_celltypes()`](https://mongi126.github.io/Sphinx/reference/annotate_celltypes.md))

**Tutorial:** [Cell
Annotation](https://mongi126.github.io/Sphinx/articles/annotation.md)

### Module 3. Spatial neighborhood and network

- Auto neighborhood network construction — kNN, Delaunay, radius, and
  window methods
  ([`build_spatial_network()`](https://mongi126.github.io/Sphinx/reference/build_spatial_network.md))
- Calculation of neighborhood features
  ([`calculate_neighborhood_features()`](https://mongi126.github.io/Sphinx/reference/calculate_neighborhood_features.md))
- Clustering and spatial visualization
- Cell–cell interaction analysis
  ([`analyze_spatial_interactions()`](https://mongi126.github.io/Sphinx/reference/analyze_spatial_interactions.md))

**Tutorial:** [Spatial Network
Analysis](https://mongi126.github.io/Sphinx/articles/spatial-network.md)

### Module 4. Functional analysis

- Differential protein analysis
  ([`perform_differential_expression()`](https://mongi126.github.io/Sphinx/reference/perform_differential_expression.md))
- Cluster-specific pathway enrichment — GO, KEGG, Reactome
  ([`perform_cluster_enrichment()`](https://mongi126.github.io/Sphinx/reference/perform_cluster_enrichment.md))
- Visualization of enrichment results
  ([`plot_enrichment_results()`](https://mongi126.github.io/Sphinx/reference/plot_enrichment_results.md))

**Tutorial:** [Functional
Analysis](https://mongi126.github.io/Sphinx/articles/functional.md)

## Getting help

- **Function reference:** [All
  functions](https://mongi126.github.io/Sphinx/reference/index.md)
- **Report issues:** [GitHub
  Issues](https://github.com/mongi126/Sphinx/issues)
- **In R:**
  [`help(package = "Sphinx")`](https://github.com/mongi126/Sphinx/reference)
