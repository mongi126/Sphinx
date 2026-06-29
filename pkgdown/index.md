# Sphinx

**Sphinx** is an R package for **spatial proteomics analysis**. It provides an end-to-end workflow from raw counts to biological interpretation — including preprocessing, cell annotation, spatial neighborhood modeling, and functional enrichment.

Supported platforms include CODEX, CyCIF, IMC, and other tabular spatial proteomics formats with x/y coordinates.

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

See the [Workflow](articles/workflow.html) article for a guided overview, or jump directly to module tutorials below.

## Workflow

### Module 1. Data preprocessing

-   Data import and format standardization ([`load_spatial_data()`](reference/load_spatial_data.html))
-   Filtering of low-quality cells and proteins ([`filter_data()`](reference/filter_data.html))
-   Normalization, scaling, and clustering ([`process_data()`](reference/process_data.html))
-   Extraction of spatial coordinates ([`extract_spatial_coordinates()`](reference/extract_spatial_coordinates.html))

**Tutorial:** [Data Preprocessing](articles/preprocessing.html)

### Module 2. Cell annotation

-   Discovery of cell-type-specific marker proteins ([`find_top_markers()`](reference/find_top_markers.html))
-   UMAP and spatial visualization
-   Automatic or semi-automatic cell type annotation ([`annotate_celltypes()`](reference/annotate_celltypes.html))

**Tutorial:** [Cell Annotation](articles/annotation.html)

### Module 3. Spatial neighborhood and network

-   Auto neighborhood network construction — kNN, Delaunay, radius, and window methods ([`build_spatial_network()`](reference/build_spatial_network.html))
-   Calculation of neighborhood features ([`calculate_neighborhood_features()`](reference/calculate_neighborhood_features.html))
-   Clustering and spatial visualization
-   Cell–cell interaction analysis ([`analyze_spatial_interactions()`](reference/analyze_spatial_interactions.html))

**Tutorial:** [Spatial Network Analysis](articles/spatial-network.html)

### Module 4. Functional analysis

-   Differential protein analysis ([`perform_differential_expression()`](reference/perform_differential_expression.html))
-   Cluster-specific pathway enrichment — GO, KEGG, Reactome ([`perform_cluster_enrichment()`](reference/perform_cluster_enrichment.html))
-   Visualization of enrichment results ([`plot_enrichment_results()`](reference/plot_enrichment_results.html))

**Tutorial:** [Functional Analysis](articles/functional.html)

## Getting help

-   **Function reference:** [All functions](reference/index.html)
-   **Report issues:** [GitHub Issues](https://github.com/mongi126/Sphinx/issues)
-   **In R:** `help(package = "Sphinx")`
