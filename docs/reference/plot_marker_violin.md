# Visualize marker expression using violin plots

Generates violin plots showing expression distribution of marker
proteins:

1.  Creates one plot per marker protein

2.  Groups cells by cluster identity

3.  Removes legend for cleaner multi-plot output

4.  Saves composite plot to PDF

## Usage

``` r
plot_marker_violin(
  seurat_obj,
  markers,
  group_by = "seurat_clusters",
  assay = "Spatial",
  ncol = 3,
  save_path = "allmarker.pdf",
  width = 12,
  height = 30
)
```

## Arguments

- seurat_obj:

  Seurat object with cluster annotations

- markers:

  Vector of marker proteins to visualize

- group_by:

  Metadata column for grouping (default: "seurat_clusters")

- assay:

  Assay containing expression data (default: "Spatial")

- ncol:

  Number of columns for multi-plot layout (default: 3)

- save_path:

  Output file path (default: "allmarker.pdf")

- width:

  Plot width in inches (default: 12)

- height:

  Plot height in inches (default: 30)

## Value

ggplot object containing violin plots
