# Create heatmap of cluster marker expression

Generates a heatmap visualizing average marker expression per cluster:

1.  Extracts unique markers from top_markers data frame

2.  Computes average expression per cluster

3.  Applies row-wise Z-score normalization

4.  Uses Red-Blue color scheme for intuitive visualization

5.  Saves high-quality PDF output

## Usage

``` r
plot_marker_heatmap(
  seurat_obj,
  top_markers,
  group_by = "seurat_clusters",
  assay = "Spatial",
  save_path = "marker_heatmap.pdf",
  width = 12,
  height = 8
)
```

## Arguments

- seurat_obj:

  Seurat object with cluster annotations

- top_markers:

  Data frame from find_top_markers

- group_by:

  Metadata column for grouping (default: "seurat_clusters")

- assay:

  Assay containing expression data (default: "Spatial")

- save_path:

  Output file path (default: "marker_heatmap.pdf")

- width:

  Plot width in inches (default: 12)

- height:

  Plot height in inches (default: 8)

## Value

Heatmap plot object
