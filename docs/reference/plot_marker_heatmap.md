# Create heatmap of cluster marker expression

Generates a heatmap visualizing average marker expression per cluster:

1.  Extracts unique markers from top\_markers data frame

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

  - seurat\_obj:
    
    Seurat object with cluster annotations

  - top\_markers:
    
    Data frame from find\_top\_markers

  - group\_by:
    
    Metadata column for grouping (default: "seurat\_clusters")

  - assay:
    
    Assay containing expression data (default: "Spatial")

  - save\_path:
    
    Output file path (default: "marker\_heatmap.pdf")

  - width:
    
    Plot width in inches (default: 12)

  - height:
    
    Plot height in inches (default: 8)

## Value

Heatmap plot object

## Examples

``` r
# \donttest{
obj <- Sphinx:::.sphinx_example_seurat(25)
mk <- find_top_markers(obj, save_path = tempfile(fileext = ".csv"),
  assay = "RNA", min.pct = 0.1, logfc.threshold = 0.1)
if (nrow(mk) > 0) {
  plot_marker_heatmap(obj, mk, assay = "RNA",
    save_path = tempfile(fileext = ".pdf"))
}
# }
```
