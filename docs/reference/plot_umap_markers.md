# Visualize marker expression in UMAP space

Generates UMAP visualization of marker expression patterns:

1.  Validates UMAP reduction existence in Seurat object

2.  Determines appropriate assay for expression data

3.  Creates individual FeaturePlot for each marker

4.  Ensures numeric expression matrix for proper visualization

5.  Combines plots using patchwork for multi-panel output

6.  Saves high-quality PDF with automatic dimension calculation

## Usage

``` r
plot_umap_markers(
  seurat_obj,
  markers = NULL,
  assay = NULL,
  ncol = 3,
  point_size = 1.2,
  alpha = 0.9,
  color_low = "lightgrey",
  color_high = "red",
  output_file = "umap_markers_plot.pdf",
  width = NULL,
  height = NULL
)
```

## Arguments

  - seurat\_obj:
    
    Seurat object with UMAP reduction

  - markers:
    
    Vector of marker genes/proteins to visualize (default: NULL)

  - assay:
    
    Assay name containing expression data (default: NULL)

  - ncol:
    
    Number of columns for multi-plot layout (default: 3)

  - point\_size:
    
    Point size for UMAP scatter plot (default: 1.2)

  - alpha:
    
    Transparency level for points (default: 0.9)

  - color\_low:
    
    Color for low expression values (default: "lightgrey")

  - color\_high:
    
    Color for high expression values (default: "red")

  - output\_file:
    
    Output file path (default: "umap\_markers\_plot.pdf")

  - width:
    
    Plot width in inches (default: NULL, auto-calculated)

  - height:
    
    Plot height in inches (default: NULL, auto-calculated)

## Value

Combined ggplot object with UMAP marker expression plots

## Examples

``` r
# \donttest{
obj <- Sphinx:::.sphinx_example_seurat(25)
p <- plot_umap_markers(obj, markers = rownames(obj)[1:2],
  output_file = tempfile(fileext = ".pdf"))
class(p)
# }
```
