# Visualize spatial distribution of cell types

Generates spatial scatter plot showing cell type distribution:

1.  Requires prior extraction of spatial coordinates

2.  Uses custom color palette for cell type distinction

3.  Maintains spatial aspect ratio for accurate representation

4.  Optimizes point size and transparency for dense tissues

5.  Enhances legend readability for complex annotations

## Usage

``` r
plot_spatial_distribution(
  seurat_obj,
  save_path = "celltype_spatial_plot.pdf",
  width = 10,
  height = 8,
  point.size = 0.5
)
```

## Arguments

- seurat_obj:

  Seurat object with spatial coordinates

- save_path:

  Output file path (default: "celltype_spatial_plot.pdf")

- width:

  Plot width in inches (default: 10)

- height:

  Plot height in inches (default: 8)

- point.size:

  numeric, size of scatter points (default 0.5)

## Value

ggplot object containing spatial distribution plot
