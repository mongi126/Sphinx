# Visualize marker expression in spatial coordinates

Generates spatial visualization of marker expression patterns:

1.  Extracts spatial coordinates from phenocycler image slot

2.  Validates marker presence in specified assay

3.  Creates individual spatial plots for each marker

4.  Combines plots using patchwork for multi-panel visualization

5.  Saves high-quality PDF output with automatic sizing

## Usage

``` r
plot_spatial_markers(
  obj,
  markers = NULL,
  assay = "Akoya",
  color_low = "lightgrey",
  color_high = "red",
  point_size = 1,
  alpha = 0.9,
  ncol = 3,
  output_file = "spatial_markers.pdf"
)
```

## Arguments

- obj:

  Seurat object with spatial coordinates

- markers:

  Vector of marker genes/proteins to visualize (default: NULL)

- assay:

  Assay name containing expression data (default: "Akoya")

- color_low:

  Color for low expression values (default: "lightgrey")

- color_high:

  Color for high expression values (default: "red")

- point_size:

  Point size for spatial scatter plot (default: 1.0)

- alpha:

  Transparency level for points (default: 0.9)

- ncol:

  Number of columns for multi-plot layout (default: 3)

- output_file:

  Output file path (default: "spatial_markers.pdf")

## Value

Combined ggplot object with spatial marker expression plots
