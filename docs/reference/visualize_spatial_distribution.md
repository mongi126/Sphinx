# Visualize spatial cell type distribution

Visualize spatial cell type distribution

## Usage

``` r
visualize_spatial_distribution(
  df,
  x_col = "X",
  y_col = "Y",
  celltype_col = "celltype",
  point_size = 1.5,
  point_alpha = 0.6,
  point_shape = 16,
  legend_point_size = 3,
  title = NULL,
  legend.position = "right",
  save_path = NULL,
  width = 10,
  height = 8
)
```

## Arguments

- df:

  Spatial data

- x_col:

  X coordinate column name (default: "X")

- y_col:

  Y coordinate column name (default: "Y")

- celltype_col:

  Cell type column name (default: "celltype")

- point_size:

  Point size (default: 1.5)

- point_alpha:

  Point transparency (default: 0.6)

- point_shape:

  Point shape (default: 16)

- legend_point_size:

  Legend point size (default: 3)

- title:

  Plot title

- legend.position:

  Legend position (default: "right")

- save_path:

  Output file path (optional)

- width:

  Plot width in inches (default: 10)

- height:

  Plot height in inches (default: 8)

## Value

ggplot object and saves plot to file if save_path provided
