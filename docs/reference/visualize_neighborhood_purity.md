# Visualize neighborhood purity

Visualize neighborhood purity

## Usage

``` r
visualize_neighborhood_purity(
  df,
  x_col = "X",
  y_col = "Y",
  max_points = 10000,
  point_size = 1.5,
  point_alpha = 0.8,
  point_shape = 16,
  title = NULL,
  legend.position = "right",
  save_path = NULL,
  width = 10,
  height = 8
)
```

## Arguments

- df:

  Spatial data with Neighborhood_Purity column

- x_col:

  X coordinate column name (default: "X")

- y_col:

  Y coordinate column name (default: "Y")

- max_points:

  Maximum points to plot (default: 10000)

- point_size:

  Point size (default: 1.5)

- point_alpha:

  Point transparency (default: 0.8)

- point_shape:

  Point shape (default: 16)

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
