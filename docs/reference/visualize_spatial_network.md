# Visualize spatial network with flexible edge selection

Visualize spatial network with flexible edge selection

## Usage

``` r
visualize_spatial_network(
  df,
  edges,
  celltype_col = "celltype",
  x_col = "X",
  y_col = "Y",
  edge_mode = c("all", "top", "random"),
  top_n = 1000,
  max_edges = 1000,
  point_size = 1.5,
  point_alpha = 0.7,
  edge_size_range = c(0.3, 2),
  edge_alpha_range = c(0.3, 0.9),
  edge_color = "darkred",
  show_points = TRUE,
  legend_point_size = 3,
  save_path = NULL,
  width = 12,
  height = 10
)
```

## Arguments

- df:

  Spatial data

- edges:

  Spatial network edges

- celltype_col:

  Cell type column name (default: "celltype")

- x_col:

  X coordinate column name (default: "X")

- y_col:

  Y coordinate column name (default: "Y")

- edge_mode:

  Edge selection mode: "all", "top", or "random" (default: "all")

- top_n:

  When edge_mode = "top", keep this many strongest edges (default: 1000)

- max_edges:

  When edge_mode = "random", sample this many edges (default: 1000)

- point_size:

  Point size (default: 1.5)

- point_alpha:

  Point transparency (default: 0.7)

- edge_size_range:

  Edge size range (default: c(0.3, 2.0))

- edge_alpha_range:

  Edge alpha range (default: c(0.3, 0.9))

- edge_color:

  Edge color (default: "darkred")

- show_points:

  Whether to show points (default: TRUE)

- legend_point_size:

  Legend point size (default: 3)

- save_path:

  Output file path (optional)

- width:

  Plot width in inches (default: 12)

- height:

  Plot height in inches (default: 10)

## Value

ggplot object and saves plot to file if save_path provided
