# Visualize spatial data using Voronoi diagrams

Visualize spatial data using Voronoi diagrams

## Usage

``` r
visualize_voronoi(
  df,
  x_col = "X",
  y_col = "Y",
  coloring = c("celltype", "neighborhood"),
  highlight_cluster = NULL,
  celltype_col = "annotation",
  neighborhood_col = "Neighborhood_Cluster",
  background_color = "gray90",
  highlight_alpha = 0.9,
  background_alpha = 0.3,
  celltype_palette = NULL,
  show_composition = TRUE,
  save_path = NULL,
  width = 12,
  height = 10
)
```

## Arguments

- df:

  Spatial data with coordinates and annotations

- x_col:

  X coordinate column name (default: "X")

- y_col:

  Y coordinate column name (default: "Y")

- coloring:

  Coloring method: "celltype" or "neighborhood" (default: "celltype")

- highlight_cluster:

  Specific cluster to highlight (optional)

- celltype_col:

  Cell type column name (default: "annotation")

- neighborhood_col:

  Neighborhood cluster column name (default: "Neighborhood_Cluster")

- background_color:

  Background color for non-highlighted cells (default: "gray90")

- highlight_alpha:

  Alpha for highlighted cells (default: 0.9)

- background_alpha:

  Alpha for background cells (default: 0.3)

- celltype_palette:

  Custom color palette for cell types (optional)

- show_composition:

  Show cell type composition in highlight (default: TRUE)

- save_path:

  Output file path (optional)

- width:

  Plot width in inches (default: 12)

- height:

  Plot height in inches (default: 10)

## Value

ggplot object and saves plot to file if save_path provided
