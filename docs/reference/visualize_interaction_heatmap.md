# Visualize cell-cell interaction matrix

Visualize cell-cell interaction matrix

## Usage

``` r
visualize_interaction_heatmap(
  interaction_matrix,
  transform = TRUE,
  color_palette = viridis::inferno,
  save_path = NULL,
  width = 10,
  height = 8
)
```

## Arguments

- interaction_matrix:

  Interaction matrix from analyze_spatial_interactions()

- transform:

  Apply log2 transformation (default: TRUE)

- color_palette:

  Color palette function (default: viridis::inferno)

- save_path:

  Output file path (optional)

- width:

  Plot width in inches (default: 10)

- height:

  Plot height in inches (default: 8)

## Value

pheatmap object and saves plot to file if save_path provided
