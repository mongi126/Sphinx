# Visualize distance heatmap between cell types

Visualize distance heatmap between cell types

## Usage

``` r
visualize_distance_heatmap(
  dist_result,
  save_path = NULL,
  width = 12,
  height = 10
)
```

## Arguments

- dist_result:

  Distance matrix result from calculate_celltype_distances()

- save_path:

  Output file path (optional)

- width:

  Plot width in inches (default: 12)

- height:

  Plot height in inches (default: 10)

## Value

ggplot object and saves plot to file if save_path provided
