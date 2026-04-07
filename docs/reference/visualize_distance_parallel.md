# Visualize distance relationships using parallel coordinates plot

Visualize distance relationships using parallel coordinates plot

## Usage

``` r
visualize_distance_parallel(
  dist_result,
  save_path = NULL,
  width = 14,
  height = 8
)
```

## Arguments

- dist_result:

  Distance matrix result from calculate_celltype_distances()

- save_path:

  Output file path (optional)

- width:

  Plot width in inches (default: 14)

- height:

  Plot height in inches (default: 8)

## Value

ggplot object and saves plot to file if save_path provided
