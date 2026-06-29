# Calculate optimal spatial analysis parameters

Calculate optimal spatial analysis parameters

## Usage

``` r
calculate_optimal_radius(
  df,
  sample_size = 1000,
  multiplier = 2,
  x_col = "X",
  y_col = "Y",
  k_nn = 10L
)
```

## Arguments

- df:

  Spatial data with coordinates

- sample_size:

  Maximum cells to sample for efficiency (default: 1000)

- multiplier:

  Factor for recommended radius (default: 2.0)

- x_col:

  Column name for X coordinates (default: "X")

- y_col:

  Column name for Y coordinates (default: "Y")

- k_nn:

  Nearest-neighbor rank used for distance summary (default: 10)

## Value

List with distance statistics and recommended parameters
