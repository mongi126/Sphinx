# Calculate optimal spatial analysis parameters

Calculate optimal spatial analysis parameters

## Usage

``` r
calculate_optimal_radius(
  df,
  sample_size = 1000,
  multiplier = 1.5,
  x_col = "X",
  y_col = "Y"
)
```

## Arguments

- df:

  Spatial data with coordinates

- sample_size:

  Maximum cells to sample for efficiency (default: 1000)

- multiplier:

  Factor for recommended radius (default: 1.5)

- x_col:

  Column name for X coordinates (default: "X")

- y_col:

  Column name for Y coordinates (default: "Y")

## Value

List with distance statistics and recommended parameters
