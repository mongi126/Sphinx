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

  - sample\_size:
    
    Maximum cells to sample for efficiency (default: 1000)

  - multiplier:
    
    Factor for recommended radius (default: 2.0)

  - x\_col:
    
    Column name for X coordinates (default: "X")

  - y\_col:
    
    Column name for Y coordinates (default: "Y")

  - k\_nn:
    
    Nearest-neighbor rank used for distance summary (default: 10)

## Value

List with distance statistics and recommended parameters

## Examples

``` r
df <- prepare_data(Sphinx:::.sphinx_example_df(50))
params <- calculate_optimal_radius(df)
names(params)
```
