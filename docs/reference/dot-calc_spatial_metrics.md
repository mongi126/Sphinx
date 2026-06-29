# Calculate spatial metrics for method selection

Calculate spatial metrics for method selection

## Usage

``` r
.calc_spatial_metrics(df, celltype_col = NULL, x_col = NULL, y_col = NULL)
```

## Arguments

- df:

  data.table with spatial data

- celltype_col:

  cell type column name

- x_col:

  optional X coordinate column name (auto-detects X/x)

- y_col:

  optional Y coordinate column name (auto-detects Y/y)

## Value

list of spatial metrics
