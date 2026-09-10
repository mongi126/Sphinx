# Calculate spatial metrics for method selection

Calculate spatial metrics for method selection

## Usage

``` r
.calc_spatial_metrics(df, celltype_col = NULL, x_col = NULL, y_col = NULL)
```

## Arguments

  - df:
    
    data.table with spatial data

  - celltype\_col:
    
    cell type column name

  - x\_col:
    
    optional X coordinate column name (auto-detects X/x)

  - y\_col:
    
    optional Y coordinate column name (auto-detects Y/y)

## Value

list of spatial metrics
