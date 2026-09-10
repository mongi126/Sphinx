# Prepare spatial data for network analysis

Prepares spatial data for network analysis by:

1.  Converting to data.table if needed

2.  Ensuring Cell\_ID column exists

3.  Converting coordinates to numeric

4.  Removing invalid coordinates and missing cell types

## Usage

``` r
prepare_data(
  df,
  cell_id_col = NULL,
  x_col = "X",
  y_col = "Y",
  celltype_col = "celltype"
)
```

## Arguments

  - df:
    
    Data frame containing spatial data

  - cell\_id\_col:
    
    Column name for cell IDs (default: NULL, auto-detect)

  - x\_col:
    
    Column name for X coordinates (default: "X")

  - y\_col:
    
    Column name for Y coordinates (default: "Y")

  - celltype\_col:
    
    Column name for cell types (default: "celltype")

## Value

data.table with standardized structure for spatial analysis

## Examples

``` r
df <- Sphinx:::.sphinx_example_df(30)
df <- prepare_data(df)
head(df[, .(Cell_ID, X, Y, celltype)])
```
