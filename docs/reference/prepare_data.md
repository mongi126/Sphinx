# Prepare spatial data for network analysis

Prepares spatial data for network analysis by:

1.  Converting to data.table if needed

2.  Ensuring Cell_ID column exists

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

- cell_id_col:

  Column name for cell IDs (default: NULL, auto-detect)

- x_col:

  Column name for X coordinates (default: "X")

- y_col:

  Column name for Y coordinates (default: "Y")

- celltype_col:

  Column name for cell types (default: "celltype")

## Value

data.table with standardized structure for spatial analysis
