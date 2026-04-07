# Calculate distances between cell types

Calculate distances between cell types

## Usage

``` r
calculate_celltype_distances(
  df,
  x_col = "X",
  y_col = "Y",
  celltype_col = "annotation"
)
```

## Arguments

- df:

  Spatial data with coordinates and cell types

- x_col:

  Column name for X coordinates (default: "X")

- y_col:

  Column name for Y coordinates (default: "Y")

- celltype_col:

  Column name for cell types (default: "annotation")

## Value

List with distance matrix and distributions
