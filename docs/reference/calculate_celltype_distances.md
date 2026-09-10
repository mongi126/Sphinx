# Calculate distances between cell types

Calculate distances between cell types

## Usage

``` r
calculate_celltype_distances(
  df,
  x_col = "X",
  y_col = "Y",
  celltype_col = "celltype"
)
```

## Arguments

  - df:
    
    Spatial data with coordinates and cell types

  - x\_col:
    
    Column name for X coordinates (default: "X")

  - y\_col:
    
    Column name for Y coordinates (default: "Y")

  - celltype\_col:
    
    Column name for cell types (default: "celltype")

## Value

List with distance matrix and distributions

## Examples

``` r
df <- prepare_data(Sphinx:::.sphinx_example_df(40))
dist_res <- calculate_celltype_distances(df, celltype_col = "celltype")
dim(dist_res$distance_matrix)
```
