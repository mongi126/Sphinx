# Calculate neighborhood composition features

Calculate neighborhood composition features

## Usage

``` r
calculate_neighborhood_features(
  df,
  edges,
  cell_id_col = "Cell_ID",
  celltype_col = "celltype"
)
```

## Arguments

- df:

  Spatial data with cell types

- edges:

  Data.frame of edges from build_spatial_network()

- cell_id_col:

  Column name for cell IDs (default: "Cell_ID")

- celltype_col:

  Column name for cell types (default: "celltype")

## Value

Enhanced data.table with neighborhood type proportions
