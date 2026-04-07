# Analyze spatial interactions between cell types

Analyze spatial interactions between cell types

## Usage

``` r
analyze_spatial_interactions(df, edges, celltype_col = "celltype")
```

## Arguments

- df:

  Spatial data with cell types

- edges:

  Spatial network edges (data.frame with "from" and "to")

- celltype_col:

  Column name for cell types in metadata.

## Value

List with interaction matrix and network graph
