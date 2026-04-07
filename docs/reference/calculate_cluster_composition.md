# Calculate cluster composition metrics

Calculate cluster composition metrics

## Usage

``` r
calculate_cluster_composition(
  df,
  cluster_col = "Neighborhood_Cluster",
  celltype_col = "celltype"
)
```

## Arguments

- df:

  Spatial data

- cluster_col:

  Cluster column name (default: "Neighborhood_Cluster")

- celltype_col:

  Cell type column name (default: "celltype")

## Value

Data frame with cluster composition statistics
