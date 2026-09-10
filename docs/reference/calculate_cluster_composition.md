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

  - cluster\_col:
    
    Cluster column name (default: "Neighborhood\_Cluster")

  - celltype\_col:
    
    Cell type column name (default: "celltype")

## Value

Data frame with cluster composition statistics

## Examples

``` r
df <- prepare_data(Sphinx:::.sphinx_example_df(40))
df$Neighborhood_Cluster <- sample(1:3, nrow(df), replace = TRUE)
comp <- calculate_cluster_composition(df)
head(comp)
```
