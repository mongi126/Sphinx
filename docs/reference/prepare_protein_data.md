# Prepare protein expression data for functional analysis

Prepare protein expression data for functional analysis

## Usage

``` r
prepare_protein_data(cluster_df, expr_df)
```

## Arguments

  - cluster\_df:
    
    Data frame with neighborhood cluster assignments and spatial
    coordinates

  - expr\_df:
    
    Data frame with protein expression data (rows = cells, columns =
    proteins)

## Value

Merged data frame containing spatial, cluster, and protein expression
data

## Examples

``` r
cluster_df <- Sphinx:::.sphinx_example_df(10)
cluster_df$Neighborhood_Cluster <- sample(1:2, 10, replace = TRUE)
expr_df <- data.frame(PD.1 = rnorm(10), CD3 = rnorm(10),
  row.names = cluster_df$Cell_ID)
merged <- prepare_protein_data(cluster_df, expr_df)
nrow(merged)
```
