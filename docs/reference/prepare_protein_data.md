# Prepare protein expression data for functional analysis

Prepare protein expression data for functional analysis

## Usage

``` r
prepare_protein_data(cluster_df, expr_df)
```

## Arguments

- cluster_df:

  Data frame with neighborhood cluster assignments and spatial
  coordinates

- expr_df:

  Data frame with protein expression data (rows = cells, columns =
  proteins)

## Value

Merged data frame containing spatial, cluster, and protein expression
data
