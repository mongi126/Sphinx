# Cluster neighborhoods using combined spatial and compositional features

Cluster neighborhoods using combined spatial and compositional features

## Usage

``` r
cluster_neighborhoods(
  feature_df,
  spatial_edges = NULL,
  method = "kmeans",
  k = 10,
  use_pca = TRUE,
  var_threshold = 0.9,
  n_components = NULL,
  min_cluster_size = 5,
  cluster_colname = "Neighborhood_Cluster"
)
```

## Arguments

- feature_df:

  Data.table with neighborhood features

- spatial_edges:

  Spatial edges from build_spatial_network()

- method:

  Clustering method ("kmeans", "hdbscan", or "louvain")

- k:

  Number of clusters (for kmeans)

- use_pca:

  Whether to use PCA for dimensionality reduction

- var_threshold:

  Variance threshold for PCA components

- n_components:

  Explicit number of PCA components (overrides var_threshold)

- min_cluster_size:

  Minimum points per cluster (for hdbscan)

- cluster_colname:

  Name for the output cluster column

## Value

Data.table with cluster assignments in cluster_colname
