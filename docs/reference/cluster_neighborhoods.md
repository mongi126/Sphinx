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

  - feature\_df:
    
    Data.table with neighborhood features

  - spatial\_edges:
    
    Spatial edges from build\_spatial\_network()

  - method:
    
    Clustering method ("kmeans", "hdbscan", or "louvain")

  - k:
    
    Number of clusters (for kmeans). This is **not** spatial neighbor
    count; use `n_neighbors` in `build_spatial_network()`.

  - use\_pca:
    
    Whether to use PCA for dimensionality reduction

  - var\_threshold:
    
    Variance threshold for PCA components

  - n\_components:
    
    Explicit number of PCA components (overrides var\_threshold)

  - min\_cluster\_size:
    
    Minimum points per cluster (for hdbscan)

  - cluster\_colname:
    
    Name for the output cluster column

## Value

Data.table with cluster assignments in cluster\_colname

## Examples

``` r
# \donttest{
df <- prepare_data(Sphinx:::.sphinx_example_df(40))
edges <- build_spatial_network(df, method = "knn", n_neighbors = 5, verbose = FALSE)
feat <- calculate_neighborhood_features(df, edges)
cl <- cluster_neighborhoods(feat, edges, method = "kmeans", k = 3)
"Neighborhood_Cluster" %in% names(cl)
# }
```
