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
    
    Data.frame of edges from build\_spatial\_network()

  - cell\_id\_col:
    
    Column name for cell IDs (default: "Cell\_ID")

  - celltype\_col:
    
    Column name for cell types (default: "celltype")

## Value

Enhanced data.table with neighborhood type proportions

## Examples

``` r
df <- prepare_data(Sphinx:::.sphinx_example_df(40))
edges <- build_spatial_network(df, method = "knn", n_neighbors = 5, verbose = FALSE)
feat <- calculate_neighborhood_features(df, edges)
names(feat)
```
