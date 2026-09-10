# Visualize cluster composition as heatmap

Visualize cluster composition as heatmap

## Usage

``` r
plot_composition_heatmap(
  composition_df,
  cluster_col = "Neighborhood_Cluster",
  celltype_col = "celltype",
  value_col = "proportion",
  save_path = NULL,
  width = 12,
  height = 10,
  cell_fontsize = 11
)
```

## Arguments

  - composition\_df:
    
    Composition data from calculate\_cluster\_composition()

  - cluster\_col:
    
    Cluster column name (default: "Neighborhood\_Cluster")

  - celltype\_col:
    
    Cell type column name (default: "celltype")

  - value\_col:
    
    Value column to visualize (default: "proportion")

  - save\_path:
    
    Output file path (optional)

  - width:
    
    Plot width in inches (default: 12)

  - height:
    
    Plot height in inches (default: 10)

  - cell\_fontsize:
    
    Font size for in-cell numbers (default: 11; minimum 8)

## Value

ComplexHeatmap object and saves plot to file if save\_path provided

## Examples

``` r
# \donttest{
if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  df <- prepare_data(Sphinx:::.sphinx_example_df(40))
  df$Neighborhood_Cluster <- sample(1:3, nrow(df), replace = TRUE)
  comp <- calculate_cluster_composition(df)
  plot_composition_heatmap(comp, save_path = tempfile(fileext = ".pdf"))
}
# }
```
