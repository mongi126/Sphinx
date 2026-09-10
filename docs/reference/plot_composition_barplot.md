# Visualize cluster composition as bar plot

Visualize cluster composition as bar plot

## Usage

``` r
plot_composition_barplot(
  composition_df,
  cluster_col = "Neighborhood_Cluster",
  celltype_col = "celltype",
  value_col = "proportion",
  save_path = NULL,
  width = 12,
  height = 8,
  base_size = 14
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
    
    Plot height in inches (default: 8)

  - base\_size:
    
    Base font size in points (default: 14; minimum 8)

## Value

ggplot object and saves plot to file if save\_path provided

## Examples

``` r
# \donttest{
df <- prepare_data(Sphinx:::.sphinx_example_df(40))
df$Neighborhood_Cluster <- sample(1:3, nrow(df), replace = TRUE)
comp <- calculate_cluster_composition(df)
p <- plot_composition_barplot(comp, save_path = tempfile(fileext = ".pdf"))
class(p)
# }
```
