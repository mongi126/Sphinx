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
  height = 8
)
```

## Arguments

- composition_df:

  Composition data from calculate_cluster_composition()

- cluster_col:

  Cluster column name (default: "Neighborhood_Cluster")

- celltype_col:

  Cell type column name (default: "celltype")

- value_col:

  Value column to visualize (default: "proportion")

- save_path:

  Output file path (optional)

- width:

  Plot width in inches (default: 12)

- height:

  Plot height in inches (default: 8)

## Value

ggplot object and saves plot to file if save_path provided
