# Visualize spatial data using Voronoi diagrams

Visualize spatial data using Voronoi diagrams

## Usage

``` r
visualize_voronoi(
  df,
  x_col = "X",
  y_col = "Y",
  coloring = c("celltype", "neighborhood"),
  highlight_cluster = NULL,
  celltype_col = "celltype",
  neighborhood_col = "Neighborhood_Cluster",
  background_color = "gray90",
  highlight_alpha = 0.9,
  background_alpha = 0.3,
  celltype_palette = NULL,
  show_composition = TRUE,
  save_path = NULL,
  width = 12,
  height = 10
)
```

## Arguments

  - df:
    
    Spatial data with coordinates and annotations

  - x\_col:
    
    X coordinate column name (default: "X")

  - y\_col:
    
    Y coordinate column name (default: "Y")

  - coloring:
    
    Coloring method: "celltype" or "neighborhood" (default: "celltype")

  - highlight\_cluster:
    
    Specific cluster to highlight (optional)

  - celltype\_col:
    
    Cell type column name (default: "celltype")

  - neighborhood\_col:
    
    Neighborhood cluster column name (default: "Neighborhood\_Cluster")

  - background\_color:
    
    Background color for non-highlighted cells (default: "gray90")

  - highlight\_alpha:
    
    Alpha for highlighted cells (default: 0.9)

  - background\_alpha:
    
    Alpha for background cells (default: 0.3)

  - celltype\_palette:
    
    Custom color palette for cell types (optional)

  - show\_composition:
    
    Show cell type composition in highlight (default: TRUE)

  - save\_path:
    
    Output file path (optional)

  - width:
    
    Plot width in inches (default: 12)

  - height:
    
    Plot height in inches (default: 10)

## Value

ggplot object and saves plot to file if save\_path provided

## Examples

``` r
# \donttest{
df <- prepare_data(Sphinx:::.sphinx_example_df(40))
df$Neighborhood_Cluster <- sample(1:3, nrow(df), replace = TRUE)
p <- visualize_voronoi(df, celltype_col = "celltype",
  save_path = tempfile(fileext = ".pdf"))
class(p)
# }
```
