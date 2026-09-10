# Visualize neighborhood purity

Visualize neighborhood purity

## Usage

``` r
visualize_neighborhood_purity(
  df,
  x_col = "X",
  y_col = "Y",
  max_points = 10000,
  point_size = 1.5,
  point_alpha = 0.9,
  point_shape = 16,
  title = NULL,
  legend.position = "right",
  base_size = 14,
  save_path = NULL,
  width = 10,
  height = 8
)
```

## Arguments

  - df:
    
    Spatial data with Neighborhood\_Purity column

  - x\_col:
    
    X coordinate column name (default: "X")

  - y\_col:
    
    Y coordinate column name (default: "Y")

  - max\_points:
    
    Maximum points to plot (default: 10000)

  - point\_size:
    
    Point size (default: 1.5)

  - point\_alpha:
    
    Point transparency (default: 0.8)

  - point\_shape:
    
    Point shape (default: 16)

  - title:
    
    Plot title

  - legend.position:
    
    Legend position (default: "right")

  - base\_size:
    
    Base font size in points (default: 14; minimum 8)

  - save\_path:
    
    Output file path (optional)

  - width:
    
    Plot width in inches (default: 10)

  - height:
    
    Plot height in inches (default: 8)

## Value

ggplot object and saves plot to file if save\_path provided

## Examples

``` r
# \donttest{
df <- prepare_data(Sphinx:::.sphinx_example_df(40))
df <- calculate_neighborhood_purity(df, method = "knn", n_neighbors = 5, verbose = FALSE)
p <- visualize_neighborhood_purity(df, save_path = tempfile(fileext = ".pdf"))
class(p)
# }
```
