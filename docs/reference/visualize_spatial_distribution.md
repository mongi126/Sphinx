# Visualize spatial cell type distribution

Visualize spatial cell type distribution

## Usage

``` r
visualize_spatial_distribution(
  df,
  x_col = "X",
  y_col = "Y",
  celltype_col = "celltype",
  point_size = 1.5,
  point_alpha = 0.85,
  point_shape = 16,
  legend_point_size = 3,
  title = NULL,
  legend.position = "right",
  base_size = 14,
  color_palette = NULL,
  all_levels = NULL,
  save_path = NULL,
  width = 10,
  height = 8
)
```

## Arguments

  - df:
    
    Spatial data

  - x\_col:
    
    X coordinate column name (default: "X")

  - y\_col:
    
    Y coordinate column name (default: "Y")

  - celltype\_col:
    
    Cell type column name (default: "celltype")

  - point\_size:
    
    Point size (default: 1.5)

  - point\_alpha:
    
    Point transparency (default: 0.6)

  - point\_shape:
    
    Point shape (default: 16)

  - legend\_point\_size:
    
    Legend point size (default: 3)

  - title:
    
    Plot title

  - legend.position:
    
    Legend position (default: "right")

  - base\_size:
    
    Base font size in points (default: 14; minimum 8)

  - color\_palette:
    
    Optional named color vector; if NULL uses assign\_celltype\_colors()

  - all\_levels:
    
    Optional full label set for stable colors when subsetting / zooming

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
p <- visualize_spatial_distribution(df, save_path = tempfile(fileext = ".pdf"))
class(p)
# }
```
