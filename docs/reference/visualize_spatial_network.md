# Visualize spatial network with flexible edge selection

Draws **real pairwise cell-cell edges** using spatial coordinates.
Supports local zoom via `xlim`/`ylim` or `zoom_center` + `zoom_radius`.
When zooming, edges with both endpoints inside the window are kept so
local connectivity is shown faithfully.

## Usage

``` r
visualize_spatial_network(
  df,
  edges,
  celltype_col = "celltype",
  x_col = "X",
  y_col = "Y",
  edge_mode = c("all", "top", "random"),
  top_n = 1000,
  max_edges = 1000,
  xlim = NULL,
  ylim = NULL,
  zoom_center = NULL,
  zoom_radius = NULL,
  point_size = 1.5,
  point_alpha = 1,
  edge_size_range = c(0.3, 2),
  edge_alpha_range = c(0.15, 0.45),
  edge_color = "grey35",
  show_points = TRUE,
  legend_point_size = 3,
  base_size = 14,
  title = NULL,
  save_path = NULL,
  width = 12,
  height = 10
)
```

## Arguments

  - df:
    
    Spatial data

  - edges:
    
    Spatial network edges

  - celltype\_col:
    
    Cell type column name (default: "celltype")

  - x\_col:
    
    X coordinate column name (default: "X")

  - y\_col:
    
    Y coordinate column name (default: "Y")

  - edge\_mode:
    
    Edge selection mode: "all", "top", or "random" (default: "all")

  - top\_n:
    
    When edge\_mode = "top", keep this many strongest edges (default:
    1000)

  - max\_edges:
    
    When edge\_mode = "random", sample this many edges (default: 1000)

  - xlim:
    
    Optional x-axis limits for local zoom, e.g. `c(xmin, xmax)`

  - ylim:
    
    Optional y-axis limits for local zoom, e.g. `c(ymin, ymax)`

  - zoom\_center:
    
    Optional center `c(x, y)` for circular/rectangular zoom

  - zoom\_radius:
    
    Half-width of zoom window around `zoom_center` (same units as
    coords)

  - point\_size:
    
    Point size (default: 1.5)

  - point\_alpha:
    
    Point transparency (default: 1; keep opaque so cell-type colors stay
    clear over edges)

  - edge\_size\_range:
    
    Edge size range (default: c(0.3, 2.0))

  - edge\_alpha\_range:
    
    Edge alpha range (default: c(0.3, 0.9))

  - edge\_color:
    
    Edge color (default: "grey35")

  - show\_points:
    
    Whether to show points (default: TRUE)

  - legend\_point\_size:
    
    Legend point size (default: 3)

  - base\_size:
    
    Base font size (default: 14)

  - title:
    
    Plot title

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
edges <- build_spatial_network(df, method = "knn", n_neighbors = 5, verbose = FALSE)
p <- visualize_spatial_network(df, edges, save_path = tempfile(fileext = ".pdf"))
# Local zoom around tissue center:
# visualize_spatial_network(df, edges, zoom_center = c(500, 500), zoom_radius = 80)
class(p)
# }
```
