# Visualize distance heatmap between cell types

Visualize distance heatmap between cell types

## Usage

``` r
visualize_distance_heatmap(
  dist_result,
  save_path = NULL,
  width = 12,
  height = 10,
  base_size = 14,
  show_values = TRUE
)
```

## Arguments

  - dist\_result:
    
    Distance matrix result from calculate\_celltype\_distances()

  - save\_path:
    
    Output file path (optional)

  - width:
    
    Plot width in inches (default: 12)

  - height:
    
    Plot height in inches (default: 10)

  - base\_size:
    
    Base font size in points (default: 14; minimum 8)

  - show\_values:
    
    Whether to print distance values in cells (default: TRUE)

## Value

ggplot object and saves plot to file if save\_path provided

## Examples

``` r
# \donttest{
df <- prepare_data(Sphinx:::.sphinx_example_df(40))
dist_res <- calculate_celltype_distances(df, celltype_col = "celltype")
p <- visualize_distance_heatmap(dist_res, save_path = tempfile(fileext = ".pdf"))
class(p)
# }
```
