# Visualize distance relationships using parallel coordinates plot

Visualize distance relationships using parallel coordinates plot

## Usage

``` r
visualize_distance_parallel(
  dist_result,
  save_path = NULL,
  width = 14,
  height = 8,
  base_size = 14
)
```

## Arguments

  - dist\_result:
    
    Distance matrix result from calculate\_celltype\_distances()

  - save\_path:
    
    Output file path (optional)

  - width:
    
    Plot width in inches (default: 14)

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
dist_res <- calculate_celltype_distances(df, celltype_col = "celltype")
p <- visualize_distance_parallel(dist_res, save_path = tempfile(fileext = ".pdf"))
class(p)
# }
```
