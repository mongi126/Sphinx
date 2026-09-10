# Generate elbow plot for dimensionality reduction

Creates elbow plot to determine optimal number of PCA dimensions. Saves
plot to specified path and returns ggplot object.

## Usage

``` r
plot_elbow(obj, save_path = "elbow_plot.pdf")
```

## Arguments

  - obj:
    
    Processed Seurat object

  - save\_path:
    
    Output file path (default: "elbow\_plot.pdf")

## Value

ggplot object containing elbow plot

## Examples

``` r
# \donttest{
obj <- Sphinx:::.sphinx_example_seurat(30)
p <- plot_elbow(obj, save_path = tempfile(fileext = ".pdf"))
class(p)
# }
```
