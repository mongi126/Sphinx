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

- save_path:

  Output file path (default: "elbow_plot.pdf")

## Value

ggplot object containing elbow plot
