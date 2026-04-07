# Visualize clustering results

Generates and saves two key visualizations:

1.  UMAP projection showing cell clusters

2.  Spatial scatter plot showing cluster distribution in tissue context

## Usage

``` r
visualize_results(obj, save_dir = "./")
```

## Arguments

- obj:

  Processed Seurat object with spatial coordinates

- save_dir:

  Output directory (default: "./")

## Value

List containing UMAP and spatial plot objects
