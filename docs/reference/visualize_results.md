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

  - save\_dir:
    
    Output directory (default: "./")

## Value

List containing UMAP and spatial plot objects

## Examples

``` r
# \donttest{
obj <- Sphinx:::.sphinx_example_seurat(25)
obj <- annotate_celltypes(obj,
  cluster_ids = levels(obj@meta.data$seurat_clusters),
  celltype_labels = paste0("Type", seq_along(levels(obj@meta.data$seurat_clusters)))
)
td <- tempdir()
plots <- visualize_results(obj, save_dir = td)
names(plots)
# }
```
