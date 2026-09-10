# Annotate cell types based on cluster IDs

Assigns biological cell type annotations to clusters:

1.  Validates equal length of cluster\_ids and celltype\_labels

2.  Creates cluster-to-celltype mapping

3.  Adds "celltype" column to object metadata

4.  Handles both numeric and character cluster IDs

## Usage

``` r
annotate_celltypes(
  seurat_obj,
  cluster_ids,
  celltype_labels,
  cluster_column = "seurat_clusters"
)
```

## Arguments

  - seurat\_obj:
    
    Seurat object with cluster assignments

  - cluster\_ids:
    
    Vector of cluster IDs to annotate

  - celltype\_labels:
    
    Vector of cell type labels corresponding to cluster\_ids

  - cluster\_column:
    
    Metadata column containing cluster IDs (default: "seurat\_clusters")

## Value

Seurat object with added "celltype" metadata

## Examples

``` r
# \donttest{
obj <- Sphinx:::.sphinx_example_seurat(20)
cls <- levels(obj@meta.data$seurat_clusters)
obj <- annotate_celltypes(obj, cls, paste0("Type", seq_along(cls)))
table(obj@meta.data$celltype)
# }
```
