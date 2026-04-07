# Annotate cell types based on cluster IDs

Assigns biological cell type annotations to clusters:

1.  Validates equal length of cluster_ids and celltype_labels

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

- seurat_obj:

  Seurat object with cluster assignments

- cluster_ids:

  Vector of cluster IDs to annotate

- celltype_labels:

  Vector of cell type labels corresponding to cluster_ids

- cluster_column:

  Metadata column containing cluster IDs (default: "seurat_clusters")

## Value

Seurat object with added "celltype" metadata
