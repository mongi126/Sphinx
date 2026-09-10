# Visualize cell types in UMAP space

Creates UMAP visualization colored by annotated cell types:

1.  Uses custom color palette for distinct cell types

2.  Labels clusters with centered text boxes

3.  Maintains aspect ratio for proper spatial representation

4.  Removes default title for cleaner presentation

## Usage

``` r
plot_annotated_umap(
  seurat_obj,
  save_path = "celltype_umap.pdf",
  width = 11,
  height = 8
)
```

## Arguments

  - seurat\_obj:
    
    Seurat object with cell type annotations

  - save\_path:
    
    Output file path (default: "celltype\_umap.pdf")

  - width:
    
    Plot width in inches (default: 11)

  - height:
    
    Plot height in inches (default: 8)

## Value

ggplot object containing annotated UMAP

## Examples

``` r
# \donttest{
obj <- Sphinx:::.sphinx_example_seurat(25)
cls <- levels(obj@meta.data$seurat_clusters)
obj <- annotate_celltypes(obj, cls, paste0("Type", seq_along(cls)))
p <- plot_annotated_umap(obj, save_path = tempfile(fileext = ".pdf"))
class(p)
# }
```
