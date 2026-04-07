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

- seurat_obj:

  Seurat object with cell type annotations

- save_path:

  Output file path (default: "celltype_umap.pdf")

- width:

  Plot width in inches (default: 11)

- height:

  Plot height in inches (default: 8)

## Value

ggplot object containing annotated UMAP
