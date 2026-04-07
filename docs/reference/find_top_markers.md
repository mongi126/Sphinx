# Identify top marker proteins for cell clusters

Performs differential expression analysis to identify cluster
biomarkers:

1.  Uses Seurat's FindAllMarkers with positive expression only

2.  Selects top 5 markers per cluster based on average log2FC

3.  Saves results to CSV file

## Usage

``` r
find_top_markers(
  seurat_obj,
  save_path = "top5proteins.csv",
  assay = "Spatial",
  min.pct = 0.25,
  logfc.threshold = 0.25,
  only.pos = TRUE
)
```

## Arguments

- seurat_obj:

  Seurat object containing clustered spatial data

- save_path:

  Output file path for marker results (default: "top5proteins.csv")

- assay:

  Assay name containing protein expression (default: "Spatial")

- min.pct:

  Minimum detection fraction threshold (default: 0.25)

- logfc.threshold:

  Minimum log-fold change threshold (default: 0.25)

- only.pos:

  Whether to return only positive markers (default: TRUE)

## Value

Data frame containing top 5 markers per cluster
