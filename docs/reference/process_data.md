# Process spatial data

Performs comprehensive data processing pipeline:

1.  Normalization (CLR recommended for proteomics)

2.  Feature selection (uses all detected features)

3.  Data scaling

4.  PCA dimensionality reduction

5.  Nearest-neighbor graph construction

6.  Cluster identification

7.  UMAP visualization

## Usage

``` r
process_data(
  obj,
  normalization.method = "CLR",
  margin = 2,
  dims = 1:10,
  resolution = 0.5
)
```

## Arguments

- obj:

  Seurat object containing filtered spatial data

- normalization.method:

  Normalization method (default: "CLR")

- margin:

  Margin for normalization (1=features, 2=cells) (default: 2)

- dims:

  Dimensions for reduction (default: 1:10)

- resolution:

  Clustering resolution (default: 0.5)

## Value

Processed Seurat object with dimensionality reduction and clustering
