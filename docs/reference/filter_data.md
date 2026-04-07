# Filter low-quality cells from spatial data

Performs two-step quality control filtering:

1.  Removes cells with outlier total counts (MAD-based)

2.  Removes cells with low detected feature counts (quantile-based)

## Usage

``` r
filter_data(obj, nCount_mad_threshold = 3, nFeature_quantile_threshold = 0.05)
```

## Arguments

- obj:

  Seurat object containing spatial data

- nCount_mad_threshold:

  MAD threshold for total count filtering (default: 3)

- nFeature_quantile_threshold:

  Quantile threshold for detected features (default: 0.05)

## Value

Filtered Seurat object
