# Extract spatial coordinates from Seurat object

Extracts spatial coordinates from reduction slot and adds them to object
metadata. Verifies successful extraction by printing coordinate sample.

## Usage

``` r
extract_spatial_coordinates(obj)
```

## Arguments

- obj:

  Seurat object with spatial data

## Value

Seurat object with coordinates added to metadata (if not already
present)
