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

## Examples

``` r
# \donttest{
obj <- Sphinx:::.sphinx_example_seurat(20)
obj <- extract_spatial_coordinates(obj)
head(obj@meta.data[, c("X", "Y")])
# }
```
