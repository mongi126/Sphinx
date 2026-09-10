# Load spatial data from various sources

Loads spatial data from either:

1.  An existing Seurat object (RDS file)

2.  Counts and coordinates CSV files

3.  Akoya platform data files

## Usage

``` r
load_spatial_data(
  filename = NULL,
  counts_file = NULL,
  coords_file = NULL,
  type = c("auto", "rds", "csv", "akoya"),
  akoya_type = "qupath",
  fov = "phenocycler",
  ...
)
```

## Arguments

  - filename:
    
    Path to input data file (RDS, CSV counts, or CSV coordinates)

  - counts\_file:
    
    Path to counts CSV file (if using separate counts and coordinates)

  - coords\_file:
    
    Path to coordinates CSV file (if using separate counts and
    coordinates)

  - type:
    
    Data source type: "auto", "rds", "csv", or "akoya" (default: "auto")

  - akoya\_type:
    
    Type of Akoya data: "qupath" or other formats (default: "qupath")

  - fov:
    
    Field of view specification (default: "phenocycler")

  - ...:
    
    Additional parameters passed to loading functions

## Value

Seurat object containing spatial data

## Examples

``` r
if (FALSE) { # \dontrun{
# From RDS or CSV files on disk:
# obj <- load_spatial_data(filename = "sample.rds")
# obj <- load_spatial_data(counts_file = "counts.csv", coords_file = "coords.csv")
} # }
```
