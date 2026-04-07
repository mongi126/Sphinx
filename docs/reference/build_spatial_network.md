# Construct spatial neighborhoods from coordinates with automatic method selection

This function builds a cell-cell spatial network using one of four
methods ("radius", "knn", "delaunay", "window"). When `method = "auto"`,
it selects the method based on robust spatial metrics computed from the
data.

## Usage

``` r
build_spatial_network(
  df,
  method = "auto",
  k = NULL,
  radius = NULL,
  window_size = NULL,
  max_edges = 1e+06,
  celltype_col = "annotation",
  verbose = TRUE
)
```

## Arguments

- df:

  data.frame/data.table with columns: X, Y, Cell_ID; optional cell type
  column

- method:

  one of c("auto","radius","knn","delaunay","window")

- k:

  integer, neighbors for knn/window methods (default chosen
  automatically)

- radius:

  numeric, search radius for radius method (auto if NULL)

- window_size:

  numeric, grid side length for window method (auto if NULL)

- max_edges:

  integer, hard cap on number of returned edges (may be downsampled)

- celltype_col:

  character, column with cell-type labels (optional)

- verbose:

  logical, print decisions and key metrics

## Value

data.table with columns: from, to, dist (and optional context columns)
