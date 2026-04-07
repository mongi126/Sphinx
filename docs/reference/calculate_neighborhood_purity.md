# Calculate neighborhood purity using flexible neighbor definitions

Calculate neighborhood purity using flexible neighbor definitions

## Usage

``` r
calculate_neighborhood_purity(
  df,
  x_col = "X",
  y_col = "Y",
  celltype_col = "celltype",
  method = c("window", "radius", "knn", "delaunay"),
  k = 30,
  radius = 30,
  min_cells = 5,
  verbose = TRUE
)
```

## Arguments

- df:

  A data.frame or data.table containing spatial coordinates and cell
  type labels.

- x_col:

  Character, name of the X-coordinate column.

- y_col:

  Character, name of the Y-coordinate column.

- celltype_col:

  Character, name of the cell type column.

- method:

  Character, neighbor definition: "window", "radius", "knn", or
  "delaunay".

- k:

  Integer, number of nearest neighbors for method = "knn" (ignored
  otherwise).

- radius:

  Numeric, distance threshold for methods "radius" and "window" (ignored
  otherwise).

- min_cells:

  Integer, minimum number of neighbors required to compute purity.

- verbose:

  Logical, print progress messages.

## Value

A data.table with an added column `Neighborhood_Purity`.
