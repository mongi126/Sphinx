# Calculate neighborhood purity using flexible neighbor definitions

Neighbor definitions match `build_spatial_network()`: knn, radius ball,
Delaunay triangulation, or Stereopy-style sliding window + local kNN.

## Usage

``` r
calculate_neighborhood_purity(
  df,
  x_col = "X",
  y_col = "Y",
  celltype_col = "celltype",
  method = c("window", "radius", "knn", "delaunay"),
  n_neighbors = 10L,
  k = NULL,
  radius = NULL,
  window_size = NULL,
  window_slide_step = NULL,
  min_cells = 5,
  verbose = TRUE
)
```

## Arguments

  - df:
    
    A data.frame or data.table containing spatial coordinates and cell
    type labels.

  - x\_col:
    
    Character, name of the X-coordinate column.

  - y\_col:
    
    Character, name of the Y-coordinate column.

  - celltype\_col:
    
    Character, name of the cell type column.

  - method:
    
    Character, neighbor definition: "window", "radius", "knn", or
    "delaunay".

  - n\_neighbors:
    
    Integer, neighbor count for `knn` and `window`.

  - k:
    
    Integer, **deprecated** neighbor count. Use `n_neighbors`.

  - radius:
    
    Numeric search radius for `method = "radius"`. `NULL` uses the same
    1-NN x 3 default as `build_spatial_network()`.

  - window\_size:
    
    Sliding-window side length for `method = "window"`. `NULL` uses the
    same data default as `build_spatial_network()`.

  - window\_slide\_step:
    
    Window stride; `NULL` uses `window_size / 2`.

  - min\_cells:
    
    Integer, minimum number of neighbors required to compute purity.

  - verbose:
    
    Logical, print progress messages.

## Value

A data.table with an added column `Neighborhood_Purity`.

## Examples

``` r
df <- prepare_data(Sphinx:::.sphinx_example_df(40))
out <- calculate_neighborhood_purity(df, method = "knn", n_neighbors = 5, verbose = FALSE)
summary(out$Neighborhood_Purity)
```
