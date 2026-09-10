# Construct spatial neighborhoods from coordinates with automatic method selection

This function builds a cell-cell spatial network using one of four
methods ("radius", "knn", "delaunay", "window"). When `method = "auto"`,
it selects the method from within-sample spatial metrics using an
empirical decision tree (multi-cohort network evaluation; see
`.select_method()`): baseline knn(k=10, mutual); radius when types are
strongly clustered and density is fairly uniform; delaunay for
continuous / large / heterogeneous layouts (optional length gate);
window only for ordered tiled-like FOVs.

## Usage

``` r
build_spatial_network(
  df,
  method = "auto",
  k = NULL,
  radius = NULL,
  window_size = NULL,
  window_slide_step = NULL,
  max_edge_length = NULL,
  max_dist = NULL,
  edge_length_mult = 1,
  edge_length_quantile = 0.95,
  n_neighbors = 10L,
  max_degree = NULL,
  require_mutual = TRUE,
  max_edges = NULL,
  n_cores = NULL,
  celltype_col = "celltype",
  verbose = TRUE
)
```

## Arguments

  - df:
    
    data.frame/data.table with columns: X, Y, Cell\_ID; optional cell
    type column

  - method:
    
    one of c("auto","radius","knn","delaunay","window")

  - k:
    
    **deprecated** neighbor count. Use `n_neighbors`. Clustering `k`
    belongs to `cluster_neighborhoods()`, not this function.

  - radius:
    
    numeric, search radius for radius method (auto if NULL)

  - window\_size:
    
    numeric, sliding-window side length for window method (auto if
    NULL). Stereopy CCD: target 30-50 cells per window.

  - window\_slide\_step:
    
    numeric, stride of the sliding window. `NULL` (default) uses half of
    `window_size`, as in Stereopy CCD.

  - max\_edge\_length:
    
    distance gate. `NULL` / `Inf` (default): no gate for
    knn/delaunay/window. Numeric: keep edges with `dist <=` this value.
    `"auto"`: use `edge_length_quantile` of per-cell k-th NN distances
    (scaled by `edge_length_mult`). For `method = "radius"`, the search
    radius itself is the gate.

  - max\_dist:
    
    numeric absolute distance cap in the same units as X/Y. If set,
    overrides `max_edge_length`.

  - edge\_length\_mult:
    
    numeric multiplier for `max_edge_length = "auto"` (default 1.0).
    Ignored unless auto gating is requested.

  - edge\_length\_quantile:
    
    quantile of per-cell k-th NN distances when `max_edge_length =
    "auto"` (default 0.95).

  - n\_neighbors:
    
    integer number of spatial neighbors (knn/window search rank and
    post-filter degree cap). Default `10` for knn/window (aligned with
    `.default_nn_rank()` and the auto baseline). For `method = "radius"`
    or `"delaunay"`, omitted/`NULL` means no degree cap (`Inf`: keep the
    ball graph / triangulation). Alias of `max_degree`. This is **not**
    the clustering `k` used by `cluster_neighborhoods()`.

  - max\_degree:
    
    integer. Alias of `n_neighbors`. Used only when `n_neighbors` is not
    supplied. If both are set, `n_neighbors` wins. Set `Inf` to disable
    degree capping.

  - require\_mutual:
    
    logical. If TRUE, keep an edge only when both cells list each other
    among their `n_neighbors` nearest contacts. Default `TRUE` for
    knn/window (auto baseline). For `method = "radius"` defaults to
    `FALSE` unless set explicitly. Delaunay edges are already mutual;
    the flag only matters if a degree cap is also set.

  - max\_edges:
    
    deprecated; ignored. Edges are filtered by distance and per-cell
    degree only (no global random downsampling).

  - n\_cores:
    
    integer parallel workers for knn / radius construction. NULL uses
    `SPHINX_N_CORES` if set, otherwise up to 8 detected cores (fork on
    Unix). Also passed as `num.threads` to kNN search when using a
    single chunk.

  - celltype\_col:
    
    character, column with cell-type labels (optional)

  - verbose:
    
    logical, print decisions and key metrics

## Value

data.table with columns: from, to, dist (and optional context columns)

## Details

By default there is **no** distance gate for knn / delaunay / window
(same idea as leaving `max_dist` unset), except that `method = "auto"`
may recommend `max_edge_length = "auto"` for delaunay. Optionally set
`max_dist` (absolute), a numeric `max_edge_length`, or `max_edge_length
= "auto"` (upper quantile of per-cell k-th NN distances) to drop long
edges. Each cell can also be limited to its `n_neighbors` / `max_degree`
nearest contacts.

## Examples

``` r
df <- prepare_data(Sphinx:::.sphinx_example_df(40))
edges <- build_spatial_network(df, method = "knn", n_neighbors = 5, verbose = FALSE)
# Optional absolute distance gate (same units as X/Y):
edges2 <- build_spatial_network(
  df, method = "knn", n_neighbors = 10, max_dist = 40, verbose = FALSE
)
# Optional adaptive gate (quantile of k-th NN distances):
edges3 <- build_spatial_network(
  df, method = "knn", n_neighbors = 10, max_edge_length = "auto", verbose = FALSE
)
head(edges)
```
