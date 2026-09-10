# Default maximum edge length from the k-th NN distance distribution

Uses the upper quantile of each cell's distance to its k-th nearest
neighbor. This matches the natural length scale of a kNN / local
triangulation graph: typical neighborhoods are kept, while only
sparse-region outliers are cut. (A cap of \~1.5xmedian(1-NN) is far too
tight once k \>= 5-10.)

## Usage

``` r
.default_max_edge_length(
  coords,
  k = .default_nn_rank(),
  quantile = 0.95,
  multiplier = 1
)
```

## Arguments

  - coords:
    
    numeric matrix with X/Y columns (or a single spacing scalar for
    backward-compatible callers - then treated as median 1-NN x 3)

  - k:
    
    neighbor rank used for the distance distribution (default 10)

  - quantile:
    
    upper quantile in (0, 1\] (default 0.95)

  - multiplier:
    
    optional scale on the quantile (default 1)

## Value

numeric max edge length
