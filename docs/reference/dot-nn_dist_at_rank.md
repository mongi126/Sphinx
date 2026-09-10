# Per-cell distance to the k-th nearest neighbor (excluding self)

Per-cell distance to the k-th nearest neighbor (excluding self)

## Usage

``` r
.nn_dist_at_rank(coords, k_nn = .default_nn_rank())
```

## Arguments

  - coords:
    
    numeric matrix with X/Y columns

  - k\_nn:
    
    neighbor rank (default: 10)

## Value

numeric vector of k\_nn-th NN distances per cell
