# Resolve parallel worker count for network builders

Order: explicit `n_cores` \> env `SPHINX_N_CORES` \> min(detected cores,
8).

## Usage

``` r
.sphinx_default_cores(n_cores = NULL)
```

## Arguments

  - n\_cores:
    
    optional integer

## Value

integer \>= 1
