# Build edges using kNN method

Exact k nearest neighbors via `BiocNeighbors::findKNN` (self excluded).
Parallelism uses `num.threads` inside the search (fork-safe). For large
point sets, query rows are processed in serial chunks via `subset=` so
peak memory stays bounded; chunk results concatenate to the same
neighbors as a single full call (no downsampling).

## Usage

``` r
.edges_knn(df, k, max_edge_length = NULL, n_cores = NULL)
```

## Arguments

  - df:
    
    data.table with spatial data

  - k:
    
    number of neighbors

  - max\_edge\_length:
    
    optional maximum edge length filter

  - n\_cores:
    
    parallel workers / search threads (NULL = auto)

## Value

data.table with edges
