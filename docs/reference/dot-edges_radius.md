# Build edges using radius method (exact ball graph)

Connects every pair of cells whose Euclidean distance is \<= `radius`
using `dbscan::frNN` as the primary search. Query points are split into
chunks and searched against the full coordinate set (`frNN(...,
query=)`). Chunks run in parallel via fork on Unix to lower peak memory
and wall time; neighbors match a single full `frNN` call after dropping
self-hits from the query interface. Undirected edges are emitted once
(`j > i`).

## Usage

``` r
.edges_radius(df, radius, n_cores = NULL)
```

## Arguments

  - df:
    
    data.table with spatial data

  - radius:
    
    search radius

  - n\_cores:
    
    parallel workers (NULL = auto)

## Value

data.table with edges
