# Build edges using Delaunay triangulation

For large point sets, `deldir` (Fortran) fails with "long vectors are
not supported in .Fortran". In that regime we tile the FOV into
overlapping spatial chunks, triangulate each chunk, and merge undirected
edges (deduplicated).

## Usage

``` r
.edges_delaunay(df, max_edge_length = NULL, k = 6L, chunk_max_points = 250000L)
```

## Arguments

  - df:
    
    data.table with spatial data

  - max\_edge\_length:
    
    optional maximum edge length filter

  - k:
    
    neighbor count used only if `deldir` is missing (kNN fallback)

  - chunk\_max\_points:
    
    max points per chunk before tiling (default 2.5e5)

## Value

data.table with edges
