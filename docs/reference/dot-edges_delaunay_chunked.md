# Overlapping spatial-chunk Delaunay for large FOVs

Each grid core is triangulated with a spatial buffer. An edge is kept if
its midpoint falls in that core (half-open), so every edge has a unique
owner tile while the buffer supplies the neighbours needed for a correct
local DT.

## Usage

``` r
.edges_delaunay_chunked(df, max_edge_length = NULL, chunk_max_points = 250000L)
```
