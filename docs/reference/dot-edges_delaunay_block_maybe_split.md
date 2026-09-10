# Triangulate one expanded block; keep edges whose midpoint is in core

Triangulate one expanded block; keep edges whose midpoint is in core

## Usage

``` r
.edges_delaunay_block_maybe_split(
  block,
  in_core_mid,
  max_edge_length = NULL,
  chunk_max_points = 250000L,
  depth = 0L
)
```
