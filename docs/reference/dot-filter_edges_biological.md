# Filter edges by optional distance and degree rules

Rules (no global random downsampling):

1.  Distance (optional): if `max_edge_length` is finite, keep edges with
    dist \<= max\_edge\_length.

2.  Neighbor cap: each cell keeps at most `max_degree` / `n_neighbors`
    nearest contacts (default 10 for non-radius methods).

3.  Optional mutual filter: if `require_mutual`, keep an edge only when
    both endpoints nominate each other (reciprocal contact).

## Usage

``` r
.filter_edges_biological(
  edges,
  max_edge_length = Inf,
  max_degree = 10L,
  require_mutual = TRUE,
  verbose = FALSE
)
```

## Arguments

  - edges:
    
    data.table with from, to, dist

  - max\_edge\_length:
    
    numeric distance cap (Inf to skip)

  - max\_degree:
    
    integer per-cell neighbor cap (Inf to skip); aka n\_neighbors

  - require\_mutual:
    
    logical; mutual nearest neighbors if TRUE

  - verbose:
    
    logical

## Value

filtered data.table with from, to, dist
