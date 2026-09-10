# Build edges with a Stereopy-style sliding window + local kNN

Windows of side `tile` (`window_size`) slide over the field with stride
`sliding_step` (default half the window). Inside each window a kNN graph
is built; overlapping windows are merged (undirected, shortest edge
kept).

## Usage

``` r
.edges_window(df, k, tile, sliding_step = NULL, max_edge_length = NULL)
```

## Arguments

  - df:
    
    data.table with spatial data

  - k:
    
    number of neighbors inside each window

  - tile:
    
    window side length (Stereopy `d`)

  - sliding\_step:
    
    window stride (Stereopy `s`; default `tile / 2`)

  - max\_edge\_length:
    
    optional maximum edge length filter

## Value

data.table with edges
