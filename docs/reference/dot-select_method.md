# Select optimal method based on spatial metrics

Empirical decision tree for `method = "auto"`, calibrated on a
multi-cohort within-sample composite (PAS, SI\_expr, speed) across \~180
FOVs. Rules are applied in order; first match wins:

1.  Tiled / blocked-like FOV (ordered Clark-Evans R \> 1.15 and uniform
    density CV \< 0.3, n \< 5e4) -\> `window` with knn(k=10, mutual).

2.  Strong type clustering (`cat_morans_I >= 0.4`) **and** fairly
    uniform density (`cv_nn_dist < 0.3`) -\> `radius` at \~3x median
    1-NN.

3.  Continuous epithelium / few holes: large n (`>= 2e5`), heterogeneous
    density (`cv >= 0.6`), or space-filling irregular layout -\>
    `delaunay` with optional `max_edge_length = "auto"`.

4.  Otherwise baseline -\> `knn` with `n_neighbors = 10`,
    `require_mutual = TRUE`.

Scoring parallelism (`n_jobs`) is an evaluation concern, not encoded
here; builders still use `n_cores`.

## Usage

``` r
.select_method(m, n)
```

## Arguments

  - m:
    
    spatial metrics

  - n:
    
    number of cells

## Value

list with `method`, `params`, and `reason`
