# Select optimal method based on spatial metrics

Decision tree (when called from `method = "auto"`):

1.  n \> 50,000 -\> window

2.  cat_morans_I \>= 0.4 -\> radius (strong type clustering)

3.  cv_nn_dist \< 0.3 -\> knn (CSR-like) or radius (mild deviation from
    CSR within the same band)

4.  otherwise -\> delaunay (fallback)

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

list with `method` (character) and `params` (list of recommended
parameters)
