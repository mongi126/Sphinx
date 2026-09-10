# Categorical Moran's I via proportion-weighted binary Moran's I per type

Categorical Moran's I via proportion-weighted binary Moran's I per type

## Usage

``` r
.cat_morans_I_knn(coords, types, k)
```

## Arguments

  - coords:
    
    numeric matrix of X/Y coordinates

  - types:
    
    character/factor vector of category labels

  - k:
    
    number of nearest neighbors (excluding self)

## Value

scalar Moran's I aggregated across categories, or NA
