# Perform differential expression analysis

Compares each neighborhood cluster to all other cells. For
CLR-normalized protein values, reports mean difference (Target -
Control) rather than log2 fold change to avoid redundant log transforms.

## Usage

``` r
perform_differential_expression(
  protein_df,
  test_level = c("spatial_block", "cell"),
  block_size = NULL,
  min_blocks = 3L
)
```

## Arguments

  - protein\_df:
    
    Data frame with protein expression and cluster information

  - test\_level:
    
    `"spatial_block"` (default) or `"cell"` (legacy; not recommended due
    to pseudo-replication)

  - block\_size:
    
    numeric grid side length for spatial blocks. NULL = auto.

  - min\_blocks:
    
    minimum number of Target and Control blocks required (default: 3)

## Value

Data frame with columns Mean\_Target, Mean\_Control, MeanDiff, p.value,
adj.p.value, Significance, and `n_target` / `n_control` (cells or
blocks)

## Details

By default, tests are run on **spatial-block means** rather than
individual cells. Treating every cell as an independent observation
overstates degrees of freedom because nearby cells are spatially
correlated (pseudo-replication). Aggregating to tissue blocks reduces
that inflation.

## Examples

``` r
df <- Sphinx:::.sphinx_example_protein_df(40)
res <- perform_differential_expression(df)
head(res[, c("Protein", "Cluster", "MeanDiff")])
```
