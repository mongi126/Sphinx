# Perform differential expression analysis

Compares each neighborhood cluster to all other cells. For
CLR-normalized protein values, reports mean difference (Target -
Control) rather than log2 fold change to avoid redundant log transforms.

## Usage

``` r
perform_differential_expression(protein_df)
```

## Arguments

- protein_df:

  Data frame with protein expression and cluster information

## Value

Data frame with columns Mean_Target, Mean_Control, MeanDiff, p.value,
adj.p.value, and Significance
