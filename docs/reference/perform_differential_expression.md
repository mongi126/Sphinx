# Perform differential expression analysis

Perform differential expression analysis

Plot volcano plots for differential proteins across all clusters

## Usage

``` r
perform_differential_expression(protein_df)

perform_differential_expression(protein_df)
```

## Arguments

- protein_df:

  Data frame with protein expression and cluster information

- diff_results:

  Differential expression results

- fc_thresh:

  Fold change threshold (default: 0.25)

- p_thresh:

  P-value threshold (default: 0.05)

- cn_cluster:

  Optional vector of clusters to plot (default: NULL = all)

- y_max:

  Optional maximum for y-axis (logP) (default: NULL = auto)

- save_plot:

  Whether to save the plot (default: FALSE)

- output_dir:

  Output directory for saving (default: "plots")

- filename:

  Output filename (default: "volcano_plots")

## Value

Data frame with differential expression results

ggplot object
