# Plot volcano plots for differential proteins across clusters

Plot volcano plots for differential proteins across clusters

## Usage

``` r
plot_volcano_cn_clusters(
  diff_results,
  fc_thresh = 0.25,
  p_thresh = 0.05,
  cn_cluster = NULL,
  y_max = NULL,
  cap_p = 1e-50,
  cap_jitter = 0.15,
  save_plot = FALSE,
  output_dir = "plots",
  filename = "volcano_plots"
)
```

## Arguments

- diff_results:

  Differential expression results (from perform_differential_expression)

- fc_thresh:

  Fold change threshold on Log2FC (default: 0.25)

- p_thresh:

  Adjusted p-value threshold (default: 0.05)

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

ggplot object
