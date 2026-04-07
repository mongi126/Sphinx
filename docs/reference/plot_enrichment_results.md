# Generate publication-quality enrichment visualization plots

Generate publication-quality enrichment visualization plots

## Usage

``` r
plot_enrichment_results(
  cluster_enrich,
  top_n = 5,
  fdr_cutoff = 0.05,
  term_trunc_length = 40,
  base_font_size = 7,
  cluster_colors = NULL,
  plot_types = c("bar", "heatmap", "dot"),
  save_plot = FALSE,
  output_dir = "enrichment_plots"
)
```

## Arguments

- cluster_enrich:

  Enrichment results data frame

- top_n:

  Number of top terms to show per cluster (default: 5)

- fdr_cutoff:

  FDR cutoff for filtering (default: 0.05)

- term_trunc_length:

  Length to truncate term names (default: 40)

- base_font_size:

  Base font size for plots (default: 7)

- cluster_colors:

  Optional vector of colors for clusters

- plot_types:

  Types of plots to generate (default: all)

- save_plot:

  Whether to save plots (default: FALSE)

- output_dir:

  Output directory for saving (default: "enrichment_plots")

## Value

List of ggplot objects
