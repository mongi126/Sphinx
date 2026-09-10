# Generate publication-quality enrichment visualization plots

Generate publication-quality enrichment visualization plots

## Usage

``` r
plot_enrichment_results(
  cluster_enrich,
  top_n = 5,
  fdr_cutoff = 0.05,
  term_trunc_length = 40,
  base_font_size = 12,
  cluster_colors = NULL,
  plot_types = c("bar", "heatmap", "dot"),
  save_plot = FALSE,
  output_dir = "enrichment_plots",
  width = 9,
  height = 10
)
```

## Arguments

  - cluster\_enrich:
    
    Enrichment results data frame

  - top\_n:
    
    Number of top terms to show per cluster (default: 5)

  - fdr\_cutoff:
    
    FDR cutoff for filtering (default: 0.05)

  - term\_trunc\_length:
    
    Length to truncate term names (default: 40)

  - base\_font\_size:
    
    Base font size for plots (default: 7)

  - cluster\_colors:
    
    Optional vector of colors for clusters

  - plot\_types:
    
    Types of plots to generate (default: all)

  - save\_plot:
    
    Whether to save plots (default: FALSE)

  - output\_dir:
    
    Output directory for saving (default: "enrichment\_plots")

  - width:
    
    Plot width in inches when saving (default: 9)

  - height:
    
    Plot height in inches when saving (default: 10)

## Value

List of ggplot objects

## Examples

``` r
if (FALSE) { # \dontrun{
# Requires enrichment results from perform_cluster_enrichment():
# plot_enrichment_results(enrich)
} # }
```
