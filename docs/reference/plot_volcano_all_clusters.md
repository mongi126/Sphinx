# Plot volcano plots for differential proteins across all clusters

Y-axis `-log10(adj.p)` is capped at `y_cap` (default 50). Points
exceeding the cap are drawn as upward triangles on the truncation line.

## Usage

``` r
plot_volcano_all_clusters(
  diff_results,
  diff_thresh = 0.25,
  p_thresh = 0.05,
  y_cap = 50,
  base_size = 14,
  ncol = 3,
  save_plot = FALSE,
  output_dir = "plots",
  filename = "volcano_plots",
  width = 12,
  height = 8
)
```

## Arguments

  - diff\_results:
    
    Differential expression results

  - diff\_thresh:
    
    Mean difference threshold, Target - Control (default: 0.25)

  - p\_thresh:
    
    P-value threshold (default: 0.05)

  - y\_cap:
    
    Truncation for -log10(adj.p.value) (default: 50)

  - base\_size:
    
    Base font size (default: 14)

  - ncol:
    
    Number of facet columns (default: 3)

  - save\_plot:
    
    Whether to save the plot (default: FALSE)

  - output\_dir:
    
    Output directory for saving (default: "plots")

  - filename:
    
    Output filename (default: "volcano\_plots")

  - width:
    
    Plot width in inches (default: 12)

  - height:
    
    Plot height in inches (default: 8)

## Value

ggplot object

## Examples

``` r
# \donttest{
df <- Sphinx:::.sphinx_example_protein_df(40)
res <- perform_differential_expression(df)
p <- plot_volcano_all_clusters(res)
class(p)
# }
```
