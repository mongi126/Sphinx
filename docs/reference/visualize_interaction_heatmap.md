# Visualize cell-cell interaction matrix

Visualize cell-cell interaction matrix

## Usage

``` r
visualize_interaction_heatmap(
  interaction_matrix,
  transform = TRUE,
  color_palette = NULL,
  main = "Cell-Cell Interaction Frequency",
  display_numbers = FALSE,
  number_format = NULL,
  cellwidth = 14,
  cellheight = 14,
  fontsize = 12,
  save_path = NULL,
  width = 10,
  height = 8
)
```

## Arguments

  - interaction\_matrix:
    
    Interaction matrix from analyze\_spatial\_interactions() (raw
    counts, enrichment, or abundance-normalized scores)

  - transform:
    
    Apply log2(x+1) transformation (default: TRUE). Set FALSE for
    enrichment / abundance-normalized matrices.

  - color\_palette:
    
    Color palette function (default: soft sequential)

  - main:
    
    Heatmap title

  - display\_numbers:
    
    Whether to show numbers in cells (default: FALSE)

  - number\_format:
    
    sprintf format for cell labels (default: auto)

  - cellwidth:
    
    Cell width in points (default: 14)

  - cellheight:
    
    Cell height in points (default: 14)

  - fontsize:
    
    Base font size (default: 12)

  - save\_path:
    
    Output file path (optional)

  - width:
    
    Plot width in inches (default: 10)

  - height:
    
    Plot height in inches (default: 8)

## Value

pheatmap object and saves plot to file if save\_path provided

## Examples

``` r
# \donttest{
df <- prepare_data(Sphinx:::.sphinx_example_df(40))
edges <- build_spatial_network(df, method = "knn", n_neighbors = 5, verbose = FALSE)
intx <- analyze_spatial_interactions(df, edges)
visualize_interaction_heatmap(intx$interaction_matrix,
  save_path = tempfile(fileext = ".pdf"))
visualize_interaction_heatmap(intx$enrichment_matrix, transform = FALSE,
  main = "Abundance-normalized contact enrichment")
# }
```
