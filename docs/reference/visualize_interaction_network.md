# Visualize spatial interaction network using graph layout

Visualize spatial interaction network using graph layout

## Usage

``` r
visualize_interaction_network(
  network,
  node_size_range = c(5, 15),
  edge_size_range = c(0.5, 3),
  label_size = 4,
  show_labels = TRUE,
  max_nodes = 50,
  save_path = NULL,
  width = 12,
  height = 10,
  layout = "fr",
  base_size = 14
)
```

## Arguments

  - network:
    
    An igraph network object from analyze\_spatial\_interactions()

  - node\_size\_range:
    
    Range of node sizes (default: c(5, 15))

  - edge\_size\_range:
    
    Range of edge sizes (default: c(0.5, 3))

  - label\_size:
    
    Text label size (default: 4)

  - show\_labels:
    
    Whether to show node labels (default: TRUE)

  - max\_nodes:
    
    Maximum number of nodes to display (default: 50)

  - save\_path:
    
    Output file path (optional)

  - width:
    
    Plot width in inches (default: 12)

  - height:
    
    Plot height in inches (default: 10)

  - layout:
    
    Network layout algorithm (default: "fr")

  - base\_size:
    
    Base font size in points (default: 14; minimum 8)

## Value

ggraph object and saves plot to file if save\_path provided

## Examples

``` r
# \donttest{
df <- prepare_data(Sphinx:::.sphinx_example_df(40))
edges <- build_spatial_network(df, method = "knn", n_neighbors = 5, verbose = FALSE)
intx <- analyze_spatial_interactions(df, edges)
p <- visualize_interaction_network(intx$network, save_path = tempfile(fileext = ".pdf"))
class(p)
# }
```
