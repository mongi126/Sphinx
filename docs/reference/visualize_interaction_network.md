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
  layout = "fr"
)
```

## Arguments

- network:

  An igraph network object from analyze_spatial_interactions()

- node_size_range:

  Range of node sizes (default: c(5, 15))

- edge_size_range:

  Range of edge sizes (default: c(0.5, 3))

- label_size:

  Text label size (default: 4)

- show_labels:

  Whether to show node labels (default: TRUE)

- max_nodes:

  Maximum number of nodes to display (default: 50)

- save_path:

  Output file path (optional)

- width:

  Plot width in inches (default: 12)

- height:

  Plot height in inches (default: 10)

- layout:

  Network layout algorithm (default: "fr")

## Value

ggraph object and saves plot to file if save_path provided
