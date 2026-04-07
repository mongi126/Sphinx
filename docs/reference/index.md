# Package index

## Data Preprocessing

- [`spatial_colors`](https://github.com/mongi126/Sphinx/reference/spatial_colors.md)
  : Custom color palette for visualizations
- [`load_spatial_data()`](https://github.com/mongi126/Sphinx/reference/load_spatial_data.md)
  : Load spatial data from various sources
- [`filter_data()`](https://github.com/mongi126/Sphinx/reference/filter_data.md)
  : Filter low-quality cells from spatial data
- [`process_data()`](https://github.com/mongi126/Sphinx/reference/process_data.md)
  : Process spatial data
- [`plot_elbow()`](https://github.com/mongi126/Sphinx/reference/plot_elbow.md)
  : Generate elbow plot for dimensionality reduction
- [`extract_spatial_coordinates()`](https://github.com/mongi126/Sphinx/reference/extract_spatial_coordinates.md)
  : Extract spatial coordinates from Seurat object
- [`visualize_results()`](https://github.com/mongi126/Sphinx/reference/visualize_results.md)
  : Visualize clustering results

## Cell Annotation

- [`celltype_colors`](https://github.com/mongi126/Sphinx/reference/celltype_colors.md)
  : Custom color palette for cell type visualizations
- [`find_top_markers()`](https://github.com/mongi126/Sphinx/reference/find_top_markers.md)
  : Identify top marker proteins for cell clusters
- [`plot_marker_violin()`](https://github.com/mongi126/Sphinx/reference/plot_marker_violin.md)
  : Visualize marker expression using violin plots
- [`plot_marker_heatmap()`](https://github.com/mongi126/Sphinx/reference/plot_marker_heatmap.md)
  : Create heatmap of cluster marker expression
- [`annotate_celltypes()`](https://github.com/mongi126/Sphinx/reference/annotate_celltypes.md)
  : Annotate cell types based on cluster IDs
- [`plot_annotated_umap()`](https://github.com/mongi126/Sphinx/reference/plot_annotated_umap.md)
  : Visualize cell types in UMAP space
- [`plot_spatial_distribution()`](https://github.com/mongi126/Sphinx/reference/plot_spatial_distribution.md)
  : Visualize spatial distribution of cell types
- [`plot_spatial_markers()`](https://github.com/mongi126/Sphinx/reference/plot_spatial_markers.md)
  : Visualize marker expression in spatial coordinates
- [`plot_umap_markers()`](https://github.com/mongi126/Sphinx/reference/plot_umap_markers.md)
  : Visualize marker expression in UMAP space

## Spatial Network Analysis

- [`prepare_data()`](https://github.com/mongi126/Sphinx/reference/prepare_data.md)
  : Prepare spatial data for network analysis
- [`calculate_optimal_radius()`](https://github.com/mongi126/Sphinx/reference/calculate_optimal_radius.md)
  : Calculate optimal spatial analysis parameters
- [`calculate_celltype_distances()`](https://github.com/mongi126/Sphinx/reference/calculate_celltype_distances.md)
  : Calculate distances between cell types
- [`visualize_distance_heatmap()`](https://github.com/mongi126/Sphinx/reference/visualize_distance_heatmap.md)
  : Visualize distance heatmap between cell types
- [`visualize_distance_parallel()`](https://github.com/mongi126/Sphinx/reference/visualize_distance_parallel.md)
  : Visualize distance relationships using parallel coordinates plot
- [`build_spatial_network()`](https://github.com/mongi126/Sphinx/reference/build_spatial_network.md)
  : Construct spatial neighborhoods from coordinates with automatic
  method selection
- [`calculate_neighborhood_features()`](https://github.com/mongi126/Sphinx/reference/calculate_neighborhood_features.md)
  : Calculate neighborhood composition features
- [`cluster_neighborhoods()`](https://github.com/mongi126/Sphinx/reference/cluster_neighborhoods.md)
  : Cluster neighborhoods using combined spatial and compositional
  features
- [`visualize_spatial_distribution()`](https://github.com/mongi126/Sphinx/reference/visualize_spatial_distribution.md)
  : Visualize spatial cell type distribution
- [`calculate_cluster_composition()`](https://github.com/mongi126/Sphinx/reference/calculate_cluster_composition.md)
  : Calculate cluster composition metrics
- [`plot_composition_barplot()`](https://github.com/mongi126/Sphinx/reference/plot_composition_barplot.md)
  : Visualize cluster composition as bar plot
- [`plot_composition_heatmap()`](https://github.com/mongi126/Sphinx/reference/plot_composition_heatmap.md)
  : Visualize cluster composition as heatmap
- [`calculate_neighborhood_purity()`](https://github.com/mongi126/Sphinx/reference/calculate_neighborhood_purity.md)
  : Calculate neighborhood purity using flexible neighbor definitions
- [`visualize_neighborhood_purity()`](https://github.com/mongi126/Sphinx/reference/visualize_neighborhood_purity.md)
  : Visualize neighborhood purity
- [`analyze_spatial_interactions()`](https://github.com/mongi126/Sphinx/reference/analyze_spatial_interactions.md)
  : Analyze spatial interactions between cell types
- [`visualize_interaction_heatmap()`](https://github.com/mongi126/Sphinx/reference/visualize_interaction_heatmap.md)
  : Visualize cell-cell interaction matrix
- [`visualize_interaction_network()`](https://github.com/mongi126/Sphinx/reference/visualize_interaction_network.md)
  : Visualize spatial interaction network using graph layout
- [`visualize_spatial_network()`](https://github.com/mongi126/Sphinx/reference/visualize_spatial_network.md)
  : Visualize spatial network with flexible edge selection
- [`visualize_voronoi()`](https://github.com/mongi126/Sphinx/reference/visualize_voronoi.md)
  : Visualize spatial data using Voronoi diagrams
- [`get_color_palette()`](https://github.com/mongi126/Sphinx/reference/get_color_palette.md)
  : Generate color palette for visualizations

## Functional Analysis

- [`prepare_protein_data()`](https://github.com/mongi126/Sphinx/reference/prepare_protein_data.md)
  : Prepare protein expression data for functional analysis
- [`perform_differential_expression()`](https://github.com/mongi126/Sphinx/reference/perform_differential_expression.md)
  : Perform differential expression analysis
- [`plot_volcano_cn_clusters()`](https://github.com/mongi126/Sphinx/reference/plot_volcano_cn_clusters.md)
  : Plot volcano plots for differential proteins across clusters
- [`perform_cluster_enrichment()`](https://github.com/mongi126/Sphinx/reference/perform_cluster_enrichment.md)
  : Perform cluster-specific protein enrichment analysis
- [`plot_enrichment_results()`](https://github.com/mongi126/Sphinx/reference/plot_enrichment_results.md)
  : Generate publication-quality enrichment visualization plots
