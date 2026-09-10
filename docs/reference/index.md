# Package index

## Data Preprocessing

<!-- end list -->

  - `spatial_colors` : Custom color palette for visualizations
  - `load_spatial_data()` : Load spatial data from various sources
  - `filter_data()` : Filter low-quality cells from spatial data
  - `process_data()` : Process spatial data
  - `plot_elbow()` : Generate elbow plot for dimensionality reduction
  - `extract_spatial_coordinates()` : Extract spatial coordinates from
    Seurat object
  - `visualize_results()` : Visualize clustering results

## Cell Annotation

<!-- end list -->

  - `celltype_colors` : Custom color palette for cell type
    visualizations
  - `find_top_markers()` : Identify top marker proteins for cell
    clusters
  - `plot_marker_violin()` : Visualize marker expression using violin
    plots
  - `plot_marker_heatmap()` : Create heatmap of cluster marker
    expression
  - `annotate_celltypes()` : Annotate cell types based on cluster IDs
  - `plot_annotated_umap()` : Visualize cell types in UMAP space
  - `plot_spatial_distribution()` : Visualize spatial distribution of
    cell types
  - `plot_spatial_markers()` : Visualize marker expression in spatial
    coordinates
  - `plot_umap_markers()` : Visualize marker expression in UMAP space

## Spatial Network Analysis

<!-- end list -->

  - `prepare_data()` : Prepare spatial data for network analysis
  - `calculate_optimal_radius()` : Calculate optimal spatial analysis
    parameters
  - `calculate_optimal_window_size()` : Calculate an adaptive window
    size from cell density
  - `calculate_celltype_distances()` : Calculate distances between cell
    types
  - `visualize_distance_heatmap()` : Visualize distance heatmap between
    cell types
  - `visualize_distance_parallel()` : Visualize distance relationships
    using parallel coordinates plot
  - `build_spatial_network()` : Construct spatial neighborhoods from
    coordinates with automatic method selection
  - `calculate_neighborhood_features()` : Calculate neighborhood
    composition features
  - `cluster_neighborhoods()` : Cluster neighborhoods using combined
    spatial and compositional features
  - `visualize_spatial_distribution()` : Visualize spatial cell type
    distribution
  - `calculate_cluster_composition()` : Calculate cluster composition
    metrics
  - `plot_composition_barplot()` : Visualize cluster composition as bar
    plot
  - `plot_composition_heatmap()` : Visualize cluster composition as
    heatmap
  - `calculate_neighborhood_purity()` : Calculate neighborhood purity
    using flexible neighbor definitions
  - `visualize_neighborhood_purity()` : Visualize neighborhood purity
  - `analyze_spatial_interactions()` : Analyze spatial interactions
    between cell types
  - `visualize_interaction_heatmap()` : Visualize cell-cell interaction
    matrix
  - `visualize_interaction_network()` : Visualize spatial interaction
    network using graph layout
  - `visualize_spatial_network()` : Visualize spatial network with
    flexible edge selection
  - `visualize_voronoi()` : Visualize spatial data using Voronoi
    diagrams
  - `get_color_palette()` : Generate color palette for visualizations
  - `assign_celltype_colors()` : Stable cell-type -\> color mapping
    (same label =\> same color across plots)

## Functional Analysis

<!-- end list -->

  - `prepare_protein_data()` : Prepare protein expression data for
    functional analysis
  - `perform_differential_expression()` : Perform differential
    expression analysis
  - `plot_volcano_all_clusters()` : Plot volcano plots for differential
    proteins across all clusters
  - `perform_cluster_enrichment()` : Perform cluster-specific protein
    enrichment analysis
  - `plot_enrichment_results()` : Generate publication-quality
    enrichment visualization plots
