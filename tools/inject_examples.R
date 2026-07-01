# Inject @examples blocks into R source files (one-time helper)
pkg_root <- "d:/Sphinx"

examples <- list(
  load_spatial_data = c(
    "#' @examples",
    "#' \\dontrun{",
    "#' # From RDS or CSV files on disk:",
    "#' # obj <- load_spatial_data(filename = \"sample.rds\")",
    "#' # obj <- load_spatial_data(counts_file = \"counts.csv\", coords_file = \"coords.csv\")",
    "#' }"
  ),
  filter_data = c(
    "#' @examples",
    "#' \\donttest{",
    "#' obj <- .sphinx_example_seurat(30)",
    "#' obj_filt <- filter_data(obj)",
    "#' ncol(obj_filt)",
    "#' }"
  ),
  process_data = c(
    "#' @examples",
    "#' \\donttest{",
    "#' obj <- .sphinx_example_seurat(30)",
    "#' obj <- process_data(obj, dims = 1:3)",
    "#' }",
    "#' \\dontrun{",
    "#' # Full-resolution clustering on real data:",
    "#' # obj <- process_data(obj, dims = 1:10, resolution = 0.5)",
    "#' }"
  ),
  plot_elbow = c(
    "#' @examples",
    "#' \\donttest{",
    "#' obj <- .sphinx_example_seurat(30)",
    "#' p <- plot_elbow(obj, save_path = tempfile(fileext = \".pdf\"))",
    "#' class(p)",
    "#' }"
  ),
  extract_spatial_coordinates = c(
    "#' @examples",
    "#' \\donttest{",
    "#' obj <- .sphinx_example_seurat(20)",
    "#' obj <- extract_spatial_coordinates(obj)",
    "#' head(obj@meta.data[, c(\"X\", \"Y\")])",
    "#' }"
  ),
  visualize_results = c(
    "#' @examples",
    "#' \\donttest{",
    "#' obj <- .sphinx_example_seurat(25)",
    "#' obj <- annotate_celltypes(obj,",
    "#'   cluster_ids = levels(obj@meta.data$seurat_clusters),",
    "#'   celltype_labels = paste0(\"Type\", seq_along(levels(obj@meta.data$seurat_clusters)))",
    "#' )",
    "#' td <- tempdir()",
    "#' plots <- visualize_results(obj, save_dir = td)",
    "#' names(plots)",
    "#' }"
  ),
  find_top_markers = c(
    "#' @examples",
    "#' \\donttest{",
    "#' obj <- .sphinx_example_seurat(30)",
    "#' markers <- find_top_markers(obj, save_path = tempfile(fileext = \".csv\"),",
    "#'   assay = \"RNA\", min.pct = 0.1, logfc.threshold = 0.1)",
    "#' head(markers)",
    "#' }"
  ),
  plot_marker_violin = c(
    "#' @examples",
    "#' \\donttest{",
    "#' obj <- .sphinx_example_seurat(25)",
    "#' markers <- rownames(obj)[1:3]",
    "#' p <- plot_marker_violin(obj, markers, save_path = tempfile(fileext = \".pdf\"))",
    "#' class(p)",
    "#' }"
  ),
  plot_marker_heatmap = c(
    "#' @examples",
    "#' \\donttest{",
    "#' obj <- .sphinx_example_seurat(25)",
    "#' mk <- find_top_markers(obj, save_path = tempfile(fileext = \".csv\"),",
    "#'   assay = \"RNA\", min.pct = 0.1, logfc.threshold = 0.1)",
    "#' if (nrow(mk) > 0) {",
    "#'   plot_marker_heatmap(obj, mk, save_path = tempfile(fileext = \".pdf\"))",
    "#' }",
    "#' }"
  ),
  annotate_celltypes = c(
    "#' @examples",
    "#' obj <- .sphinx_example_seurat(20)",
    "#' cls <- levels(obj@meta.data$seurat_clusters)",
    "#' obj <- annotate_celltypes(obj, cls, paste0(\"Type\", seq_along(cls)))",
    "#' table(obj@meta.data$celltype)"
  ),
  plot_annotated_umap = c(
    "#' @examples",
    "#' \\donttest{",
    "#' obj <- .sphinx_example_seurat(25)",
    "#' cls <- levels(obj@meta.data$seurat_clusters)",
    "#' obj <- annotate_celltypes(obj, cls, paste0(\"Type\", seq_along(cls)))",
    "#' p <- plot_annotated_umap(obj, save_path = tempfile(fileext = \".pdf\"))",
    "#' class(p)",
    "#' }"
  ),
  plot_spatial_distribution = c(
    "#' @examples",
    "#' \\donttest{",
    "#' obj <- .sphinx_example_seurat(25)",
    "#' cls <- levels(obj@meta.data$seurat_clusters)",
    "#' obj <- annotate_celltypes(obj, cls, paste0(\"Type\", seq_along(cls)))",
    "#' p <- plot_spatial_distribution(obj, save_path = tempfile(fileext = \".pdf\"))",
    "#' class(p)",
    "#' }"
  ),
  plot_spatial_markers = c(
    "#' @examples",
    "#' \\dontrun{",
    "#' # Requires Akoya/phenocycler image slot in the Seurat object:",
    "#' # plot_spatial_markers(obj, markers = c(\"CD3\", \"CD8\"))",
    "#' }"
  ),
  plot_umap_markers = c(
    "#' @examples",
    "#' \\donttest{",
    "#' obj <- .sphinx_example_seurat(25)",
    "#' p <- plot_umap_markers(obj, markers = rownames(obj)[1:2],",
    "#'   output_file = tempfile(fileext = \".pdf\"))",
    "#' class(p)",
    "#' }"
  ),
  prepare_protein_data = c(
    "#' @examples",
    "#' cluster_df <- .sphinx_example_df(10)",
    "#' cluster_df$Neighborhood_Cluster <- sample(1:2, 10, replace = TRUE)",
    "#' expr_df <- data.frame(PD.1 = rnorm(10), CD3 = rnorm(10),",
    "#'   row.names = cluster_df$Cell_ID)",
    "#' merged <- prepare_protein_data(cluster_df, expr_df)",
    "#' nrow(merged)"
  ),
  perform_differential_expression = c(
    "#' @examples",
    "#' df <- .sphinx_example_protein_df(40)",
    "#' res <- perform_differential_expression(df)",
    "#' head(res[, c(\"Protein\", \"Cluster\", \"MeanDiff\")])"
  ),
  plot_volcano_all_clusters = c(
    "#' @examples",
    "#' \\donttest{",
    "#' df <- .sphinx_example_protein_df(40)",
    "#' res <- perform_differential_expression(df)",
    "#' p <- plot_volcano_all_clusters(res)",
    "#' class(p)",
    "#' }"
  ),
  perform_cluster_enrichment = c(
    "#' @examples",
    "#' \\dontrun{",
    "#' # Requires clusterProfiler and org.Hs.eg.db:",
    "#' # res <- perform_differential_expression(.sphinx_example_protein_df(40))",
    "#' # enrich <- perform_cluster_enrichment(res)",
    "#' }"
  ),
  plot_enrichment_results = c(
    "#' @examples",
    "#' \\dontrun{",
    "#' # Requires enrichment results from perform_cluster_enrichment():",
    "#' # plot_enrichment_results(enrich)",
    "#' }"
  ),
  prepare_data = c(
    "#' @examples",
    "#' df <- .sphinx_example_df(30)",
    "#' df <- prepare_data(df)",
    "#' head(df[, .(Cell_ID, X, Y, celltype)])"
  ),
  calculate_optimal_radius = c(
    "#' @examples",
    "#' df <- prepare_data(.sphinx_example_df(50))",
    "#' params <- calculate_optimal_radius(df)",
    "#' names(params)"
  ),
  calculate_optimal_window_size = c(
    "#' @examples",
    "#' df <- prepare_data(.sphinx_example_df(50))",
    "#' ws <- calculate_optimal_window_size(df)",
    "#' ws$window_size"
  ),
  calculate_celltype_distances = c(
    "#' @examples",
    "#' df <- prepare_data(.sphinx_example_df(40))",
    "#' dist_res <- calculate_celltype_distances(df, celltype_col = \"celltype\")",
    "#' dim(dist_res$distance_matrix)"
  ),
  build_spatial_network = c(
    "#' @examples",
    "#' df <- prepare_data(.sphinx_example_df(40))",
    "#' edges <- build_spatial_network(df, method = \"knn\", k = 5, verbose = FALSE)",
    "#' head(edges)"
  ),
  calculate_neighborhood_features = c(
    "#' @examples",
    "#' df <- prepare_data(.sphinx_example_df(40))",
    "#' edges <- build_spatial_network(df, method = \"knn\", k = 5, verbose = FALSE)",
    "#' feat <- calculate_neighborhood_features(df, edges)",
    "#' names(feat)"
  ),
  cluster_neighborhoods = c(
    "#' @examples",
    "#' \\donttest{",
    "#' df <- prepare_data(.sphinx_example_df(40))",
    "#' edges <- build_spatial_network(df, method = \"knn\", k = 5, verbose = FALSE)",
    "#' feat <- calculate_neighborhood_features(df, edges)",
    "#' cl <- cluster_neighborhoods(feat, edges, method = \"kmeans\", k = 3)",
    "#' \"Neighborhood_Cluster\" %in% names(cl)",
    "#' }"
  ),
  calculate_neighborhood_purity = c(
    "#' @examples",
    "#' df <- prepare_data(.sphinx_example_df(40))",
    "#' out <- calculate_neighborhood_purity(df, method = \"knn\", k = 5, verbose = FALSE)",
    "#' summary(out$Neighborhood_Purity)"
  ),
  analyze_spatial_interactions = c(
    "#' @examples",
    "#' df <- prepare_data(.sphinx_example_df(40))",
    "#' edges <- build_spatial_network(df, method = \"knn\", k = 5, verbose = FALSE)",
    "#' intx <- analyze_spatial_interactions(df, edges)",
    "#' dim(intx$interaction_matrix)"
  ),
  get_color_palette = c(
    "#' @examples",
    "#' get_color_palette(5)"
  ),
  visualize_spatial_distribution = c(
    "#' @examples",
    "#' \\donttest{",
    "#' df <- prepare_data(.sphinx_example_df(40))",
    "#' p <- visualize_spatial_distribution(df, save_path = tempfile(fileext = \".pdf\"))",
    "#' class(p)",
    "#' }"
  ),
  visualize_distance_heatmap = c(
    "#' @examples",
    "#' \\donttest{",
    "#' df <- prepare_data(.sphinx_example_df(40))",
    "#' dist_res <- calculate_celltype_distances(df, celltype_col = \"celltype\")",
    "#' p <- visualize_distance_heatmap(dist_res, save_path = tempfile(fileext = \".pdf\"))",
    "#' class(p)",
    "#' }"
  ),
  visualize_distance_parallel = c(
    "#' @examples",
    "#' \\donttest{",
    "#' df <- prepare_data(.sphinx_example_df(40))",
    "#' dist_res <- calculate_celltype_distances(df, celltype_col = \"celltype\")",
    "#' p <- visualize_distance_parallel(dist_res, save_path = tempfile(fileext = \".pdf\"))",
    "#' class(p)",
    "#' }"
  ),
  visualize_neighborhood_purity = c(
    "#' @examples",
    "#' \\donttest{",
    "#' df <- prepare_data(.sphinx_example_df(40))",
    "#' df <- calculate_neighborhood_purity(df, method = \"knn\", k = 5, verbose = FALSE)",
    "#' p <- visualize_neighborhood_purity(df, save_path = tempfile(fileext = \".pdf\"))",
    "#' class(p)",
    "#' }"
  ),
  calculate_cluster_composition = c(
    "#' @examples",
    "#' df <- prepare_data(.sphinx_example_df(40))",
    "#' df$Neighborhood_Cluster <- sample(1:3, nrow(df), replace = TRUE)",
    "#' comp <- calculate_cluster_composition(df)",
    "#' head(comp)"
  ),
  plot_composition_heatmap = c(
    "#' @examples",
    "#' \\donttest{",
    "#' df <- prepare_data(.sphinx_example_df(40))",
    "#' df$Neighborhood_Cluster <- sample(1:3, nrow(df), replace = TRUE)",
    "#' comp <- calculate_cluster_composition(df)",
    "#' plot_composition_heatmap(comp, save_path = tempfile(fileext = \".pdf\"))",
    "#' }"
  ),
  plot_composition_barplot = c(
    "#' @examples",
    "#' \\donttest{",
    "#' df <- prepare_data(.sphinx_example_df(40))",
    "#' df$Neighborhood_Cluster <- sample(1:3, nrow(df), replace = TRUE)",
    "#' comp <- calculate_cluster_composition(df)",
    "#' p <- plot_composition_barplot(comp, save_path = tempfile(fileext = \".pdf\"))",
    "#' class(p)",
    "#' }"
  ),
  visualize_interaction_heatmap = c(
    "#' @examples",
    "#' \\donttest{",
    "#' df <- prepare_data(.sphinx_example_df(40))",
    "#' edges <- build_spatial_network(df, method = \"knn\", k = 5, verbose = FALSE)",
    "#' intx <- analyze_spatial_interactions(df, edges)",
    "#' visualize_interaction_heatmap(intx$interaction_matrix,",
    "#'   save_path = tempfile(fileext = \".pdf\"))",
    "#' }"
  ),
  visualize_spatial_network = c(
    "#' @examples",
    "#' \\donttest{",
    "#' df <- prepare_data(.sphinx_example_df(40))",
    "#' edges <- build_spatial_network(df, method = \"knn\", k = 5, verbose = FALSE)",
    "#' p <- visualize_spatial_network(df, edges, save_path = tempfile(fileext = \".pdf\"))",
    "#' class(p)",
    "#' }"
  ),
  visualize_interaction_network = c(
    "#' @examples",
    "#' \\donttest{",
    "#' df <- prepare_data(.sphinx_example_df(40))",
    "#' edges <- build_spatial_network(df, method = \"knn\", k = 5, verbose = FALSE)",
    "#' intx <- analyze_spatial_interactions(df, edges)",
    "#' p <- visualize_interaction_network(intx$network, save_path = tempfile(fileext = \".pdf\"))",
    "#' class(p)",
    "#' }"
  ),
  visualize_voronoi = c(
    "#' @examples",
    "#' \\donttest{",
    "#' df <- prepare_data(.sphinx_example_df(40))",
    "#' p <- visualize_voronoi(df, save_path = tempfile(fileext = \".pdf\"))",
    "#' class(p)",
    "#' }"
  )
)

inject_examples <- function(file) {
  lines <- readLines(file, warn = FALSE)
  fn_lines <- grep("^[a-zA-Z._][a-zA-Z0-9._]*\\s*<-\\s*function", lines)
  for (idx in sort(fn_lines, decreasing = TRUE)) {
    fn <- sub("\\s*<-.*", "", lines[idx])
    if (!fn %in% names(examples)) next
    block_start <- idx - 1L
    while (block_start > 0L && grepl("^#'", lines[block_start])) {
      block_start <- block_start - 1L
    }
    block_start <- block_start + 1L
    if (block_start >= idx) next
    block <- lines[block_start:(idx - 1L)]
    if (any(grepl("^#' @examples", block))) next
    export_idx <- which(grepl("^#' @export\\s*$", block))
    if (length(export_idx) == 0) next
    insert_at <- block_start + export_idx[1] - 1L
    lines <- c(
      lines[1:(insert_at - 1L)],
      examples[[fn]],
      lines[insert_at:length(lines)]
    )
  }
  writeLines(lines, file, useBytes = TRUE)
}

for (f in list.files(file.path(pkg_root, "R"), pattern = "\\.R$", full.names = TRUE)) {
  if (basename(f) %in% c("example-utils.R")) next
  inject_examples(f)
}
cat("Examples injected into", length(list.files(file.path(pkg_root, "R"), pattern = "\\.R$")) - 1L, "files.\n")
