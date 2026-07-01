# R/annotation.R

#' Custom color palette for cell type visualizations
#'
#' A predefined vector of 36 distinct colors for consistent cell type coloring
celltype_colors <- c(
  "#E5D2DD", "#53A85F", "#F1BB72", "#F3B1A0", "#D6E7A3", "#57C3F3",
  "#476D87", "#E95C59", "#E59CC4", "#AB3282", "#23452F", "#BD956A",
  "#8C549C", "#585658", "#9FA3A8", "#E0D4CA", "#5F3D69", "#C5DEBA",
  "#58A4C3", "#E4C755", "#F7F398", "#AA9A59", "#E63863", "#E39A35",
  "#C1E6F3", "#6778AE", "#91D0BE", "#B53E2B", "#712820", "#DCC1DD",
  "#CCE0F5", "#CCC9E6", "#625D9E", "#68A180", "#3A6963", "#968175"
)

#' Identify top marker proteins for cell clusters
#'
#' @param seurat_obj Seurat object containing clustered spatial data
#' @param save_path Output file path for marker results (default: "top5proteins.csv")
#' @param assay Assay name containing protein expression (default: "Spatial")
#' @param min.pct Minimum detection fraction threshold (default: 0.25)
#' @param logfc.threshold Minimum log-fold change threshold (default: 0.25)
#' @param only.pos Whether to return only positive markers (default: TRUE)
#' @return Data frame containing top 5 markers per cluster
#' @examples
#' \donttest{
#' obj <- Sphinx:::.sphinx_example_seurat(30)
#' markers <- find_top_markers(obj, save_path = tempfile(fileext = ".csv"),
#'   assay = "RNA", min.pct = 0.1, logfc.threshold = 0.1)
#' head(markers)
#' }
#' @export
#'
#' @description
#' Performs differential expression analysis to identify cluster biomarkers:
#' 1. Uses Seurat's FindAllMarkers with positive expression only
#' 2. Selects top 5 markers per cluster based on average log2FC
#' 3. Saves results to CSV file
find_top_markers <- function(seurat_obj,
                             save_path = "top5proteins.csv",
                             assay = "Spatial",
                             min.pct = 0.25,
                             logfc.threshold = 0.25,
                             only.pos = TRUE) {

  # Validate input object
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object")
  }

  # Check if clustering has been performed
  if (!"seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
    stop("No cluster annotations found. Run FindClusters() first.")
  }

  # Check if assay exists
  if (!assay %in% names(seurat_obj@assays)) {
    stop("Assay '", assay, "' not found in Seurat object")
  }

  # Identify cluster markers using Seurat's differential expression
  markers <- Seurat::FindAllMarkers(
    object = seurat_obj,
    assay = assay,
    only.pos = only.pos,           # Consider only positively expressed markers
    min.pct = min.pct,             # Minimum fraction of cells expressing marker
    logfc.threshold = logfc.threshold  # Minimum log-fold change
  )

  # Check if any markers were found
  if (nrow(markers) == 0) {
    warning("No significant markers found with current thresholds")
    return(data.frame())
  }

  # Select top 5 markers per cluster by expression fold change
  top5 <- markers %>%
    dplyr::group_by(cluster) %>%
    dplyr::top_n(n = 5, wt = avg_log2FC)
  utils::write.csv(top5, save_path, row.names = FALSE)

  # Save marker results to CSV
  utils::write.csv(top5, save_path, row.names = FALSE)
  message("Top markers saved to: ", save_path)

  return(top5)
}

#' Visualize marker expression using violin plots
#'
#' @param seurat_obj Seurat object with cluster annotations
#' @param markers Vector of marker proteins to visualize
#' @param group_by Metadata column for grouping (default: "seurat_clusters")
#' @param assay Assay containing expression data (default: "Spatial")
#' @param ncol Number of columns for multi-plot layout (default: 3)
#' @param save_path Output file path (default: "allmarker.pdf")
#' @param width Plot width in inches (default: 12)
#' @param height Plot height in inches (default: 30)
#' @return ggplot object containing violin plots
#' @examples
#' \donttest{
#' obj <- Sphinx:::.sphinx_example_seurat(25)
#' markers <- rownames(obj)[1:3]
#' p <- plot_marker_violin(obj, markers, assay = "RNA",
#'   save_path = tempfile(fileext = ".pdf"))
#' class(p)
#' }
#' @export
#'
#' @description
#' Generates violin plots showing expression distribution of marker proteins:
#' 1. Creates one plot per marker protein
#' 2. Groups cells by cluster identity
#' 3. Removes legend for cleaner multi-plot output
#' 4. Saves composite plot to PDF
plot_marker_violin <- function(seurat_obj,
                               markers,
                               group_by = "seurat_clusters",
                               assay = "Spatial",
                               ncol = 3,
                               save_path = "allmarker.pdf",
                               width = 12, height = 30) {

  # Validate marker input
  if (length(markers) == 0) {
    stop("No markers provided for visualization")
  }

  # Check if markers exist in the dataset
  missing_markers <- setdiff(markers, rownames(seurat_obj))
  if (length(missing_markers) > 0) {
    warning("The following markers not found in dataset: ",
            paste(missing_markers, collapse = ", "))
    markers <- base::intersect(markers, rownames(seurat_obj))
    if (length(markers) == 0) {
      stop("No valid markers found in dataset")
    }
  }

  # Create violin plots showing expression distributions
  p <- Seurat::VlnPlot(
    object = seurat_obj,
    features = markers,      # Proteins to visualize
    group.by = group_by,     # Grouping variable (clusters)
    assay = assay,           # Expression data source
    pt.size = 0,             # Hide individual points
    ncol = ncol              # Multi-plot layout
  ) +
    ggplot2::theme(legend.position = "none")  # Remove legend

  # Save composite plot to PDF
  ggplot2::ggsave(save_path, plot = p, width = width, height = height)
  message("Violin plot saved to: ", save_path)

  return(p)
}

#' Create heatmap of cluster marker expression
#'
#' @param seurat_obj Seurat object with cluster annotations
#' @param top_markers Data frame from find_top_markers
#' @param group_by Metadata column for grouping (default: "seurat_clusters")
#' @param assay Assay containing expression data (default: "Spatial")
#' @param save_path Output file path (default: "marker_heatmap.pdf")
#' @param width Plot width in inches (default: 12)
#' @param height Plot height in inches (default: 8)
#' @return Heatmap plot object
#' @examples
#' \donttest{
#' obj <- Sphinx:::.sphinx_example_seurat(25)
#' mk <- find_top_markers(obj, save_path = tempfile(fileext = ".csv"),
#'   assay = "RNA", min.pct = 0.1, logfc.threshold = 0.1)
#' if (nrow(mk) > 0) {
#'   plot_marker_heatmap(obj, mk, assay = "RNA",
#'     save_path = tempfile(fileext = ".pdf"))
#' }
#' }
#' @export
#'
#' @description
#' Generates a heatmap visualizing average marker expression per cluster:
#' 1. Extracts unique markers from top_markers data frame
#' 2. Computes average expression per cluster
#' 3. Applies row-wise Z-score normalization
#' 4. Uses Red-Blue color scheme for intuitive visualization
#' 5. Saves high-quality PDF output
plot_marker_heatmap <- function(seurat_obj,
                                top_markers,
                                group_by = "seurat_clusters",
                                assay = "Spatial",
                                save_path = "marker_heatmap.pdf",
                                width = 12, height = 8) {

  # Validate input
  if (nrow(top_markers) == 0) {
    stop("No marker data provided")
  }

  # Extract unique marker proteins
  top_genes <- unique(top_markers$gene)

  # Check if markers exist in dataset
  valid_genes <- base::intersect(top_genes, rownames(seurat_obj))
  if (length(valid_genes) == 0) {
    stop("None of the provided markers found in the dataset")
  }
  if (length(valid_genes) < length(top_genes)) {
    warning("Some markers not found in dataset: ",
            length(top_genes) - length(valid_genes), " markers omitted")
  }

  # Calculate average expression per cluster
  avg_expr <- Seurat::AverageExpression(
    object = seurat_obj,
    assays = assay,
    features = valid_genes,
    group.by = group_by
  )[[assay]]

  # Check if average expression calculation was successful
  if (nrow(avg_expr) == 0) {
    stop("Failed to calculate average expression")
  }

  # Create heatmap with hierarchical clustering
  heatmap_plot <- pheatmap::pheatmap(
    mat = avg_expr,
    scale = "row",  # Row-wise Z-score normalization
    clustering_method = "complete",  # Hierarchical clustering method
    color = grDevices::colorRampPalette(
      rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
    )(100),
    border_color = NA,  # Remove cell borders
    fontsize_row = 8,   # Adjust row label size
    filename = save_path,
    width = width,
    height = height
  )

  message("Marker heatmap saved to: ", save_path)
  return(heatmap_plot)
}

#' Annotate cell types based on cluster IDs
#'
#' @param seurat_obj Seurat object with cluster assignments
#' @param cluster_ids Vector of cluster IDs to annotate
#' @param celltype_labels Vector of cell type labels corresponding to cluster_ids
#' @param cluster_column Metadata column containing cluster IDs (default: "seurat_clusters")
#' @return Seurat object with added "celltype" metadata
#' @examples
#' \donttest{
#' obj <- Sphinx:::.sphinx_example_seurat(20)
#' cls <- levels(obj@meta.data$seurat_clusters)
#' obj <- annotate_celltypes(obj, cls, paste0("Type", seq_along(cls)))
#' table(obj@meta.data$celltype)
#' }
#' @export
#'
#' @description
#' Assigns biological cell type annotations to clusters:
#' 1. Validates equal length of cluster_ids and celltype_labels
#' 2. Creates cluster-to-celltype mapping
#' 3. Adds "celltype" column to object metadata
#' 4. Handles both numeric and character cluster IDs
annotate_celltypes <- function(seurat_obj,
                               cluster_ids,
                               celltype_labels,
                               cluster_column = "seurat_clusters") {
  # Validate input lengths
  if (length(cluster_ids) != length(celltype_labels)) {
    stop("cluster_ids and celltype_labels must have the same length")
  }

  # Validate cluster column exists
  if (!cluster_column %in% colnames(seurat_obj@meta.data)) {
    stop("Cluster column '", cluster_column, "' not found in metadata")
  }

  # Check if all cluster_ids exist in the data
  unique_clusters <- unique(seurat_obj@meta.data[[cluster_column]])
  missing_clusters <- setdiff(cluster_ids, unique_clusters)
  if (length(missing_clusters) > 0) {
    warning("The following cluster IDs not found in data: ",
            paste(missing_clusters, collapse = ", "))
  }

  # Create mapping vector
  cluster_mapping <- setNames(celltype_labels, cluster_ids)

  # Apply mapping to create celltype annotations
  current_clusters <- as.character(seurat_obj@meta.data[[cluster_column]])
  seurat_obj@meta.data$celltype <- cluster_mapping[current_clusters]

  # Report annotation statistics
  annotated_cells <- sum(!is.na(seurat_obj@meta.data$celltype))
  message("Cell types annotated: ",
          length(unique(na.omit(seurat_obj@meta.data$celltype))))
  message("Annotated cells: ", annotated_cells, "/", ncol(seurat_obj))

  return(seurat_obj)
}

#' Visualize cell types in UMAP space
#'
#' @param seurat_obj Seurat object with cell type annotations
#' @param save_path Output file path (default: "celltype_umap.pdf")
#' @param width Plot width in inches (default: 11)
#' @param height Plot height in inches (default: 8)
#' @return ggplot object containing annotated UMAP
#' @examples
#' \donttest{
#' obj <- Sphinx:::.sphinx_example_seurat(25)
#' cls <- levels(obj@meta.data$seurat_clusters)
#' obj <- annotate_celltypes(obj, cls, paste0("Type", seq_along(cls)))
#' p <- plot_annotated_umap(obj, save_path = tempfile(fileext = ".pdf"))
#' class(p)
#' }
#' @export
#'
#' @description
#' Creates UMAP visualization colored by annotated cell types:
#' 1. Uses custom color palette for distinct cell types
#' 2. Labels clusters with centered text boxes
#' 3. Maintains aspect ratio for proper spatial representation
#' 4. Removes default title for cleaner presentation
plot_annotated_umap <- function(seurat_obj,
                                save_path = "celltype_umap.pdf",
                                width = 11, height = 8) {

  # Validate cell type annotations
  if (!"celltype" %in% colnames(seurat_obj@meta.data)) {
    stop("Cell type annotations not found - run annotate_celltypes() first")
  }

  # Check if UMAP reduction exists
  if (!"umap" %in% names(seurat_obj@reductions)) {
    stop("UMAP reduction not found - run RunUMAP() first")
  }

  # Get unique cell types for color assignment
  unique_celltypes <- unique(na.omit(seurat_obj@meta.data$celltype))

  # Create UMAP plot with cell type coloring
  p <- Seurat::DimPlot(
    object = seurat_obj,
    reduction = "umap",
    cols = celltype_colors[1:length(unique_celltypes)],  # Use appropriate number of colors
    pt.size = 0.5,           # Point size optimization
    group.by = "celltype"
    ) +
    ggplot2::theme(aspect.ratio = 1) +  # Maintain aspect ratio
    ggplot2::theme(plot.title = ggplot2::element_blank())  # Remove default title

  # Save visualization to PDF
  ggplot2::ggsave(save_path, plot = p, width = width, height = height)
  message("Annotated UMAP saved to: ", save_path)

  return(p)
}

#' Visualize spatial distribution of cell types
#'
#' @param seurat_obj Seurat object with spatial coordinates
#' @param save_path Output file path (default: "celltype_spatial_plot.pdf")
#' @param width Plot width in inches (default: 10)
#' @param height Plot height in inches (default: 8)
#' @param point.size numeric, size of scatter points (default 0.5)
#' @return ggplot object containing spatial distribution plot
#' @examples
#' \donttest{
#' obj <- Sphinx:::.sphinx_example_seurat(25)
#' cls <- levels(obj@meta.data$seurat_clusters)
#' obj <- annotate_celltypes(obj, cls, paste0("Type", seq_along(cls)))
#' p <- plot_spatial_distribution(obj, save_path = tempfile(fileext = ".pdf"))
#' class(p)
#' }
#' @export
#'
#' @description
#' Generates spatial scatter plot showing cell type distribution:
#' 1. Requires prior extraction of spatial coordinates
#' 2. Uses custom color palette for cell type distinction
#' 3. Maintains spatial aspect ratio for accurate representation
#' 4. Optimizes point size and transparency for dense tissues
#' 5. Enhances legend readability for complex annotations
plot_spatial_distribution <- function(seurat_obj,
                                      save_path = "celltype_spatial_plot.pdf",
                                      width = 10, height = 8,
                                      point.size=0.5) {

  # Validate spatial coordinate existence
  if (!all(c("X", "Y") %in% colnames(seurat_obj@meta.data))) {
    stop("Spatial coordinates not found - run extract_spatial_coordinates() first")
  }

  # Validate cell type annotations
  if (!"celltype" %in% colnames(seurat_obj@meta.data)) {
    stop("Cell type annotations not found - run annotate_celltypes() first")
  }

  # Get unique cell types for color assignment
  unique_celltypes <- unique(na.omit(seurat_obj@meta.data$celltype))

  # Create spatial distribution plot
  p <- ggplot2::ggplot(seurat_obj@meta.data, ggplot2::aes(x = X, y = Y)) +
    ggplot2::geom_point(ggplot2::aes(color = celltype), size = point.size, alpha = 0.8) +
    ggplot2::theme_classic() +
    ggplot2::scale_color_manual(
      values = celltype_colors[1:length(unique_celltypes)],
      name = "Cell Type"
    ) +
    ggplot2::coord_fixed() +  # Maintain spatial aspect ratio
    ggplot2::labs(x = "X (um)", y = "Y (um)", title = "Spatial Cell Type Annotation") +
    ggplot2::theme(
      legend.text = ggplot2::element_text(size = 14),        # Increase legend text size
      legend.title = ggplot2::element_text(size = 16)         # Increase legend title size
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 2)))  # Increase legend point size

  # Save spatial plot to PDF
  ggplot2::ggsave(save_path, plot = p, width = width, height = height)
  message("Spatial distribution plot saved to: ", save_path)

  return(p)
}

#' Visualize marker expression in spatial coordinates
#'
#' @param obj Seurat object with spatial coordinates
#' @param markers Vector of marker genes/proteins to visualize (default: NULL)
#' @param assay Assay name containing expression data (default: "Akoya")
#' @param color_low Color for low expression values (default: "lightgrey")
#' @param color_high Color for high expression values (default: "red")
#' @param point_size Point size for spatial scatter plot (default: 1.0)
#' @param alpha Transparency level for points (default: 0.9)
#' @param ncol Number of columns for multi-plot layout (default: 3)
#' @param output_file Output file path (default: "spatial_markers.pdf")
#' @return Combined ggplot object with spatial marker expression plots
#' @examples
#' \dontrun{
#' # Requires Akoya/phenocycler image slot in the Seurat object:
#' # plot_spatial_markers(obj, markers = c("CD3", "CD8"))
#' }
#' @export
#'
#' @description
#' Generates spatial visualization of marker expression patterns:
#' 1. Extracts spatial coordinates from phenocycler image slot
#' 2. Validates marker presence in specified assay
#' 3. Creates individual spatial plots for each marker
#' 4. Combines plots using patchwork for multi-panel visualization
#' 5. Saves high-quality PDF output with automatic sizing
plot_spatial_markers <- function(obj,
                                 markers = NULL,
                                 assay = "Akoya",
                                 color_low = "lightgrey",
                                 color_high = "red",
                                 point_size = 1.0,
                                 alpha = 0.9,
                                 ncol = 3,
                                 output_file = "spatial_markers.pdf") {

  # Validate input Seurat object
  if (!inherits(obj, "Seurat")) stop("Input must be a Seurat object.")
  if (!assay %in% names(obj@assays)) stop("Assay not found in Seurat object.")
  if (is.null(obj@images$phenocycler)) stop("No 'phenocycler' image slot found.")

  message("Generating spatial marker plots...")

  # Extract spatial coordinates from phenocycler image
  coords <- obj@images$phenocycler@boundaries$centroids@coords
  if (is.null(coords) || nrow(coords) == 0)
    stop("No spatial coordinates found in obj@images$phenocycler@boundaries$centroids@coords")

  # Determine markers to visualize
  all_features <- obj@assays[[assay]]@meta.data$var.features
  if (is.null(markers)) {
    markers <- all_features
    message("No markers provided, using all ", length(markers), " markers.")
  } else {
    missing <- setdiff(markers, all_features)
    if (length(missing) > 0)
      warning("Some specified markers not found: ", paste(missing, collapse = ", "))
    markers <- intersect(markers, all_features)
  }
  if (length(markers) == 0) stop("No valid markers found for visualization.")

  # Helper function to create individual spatial plots
  plot_df <- function(df,
                      low_color = color_low,
                      high_color = color_high,
                      title = NULL,
                      point_size = 1.2) {
    ggplot2::ggplot(df, ggplot2::aes(x, y, color = expr)) +
      ggplot2::geom_point(size = point_size, alpha = alpha) +
      ggplot2::scale_color_gradient(low = low_color, high = high_color) +
      ggplot2::theme_void() +
      ggplot2::theme(
        legend.position = "right",
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold")
      ) +
      ggplot2::coord_fixed() +
      ggplot2::labs(title = title, color = "Expression")
  }

  # Generate individual plots for each marker
  plots <- list()
  for (marker in markers) {
    meta <- obj@assays[[assay]]@meta.data
    feat_idx <- which(meta$var.features == marker)
    if (length(feat_idx) == 0) {
      warning("Marker ", marker, " not found, skipping.")
      next
    }
    expr <- as.numeric(obj@assays[[assay]]@layers$data[feat_idx, ])
    if (length(expr) != nrow(coords)) {
      warning("Length mismatch for marker ", marker, ", skipping.")
      next
    }
    df <- data.frame(x = coords[, 1], y = coords[, 2], expr = expr)
    p <- plot_df(df, low_color = color_low, high_color = color_high,
                 title = marker, point_size = point_size)
    plots[[marker]] <- p
  }

  if (length(plots) == 0) stop("No valid plots generated.")

  # Combine all plots into multi-panel figure
  combined_plot <- patchwork::wrap_plots(plots, ncol = ncol) +
    patchwork::plot_annotation(
      title = paste0("Spatial Marker Expression (", assay, ")"),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "bold")
      )
    )

  # Determine output dimensions
  width <- 4 * ncol
  height <- 4 * ceiling(length(plots) / ncol)

  # Save combined plot to PDF
  ggplot2::ggsave(output_file, combined_plot, width = width, height = height)
  message("Combined spatial marker plot saved to: ", output_file)

  return(combined_plot)
}

#' Visualize marker expression in UMAP space
#'
#' @param seurat_obj Seurat object with UMAP reduction
#' @param markers Vector of marker genes/proteins to visualize (default: NULL)
#' @param assay Assay name containing expression data (default: NULL)
#' @param ncol Number of columns for multi-plot layout (default: 3)
#' @param point_size Point size for UMAP scatter plot (default: 1.2)
#' @param alpha Transparency level for points (default: 0.9)
#' @param color_low Color for low expression values (default: "lightgrey")
#' @param color_high Color for high expression values (default: "red")
#' @param output_file Output file path (default: "umap_markers_plot.pdf")
#' @param width Plot width in inches (default: NULL, auto-calculated)
#' @param height Plot height in inches (default: NULL, auto-calculated)
#' @return Combined ggplot object with UMAP marker expression plots
#' @examples
#' \donttest{
#' obj <- Sphinx:::.sphinx_example_seurat(25)
#' p <- plot_umap_markers(obj, markers = rownames(obj)[1:2],
#'   output_file = tempfile(fileext = ".pdf"))
#' class(p)
#' }
#' @export
#'
#' @description
#' Generates UMAP visualization of marker expression patterns:
#' 1. Validates UMAP reduction existence in Seurat object
#' 2. Determines appropriate assay for expression data
#' 3. Creates individual FeaturePlot for each marker
#' 4. Ensures numeric expression matrix for proper visualization
#' 5. Combines plots using patchwork for multi-panel output
#' 6. Saves high-quality PDF with automatic dimension calculation
plot_umap_markers <- function(seurat_obj,
                              markers = NULL,
                              assay = NULL,
                              ncol = 3,
                              point_size = 1.2,
                              alpha = 0.9,
                              color_low = "lightgrey",
                              color_high = "red",
                              output_file = "umap_markers_plot.pdf",
                              width = NULL,
                              height = NULL) {

  message("Generating UMAP marker expression visualization...")

  # ---- Check Seurat object ----
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object.")
  }

  # ---- Determine assay ----
  if (is.null(assay)) assay <- DefaultAssay(seurat_obj)
  if (!assay %in% names(seurat_obj@assays)) {
    stop("Assay '", assay, "' not found in Seurat object.")
  }

  # ---- Get all features ----
  all_features <- rownames(seurat_obj)
  if (is.null(markers)) {
    markers <- all_features
    message("No specific markers provided. Using all available features (", length(markers), ").")
  } else {
    missing <- setdiff(markers, all_features)
    if (length(missing) > 0) {
      warning("Some specified markers not found in dataset: ", paste(missing, collapse = ", "))
      markers <- intersect(markers, all_features)
    }
  }

  if (length(markers) == 0) stop("No valid markers found for visualization.")

  # ---- Ensure marker expression is numeric ----
  expr_mat <- GetAssayData(seurat_obj, assay = assay, slot = "data")
  expr_mat <- as.matrix(expr_mat)
  storage.mode(expr_mat) <- "numeric"

  # ---- Generate plots ----
  plots <- lapply(markers, function(marker) {
    if (!marker %in% rownames(expr_mat)) return(NULL)
    Seurat::FeaturePlot(
      object = seurat_obj,
      features = marker,
      cols = c(color_low, color_high),
      pt.size = point_size,
      order = TRUE
    ) +
      ggplot2::ggtitle(marker) +
      ggplot2::theme(
        legend.position = "right",
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 12)
      )
  })

  plots <- plots[!sapply(plots, is.null)]

  # ---- Combine all plots ----
  combined_plot <- patchwork::wrap_plots(plots, ncol = ncol) +
    patchwork::plot_annotation(
      title = "UMAP Marker Expression Map",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16)
      )
    )

  # ---- Determine figure size ----
  if (is.null(width)) width <- 4 * ncol
  if (is.null(height)) height <- 4 * ceiling(length(markers) / ncol)

  ggplot2::ggsave(output_file, combined_plot, width = width, height = height)
  message("Combined UMAP marker plot saved to: ", output_file)

  return(combined_plot)
}
