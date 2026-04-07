# R/preprocessing.R

#' Custom color palette for visualizations
#'
#' A predefined vector of 36 distinct colors for consistent plotting
spatial_colors <- c(
  "#E5D2DD", "#53A85F", "#F1BB72", "#F3B1A0", "#D6E7A3", "#57C3F3",
  "#476D87", "#E95C59", "#E59CC4", "#AB3282", "#23452F", "#BD956A",
  "#8C549C", "#585658", "#9FA3A8", "#E0D4CA", "#5F3D69", "#C5DEBA",
  "#58A4C3", "#E4C755", "#F7F398", "#AA9A59", "#E63863", "#E39A35",
  "#C1E6F3", "#6778AE", "#91D0BE", "#B53E2B", "#712820", "#DCC1DD",
  "#CCE0F5", "#CCC9E6", "#625D9E", "#68A180", "#3A6963", "#968175"
)

#' Load spatial data from various sources
#'
#' @param filename Path to input data file (RDS, CSV counts, or CSV coordinates)
#' @param counts_file Path to counts CSV file (if using separate counts and coordinates)
#' @param coords_file Path to coordinates CSV file (if using separate counts and coordinates)
#' @param type Data source type: "auto", "rds", "csv", or "akoya" (default: "auto")
#' @param akoya_type Type of Akoya data: "qupath" or other formats (default: "qupath")
#' @param fov Field of view specification (default: "phenocycler")
#' @param ... Additional parameters passed to loading functions
#' @return Seurat object containing spatial data
#' @export
#'
#' @description
#' Loads spatial data from either:
#' 1. An existing Seurat object (RDS file)
#' 2. Counts and coordinates CSV files
#' 3. Akoya platform data files
load_spatial_data <- function(filename = NULL,
                              counts_file = NULL,
                              coords_file = NULL,
                              type = c("auto", "rds", "csv", "akoya"),
                              akoya_type = "qupath",
                              fov = "phenocycler",
                              ...) {

  type <- match.arg(type, several.ok = FALSE)

  ## 1. Auto-detect file type if set to "auto" --------------------------------
  if (type == "auto") {
    if (!is.null(filename) && file.exists(filename)) {
      if (grepl("\\.rds$", filename, ignore.case = TRUE)) {
        type <- "rds"
      } else if (grepl("\\.csv$", filename, ignore.case = TRUE)) {
        type <- "akoya"   # Default single CSV treated as Akoya output
      } else {
        stop("Cannot auto-detect file type from 'filename' extension.")
      }
    } else if (!is.null(counts_file) && !is.null(coords_file)) {
      type <- "csv"
    } else {
      stop("Please provide either an RDS file, Akoya csv, or both counts and coordinates CSV files.")
    }
  }

  ## 2. Process data based on detected type -----------------------------------
  if (type == "rds") {
    if (is.null(filename) || !file.exists(filename))
      stop("RDS file not found.")
    obj <- readRDS(filename)
    message("Loaded Seurat object from RDS file. Cell count: ", ncol(obj))
    return(obj)
  }

  if (type == "akoya") {
    if (is.null(filename) || !file.exists(filename))
      stop("Akoya csv file not found.")
    obj <- LoadAkoya(filename = filename,
                                type = akoya_type,
                                fov = fov)
    message("Loaded Akoya data into Seurat object. Cell count: ", ncol(obj))
    return(obj)
  }

  if (type == "csv") {
    if (is.null(counts_file) || is.null(coords_file))
      stop("Both 'counts_file' and 'coords_file' must be provided for type='csv'.")
    if (!file.exists(counts_file) || !file.exists(coords_file))
      stop("One or both input CSV files do not exist.")

    counts <- read.csv(counts_file, row.names = 1, check.names = FALSE)
    coords <- read.csv(coords_file, row.names = 1, check.names = FALSE)

    # Check if coordinates have required columns
    if (!all(c("X", "Y") %in% colnames(coords))) {
      if (ncol(coords) >= 2) {
        # Assume first two columns are coordinates
        colnames(coords)[1:2] <- c("X", "Y")
        message("Assuming first two columns of coordinates file are X and Y coordinates")
      } else {
        stop("Coordinates file must contain at least 2 columns for X and Y coordinates")
      }
    }

    common_cells <- base::intersect(colnames(counts), rownames(coords))
    if (length(common_cells) == 0)
      stop("No common cell names found between counts and coordinates files.")

    counts <- counts[, common_cells, drop = FALSE]
    coords <- coords[common_cells, , drop = FALSE]

    obj <- Seurat::CreateSeuratObject(counts = counts, assay = "Spatial")

    # Add coordinates to metadata
    obj$X <- coords$X
    obj$Y <- coords$Y

    # Create spatial reduction object
    spatial_coords <- as.matrix(coords[, c("X", "Y")])
    colnames(spatial_coords) <- c("X", "Y")
    obj[["spatial"]] <- Seurat::CreateDimReducObject(
      embeddings = spatial_coords,
      assay = "Spatial",
      key = "SPATIAL_"
    )

    message("Created Seurat object from CSV files. Cell count: ", ncol(obj))
    return(obj)
  }

  stop("Unknown 'type' value.")
}

#' Filter low-quality cells from spatial data
#'
#' @param obj Seurat object containing spatial data
#' @param nCount_mad_threshold MAD threshold for total count filtering (default: 3)
#' @param nFeature_quantile_threshold Quantile threshold for detected features (default: 0.05)
#' @return Filtered Seurat object
#' @export
#'
#' @description
#' Performs two-step quality control filtering:
#' 1. Removes cells with outlier total counts (MAD-based)
#' 2. Removes cells with low detected feature counts (quantile-based)
filter_data <- function(obj,
                        nCount_mad_threshold = 3,
                        nFeature_quantile_threshold = 0.05) {

  # Validate input object type
  if (!inherits(obj, "Seurat")) {
    stop("Input must be a Seurat object")
  }

  # Determine assay name (default to first assay)
  assay_name <- Seurat::DefaultAssay(obj)

  # Calculate dynamic thresholds using Median Absolute Deviation
  nCount_name <- paste0("nCount_", assay_name)
  nFeature_name <- paste0("nFeature_", assay_name)

  # Check if required metadata columns exist
  if (!nCount_name %in% colnames(obj@meta.data)) {
    stop("Column ", nCount_name, " not found in object metadata")
  }
  if (!nFeature_name %in% colnames(obj@meta.data)) {
    stop("Column ", nFeature_name, " not found in object metadata")
  }

  median_count <- median(obj[[nCount_name]][, 1])
  mad_count <- mad(obj[[nCount_name]][, 1], constant = 1)
  low_bound <- median_count - nCount_mad_threshold * mad_count
  high_bound <- median_count + nCount_mad_threshold * mad_count

  # Calculate feature quantile threshold
  feature_threshold <- quantile(obj[[nFeature_name]][, 1], nFeature_quantile_threshold)

  # Apply filtering using subset
  cells_keep <- obj[[nCount_name]][, 1] > low_bound &
    obj[[nCount_name]][, 1] < high_bound &
    obj[[nFeature_name]][, 1] > feature_threshold

  # Remove NA values
  cells_keep[is.na(cells_keep)] <- FALSE

  obj_filtered <- obj[, cells_keep]

  # Report post-filtering cell count
  message("Initial cell count: ", ncol(obj))
  message("Filtered cell count: ", ncol(obj_filtered))
  message("Cells removed: ", ncol(obj) - ncol(obj_filtered))

  return(obj_filtered)
}

#' Process spatial data
#'
#' @param obj Seurat object containing filtered spatial data
#' @param normalization.method Normalization method (default: "CLR")
#' @param margin Margin for normalization (1=features, 2=cells) (default: 2)
#' @param dims Dimensions for reduction (default: 1:10)
#' @param resolution Clustering resolution (default: 0.5)
#' @return Processed Seurat object with dimensionality reduction and clustering
#' @export
#'
#' @description
#' Performs comprehensive data processing pipeline:
#' 1. Normalization (CLR recommended for proteomics)
#' 2. Feature selection (uses all detected features)
#' 3. Data scaling
#' 4. PCA dimensionality reduction
#' 5. Nearest-neighbor graph construction
#' 6. Cluster identification
#' 7. UMAP visualization
process_data <- function(obj,
                         normalization.method = "CLR",
                         margin = 2,
                         dims = 1:10,
                         resolution = 0.5) {

  # Determine assay name
  assay_name <- Seurat::DefaultAssay(obj)

  # Fix meta.features if they don't match data features
  if (!identical(rownames(obj[[assay_name]]@meta.features), rownames(obj))) {
    message("Fixing meta.features to match data features...")
    obj[[assay_name]]@meta.features <- data.frame(
      row.names = rownames(obj),
      name = rownames(obj)
    )
  }

  # Data normalization - CLR recommended for protein data
  tryCatch({
    obj <- Seurat::NormalizeData(
      obj,
      normalization.method = normalization.method,
      margin = margin,
      assay = assay_name,
      verbose = FALSE
    )
    message("Data normalization completed successfully")
  }, error = function(e) {
    message("NormalizeData failed with error: ", e$message)
    message("Using manual CLR normalization instead...")

    # Manual CLR normalization implementation
    counts_data <- Seurat::GetAssayData(obj, assay = assay_name, slot = "counts")

    # CLR normalization: log(1 + x / (exp(mean(log(1 + x)))))
    clr_normalize <- function(x) {
      x_nz <- x[x > 0]  # Consider only non-zero values
      if (length(x_nz) > 0) {
        mean_log <- mean(log1p(x_nz))
        log1p(x / exp(mean_log))
      } else {
        x  # If all values are zero, keep as is
      }
    }

    # Apply CLR normalization
    if (margin == 2) {  # Normalize by cells
      clr_data <- apply(counts_data, 2, clr_normalize)
    } else {  # Normalize by features
      clr_data <- apply(counts_data, 1, clr_normalize)
      clr_data <- t(clr_data)  # Transpose to maintain original dimensions
    }

    # Add normalized data to object
    obj <- Seurat::SetAssayData(
      obj,
      assay = assay_name,
      slot = "data",
      new.data = as(clr_data, "dgCMatrix")
    )
    message("Manual CLR normalization completed")
  })

  # Feature selection - use all detected features
  Seurat::VariableFeatures(obj) <- rownames(obj)

  # Feature scaling
  obj <- Seurat::ScaleData(obj)

  # Dimensionality reduction using PCA
  obj <- Seurat::RunPCA(obj)

  # Cell clustering workflow
  obj <- Seurat::FindNeighbors(obj, dims = dims)    # Build kNN graph
  obj <- Seurat::FindClusters(obj, resolution = resolution)  # Identify clusters

  # Non-linear dimensionality reduction for visualization
  obj <- Seurat::RunUMAP(obj, dims = dims)

  message("Data processing completed successfully")
  return(obj)
}

#' Generate elbow plot for dimensionality reduction
#'
#' @param obj Processed Seurat object
#' @param save_path Output file path (default: "elbow_plot.pdf")
#' @return ggplot object containing elbow plot
#' @export
#'
#' @description
#' Creates elbow plot to determine optimal number of PCA dimensions.
#' Saves plot to specified path and returns ggplot object.
plot_elbow <- function(obj, save_path = "elbow_plot.pdf") {
  # Generate elbow plot of PCA standard deviations
  p <- Seurat::ElbowPlot(obj)

  # Save plot to PDF
  suppressMessages(ggplot2::ggsave(save_path, p, width = 8, height = 6))
  message("Elbow plot saved to: ", save_path)

  return(p)
}

#' Extract spatial coordinates from Seurat object
#'
#' @param obj Seurat object with spatial data
#' @return Seurat object with coordinates added to metadata (if not already present)
#' @export
#'
#' @description
#' Extracts spatial coordinates from reduction slot and adds them to object metadata.
#' Verifies successful extraction by printing coordinate sample.
extract_spatial_coordinates <- function(obj) {
  # Check if spatial reduction exists
  if ("spatial" %in% names(obj@reductions)) {
    # Extract coordinates from spatial reduction
    coords <- Seurat::Embeddings(obj, reduction = "spatial")

    # Add coordinates to metadata if not already present
    if (!all(c("X", "Y") %in% colnames(obj@meta.data))) {
      obj$X <- coords[, 1]
      obj$Y <- coords[, 2]
      message("Spatial coordinates added to metadata from reduction.")
    } else {
      message("Spatial coordinates already present in metadata.")
    }
  } else if (all(c("X", "Y") %in% colnames(obj@meta.data))) {
    # Create spatial reduction from metadata coordinates
    spatial_coords <- as.matrix(obj@meta.data[, c("X", "Y")])
    colnames(spatial_coords) <- c("SPATIAL_1", "SPATIAL_2")

    obj[["spatial"]] <- Seurat::CreateDimReducObject(
      embeddings = spatial_coords,
      assay = Seurat::DefaultAssay(obj),
      key = "SPATIAL_"
    )
    message("Spatial reduction created from metadata coordinates.")
  } else {
    warning("No spatial coordinates found in reductions or metadata.")
    return(obj)
  }

  # Verification message with coordinate sample
  message("Sample coordinates:")
  print(utils::head(obj@meta.data[, c("X", "Y")]))

  return(obj)
}

#' Visualize clustering results
#'
#' @param obj Processed Seurat object with spatial coordinates
#' @param save_dir Output directory (default: "./")
#' @return List containing UMAP and spatial plot objects
#' @export
#'
#' @description
#' Generates and saves two key visualizations:
#' 1. UMAP projection showing cell clusters
#' 2. Spatial scatter plot showing cluster distribution in tissue context
visualize_results <- function(obj, save_dir = "./") {

  # Check if seurat_clusters exists
  if (!"seurat_clusters" %in% colnames(obj@meta.data)) {
    stop("No clustering results found. Run FindClusters() first.")
  }

  # Check if coordinates exist
  if (!all(c("X", "Y") %in% colnames(obj@meta.data))) {
    stop("Spatial coordinates not found in metadata. Run extract_spatial_coordinates() first.")
  }

  # Create UMAP cluster visualization
  p1 <- Seurat::DimPlot(
    obj,
    reduction = "umap",
    cols = spatial_colors,
    pt.size = 0.5,
    label = TRUE,
    label.box = TRUE
  ) +
    ggplot2::theme(aspect.ratio = 1) +
    ggplot2::ggtitle("Cluster UMAP")

  # Create spatial cluster distribution plot
  p2 <- ggplot2::ggplot(obj@meta.data, ggplot2::aes(x = X, y = Y)) +
    ggplot2::geom_point(ggplot2::aes(color = seurat_clusters), size = 0.5, alpha = 0.8) +
    ggplot2::theme_classic() +
    ggplot2::scale_color_manual(values = spatial_colors, name = "Cluster") +
    ggplot2::coord_fixed() +
    ggplot2::labs(x = "X (um)", y = "Y (um)", title = "Spatial Clusters") +
    ggplot2::theme(
      legend.text = ggplot2::element_text(size = 12, face = "bold"),
      legend.title = ggplot2::element_text(size = 13)
    )

  # Ensure output directory exists
  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
    message("Created output directory: ", save_dir)
  }

  # Save visualizations
  ggplot2::ggsave(file.path(save_dir, "cluster_umap.pdf"), p1, width = 10, height = 8)
  ggplot2::ggsave(file.path(save_dir, "spatial_clusters.pdf"), p2, width = 10, height = 8)

  message("Visualizations saved to: ", save_dir)

  return(list(umap_plot = p1, spatial_plot = p2))
}
