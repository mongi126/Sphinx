# R/spatial_network.R

#' Prepare spatial data for network analysis
#'
#' @param df Data frame containing spatial data
#' @param cell_id_col Column name for cell IDs (default: NULL, auto-detect)
#' @param x_col Column name for X coordinates (default: "X")
#' @param y_col Column name for Y coordinates (default: "Y")
#' @param celltype_col Column name for cell types (default: "celltype")
#' @return data.table with standardized structure for spatial analysis
#' @export
#'
#' @description
#' Prepares spatial data for network analysis by:
#' 1. Converting to data.table if needed
#' 2. Ensuring Cell_ID column exists
#' 3. Converting coordinates to numeric
#' 4. Removing invalid coordinates and missing cell types
prepare_data <- function(df,
                         cell_id_col = NULL,
                         x_col = "X",
                         y_col = "Y",
                         celltype_col = "celltype") {
  # Ensure data.table package is loaded
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Please install data.table package: install.packages('data.table')")
  }

  # Convert to data.table if not already
  if (!data.table::is.data.table(df)) {
    df <- data.table::as.data.table(df)
  } else {
    df <- data.table::copy(df)  # Avoid modifying original data
  }

  # 1. Standardize Cell_ID column name
  if (!is.null(cell_id_col)) {
    if (!cell_id_col %in% names(df)) {
      stop("Specified cell_id_col '", cell_id_col, "' does not exist in the data")
    }
    # Only rename if different from target name
    if (cell_id_col != "Cell_ID") {
      data.table::setnames(df, cell_id_col, "Cell_ID")
    }
  } else if ("V1" %in% names(df)) {
    data.table::setnames(df, "V1", "Cell_ID")
  } else if (!"Cell_ID" %in% names(df)) {
    if ("cell" %in% names(df)) {
      data.table::setnames(df, "cell", "Cell_ID")
    } else {
      # No suitable column found, create new one named Cell_ID
      df[, Cell_ID := .I]
    }
  }

  # 2. Standardize coordinate column names
  if (!"X" %in% names(df) && x_col == "X" && "x" %in% names(df)) {
    x_col <- "x"
  }
  if (!"Y" %in% names(df) && y_col == "Y" && "y" %in% names(df)) {
    y_col <- "y"
  }
  if (x_col != "X" && x_col %in% names(df)) {
    if ("X" %in% names(df)) {
      df[, X := NULL]
    }
    data.table::setnames(df, x_col, "X")
  } else if (!"X" %in% names(df)) {
    stop("X coordinate column not found")
  }

  if (y_col != "Y" && y_col %in% names(df)) {
    if ("Y" %in% names(df)) {
      df[, Y := NULL]  # Remove existing Y to avoid Y.1
    }
    data.table::setnames(df, y_col, "Y")
  } else if (!"Y" %in% names(df)) {
    stop("Y coordinate column not found")
  }

  df[, c("X", "Y") := lapply(.SD, function(z) as.numeric(as.character(z))),
     .SDcols = c("X", "Y")]

  # 3. Filter invalid rows
  df <- df[is.finite(X) & is.finite(Y)]

  # Check if celltype column exists and filter NAs
  if (celltype_col %in% names(df)) {
    df <- df[!is.na(get(celltype_col))]
  }

  return(df)
}

#' Calculate optimal spatial analysis parameters
#'
#' @param df Spatial data with coordinates
#' @param sample_size Maximum cells to sample for efficiency (default: 1000)
#' @param multiplier Factor for recommended radius (default: 2.0)
#' @param x_col Column name for X coordinates (default: "X")
#' @param y_col Column name for Y coordinates (default: "Y")
#' @param k_nn Nearest-neighbor rank used for distance summary (default: 10)
#' @return List with distance statistics and recommended parameters
#' @export
calculate_optimal_radius <- function(df,
                                     sample_size = 1000,
                                     multiplier = 2.0,
                                     x_col = "X",
                                     y_col = "Y",
                                     k_nn = 10L) {
  if (!requireNamespace("RANN", quietly = TRUE)) {
    stop("Package 'RANN' is required. Please install it.")
  }

  # Ensure we get a data.frame, not a vector
  coords <- data.frame(
    x = as.numeric(as.character(df[[x_col]])),
    y = as.numeric(as.character(df[[y_col]]))
  )

  # Remove NAs
  coords <- coords[stats::complete.cases(coords), ]

  if (nrow(coords) < 2) {
    stop("Not enough valid coordinates after removing NAs. Got only ",
         nrow(coords), " valid rows.")
  }

  coords <- as.matrix(coords)

  k_nn <- as.integer(k_nn)
  if (k_nn < 1L) stop("k_nn must be at least 1")
  k_query <- min(k_nn + 1L, nrow(coords))
  nn <- RANN::nn2(coords, k = k_query)
  nn_dists <- nn$nn.dists[, k_query]

  # Sample to calculate all pairwise distances
  set.seed(123)
  sample_size <- min(sample_size, nrow(coords))
  sample_idx <- sample(nrow(coords), sample_size)
  all_dists <- as.vector(stats::dist(coords[sample_idx, ]))
  all_dists <- all_dists[all_dists > 0]

  # Recommended radius: k_nn-th NN distance (90th percentile) scaled by multiplier
  recommended_radius <- as.numeric(stats::quantile(nn_dists, 0.9, na.rm = TRUE)) * multiplier

  return(list(
    recommended_radius = recommended_radius,
    min_radius = min(nn_dists, na.rm = TRUE),
    max_radius = max(nn_dists, na.rm = TRUE),
    min_dists = nn_dists,
    k_nn = k_query - 1L,
    all_dists = all_dists
  ))
}

#' Calculate optimal window size (Stereopy CCD-style)
#'
#' @param df Spatial data with coordinates
#' @param x_col Column name for X coordinates (default: "X")
#' @param y_col Column name for Y coordinates (default: "Y")
#' @param min_cells Target minimum mean cells per window (default: 30)
#' @param max_cells Target maximum mean cells per window (default: 50)
#' @param range_divisor Divisor for initial window size from spatial range (default: 100)
#' @param max_iter Maximum adjustment iterations (default: 100)
#' @return List with window_size, sliding_step, mean_cells, x_range, y_range
#' @export
calculate_optimal_window_size <- function(df,
                                          x_col = "X",
                                          y_col = "Y",
                                          min_cells = 30L,
                                          max_cells = 50L,
                                          range_divisor = 100L,
                                          max_iter = 100L) {
  x <- as.numeric(as.character(df[[x_col]]))
  y <- as.numeric(as.character(df[[y_col]]))
  ok <- stats::complete.cases(x, y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 1L) {
    stop("No valid coordinates after removing NAs.")
  }
  res <- .calc_stereopy_window_size(
    x, y,
    min_cells = min_cells,
    max_cells = max_cells,
    range_divisor = range_divisor,
    max_iter = max_iter
  )
  c(
    res,
    list(
      x_range = diff(range(x, na.rm = TRUE)),
      y_range = diff(range(y, na.rm = TRUE))
    )
  )
}

#' Calculate distances between cell types
#'
#' @param df Spatial data with coordinates and cell types
#' @param x_col Column name for X coordinates (default: "X")
#' @param y_col Column name for Y coordinates (default: "Y")
#' @param celltype_col Column name for cell types (default: "annotation")
#' @return List with distance matrix and distributions
#' @export
calculate_celltype_distances <- function(df,
                                         x_col = "X",
                                         y_col = "Y",
                                         celltype_col = "annotation") {
  if (!requireNamespace("RANN", quietly = TRUE)) {
    stop("Package 'RANN' is required.")
  }

  # Convert to data.table for faster processing
  if (!data.table::is.data.table(df)) {
    df <- data.table::as.data.table(df)
  }

  # Initialize result storage
  result <- list()

  # Calculate global distances
  coords <- as.matrix(df[, c(x_col, y_col), with = FALSE])
  nn_all <- RANN::nn2(coords, k = 2)

  result$global <- list(
    min_dists = nn_all$nn.dists[, 2],
    all_dists = as.vector(stats::dist(coords[sample(nrow(coords), min(10000, nrow(coords))), ])),
    recommended_radius = stats::quantile(nn_all$nn.dists[, 2], 0.9),
    min_radius = min(nn_all$nn.dists[, 2]),
    max_radius = max(nn_all$nn.dists[, 2])
  )

  # Get all cell types
  cell_types <- unique(df[[celltype_col]])

  # Initialize distance matrix
  dist_matrix <- matrix(NA, nrow = length(cell_types), ncol = length(cell_types))
  rownames(dist_matrix) <- colnames(dist_matrix) <- cell_types

  # Store all distance distributions
  dist_distributions <- list()

  # Calculate pairwise distances between cell types
  for (i in seq_along(cell_types)) {
    type1 <- cell_types[i]
    idx1  <- which(df[[celltype_col]] == type1)

    for (j in seq_along(cell_types)) {
      type2 <- cell_types[j]
      idx2  <- which(df[[celltype_col]] == type2)

      if (length(idx1) == 0 || length(idx2) == 0) next  # Skip if either is empty

      if (type1 == type2) {
        if (length(idx1) > 1) {
          sub_coords <- coords[idx1, , drop = FALSE]
          nn <- RANN::nn2(sub_coords, k = 2)
          distances <- nn$nn.dists[, 2]
        } else {
          distances <- NA_real_
        }
      } else {
        query_coords  <- coords[idx1, , drop = FALSE]
        target_coords <- coords[idx2, , drop = FALSE]
        k_safe <- min(1, nrow(target_coords))
        nn <- RANN::nn2(target_coords, query = query_coords, k = k_safe)
        distances <- nn$nn.dists[, k_safe]
      }

      dist_matrix[i, j] <- ifelse(length(distances) > 0, mean(distances, na.rm = TRUE), NA_real_)
      dist_distributions[[paste(type1, type2, sep = "_")]] <- distances
    }
  }

  result$distance_matrix <- dist_matrix
  result$distance_distributions <- dist_distributions

  return(result)
}

#' Construct spatial neighborhoods from coordinates with automatic method selection
#'
#' This function builds a cell-cell spatial network using one of four methods
#' ("radius", "knn", "delaunay", "window"). When `method = "auto"`, it
#' selects the method based on robust spatial metrics computed from the data.
#'
#' @param df data.frame/data.table with columns: X, Y, Cell_ID; optional cell type column
#' @param method one of c("auto","radius","knn","delaunay","window")
#' @param k integer, neighbors for knn/window methods (default chosen automatically)
#' @param radius numeric, search radius for radius method (auto if NULL)
#' @param window_size numeric, grid side length for window method (auto if NULL)
#' @param max_edges integer, hard cap on number of returned edges (may be downsampled)
#' @param celltype_col character, column with cell-type labels (optional)
#' @param verbose logical, print decisions and key metrics
#' @return data.table with columns: from, to, dist (and optional context columns)
#' @export
build_spatial_network <- function(
    df,
    method = "auto",
    k = NULL,
    radius = NULL,
    window_size = NULL,
    max_edges = 1e6,
    celltype_col = "celltype",
    verbose = TRUE
) {
  # ---- Check required packages ----
  req <- c("data.table","RANN")
  miss <- req[!vapply(req, function(p) requireNamespace(p, quietly = TRUE), logical(1))]
  if (length(miss)) stop("Missing packages: ", paste(miss, collapse = ", "))

  if (!data.table::is.data.table(df)) {
    df <- data.table::as.data.table(df)
  }

  if (!all(c("X", "Y", "Cell_ID") %in% names(df))) {
    df <- prepare_data(df, celltype_col = celltype_col)
  }
  if (!all(c("X", "Y", "Cell_ID") %in% names(df))) {
    stop("`df` must contain coordinates (X/Y or x/y) and cell IDs (Cell_ID or cell)")
  }

  if (nrow(df) < 2) {
    return(data.table::data.table(from=character(), to=character(), dist=numeric()))
  }

  # ---- Calculate spatial metrics (guide auto choices + sensible defaults) ----
  metrics <- .calc_spatial_metrics(df, celltype_col = celltype_col)

  # ---- Auto method selection (returns method + recommended params) ----
  chosen <- method
  auto_params <- list()
  if (method == "auto") {
    sel <- .select_method(metrics, nrow(df))
    chosen <- sel$method
    auto_params <- sel$params
  }

  # Resolve parameters: user override > auto recommendation > generic default
  if (is.null(k)) {
    k <- if (!is.null(auto_params$k)) {
      auto_params$k
    } else {
      max(4, min(30, round(8 + log10(max(10, nrow(df))) * 4)))
    }
  }
  if (is.null(radius)) {
    radius <- if (!is.null(auto_params$radius)) {
      auto_params$radius
    } else {
      .default_search_radius(metrics$mean_nn_dist)
    }
  }
  if (is.null(window_size)) {
    if (!is.null(auto_params$window_size)) {
      window_size <- auto_params$window_size
      window_slide_step <- auto_params$window_slide_step
    } else {
      opt_win <- .calc_stereopy_window_size(df$X, df$Y)
      window_size <- opt_win$window_size
      window_slide_step <- opt_win$sliding_step
    }
  } else {
    window_slide_step <- max(1L, as.integer(window_size) %/% 2L)
  }
  if (is.null(window_slide_step)) {
    window_slide_step <- max(1L, as.integer(window_size) %/% 2L)
  }
  max_edge_length <- auto_params$max_edge_length

  if (verbose) {
    msg <- paste0(
      "Method=", chosen,
      " | n=", nrow(df),
      " | radius=", round(radius,3),
      " | window_size=", round(window_size,3),
      " | meanNN=", round(metrics$mean_nn_dist,3),
      " | cvNN=", round(metrics$cv_nn_dist,3),
      " | CE_R=", round(metrics$clark_evans_R,3),
      " | CatMoranI=", round(metrics$cat_morans_I,3),
      " | Het(GS)=", round(metrics$het_gini_simpson,3)
    )
    message(msg)
  }

  # ---- Build base edges ----
  edges <- switch(
    chosen,
    knn      = .edges_knn(df, k = k),
    radius   = .edges_radius(df, radius = radius, k_cap = max(20, k * 3)),
    delaunay = .edges_delaunay(df, max_edge_length = max_edge_length),
    window   = .edges_window(df, k = min(k, 10), tile = window_size),
    stop("Unknown method: ", chosen)
  )

  if (nrow(edges) == 0L) {
    warning("No edges produced; try increasing k/radius/window_size.")
    return(edges)
  }

  # ---- Optional context (cheap) ----
  if (!is.null(celltype_col) && celltype_col %in% names(df)) {
    df_small <- df[, .(Cell_ID, celltype = get(celltype_col))]
    edges[df_small, on = .(from = Cell_ID), from_type := i.celltype]
    edges[df_small, on = .(to   = Cell_ID), to_type   := i.celltype]
    edges[, same_type := (from_type == to_type)]
  }

  # ---- Cap edges (distance-weighted sampling) ----
  if (nrow(edges) > max_edges) {
    if (verbose) message("Downsampling edges: ", nrow(edges), " -> ", max_edges)
    w <- 1 / (edges$dist + .Machine$double.eps)
    w <- w / sum(w)
    idx <- sample.int(nrow(edges), size = max_edges, replace = FALSE, prob = w)
    edges <- edges[idx]
  }

  data.table::setattr(edges, "metrics", metrics)
  data.table::setattr(edges, "method", chosen)
  data.table::setattr(edges, "parameters", list(
    k = k,
    radius = radius,
    window_size = window_size,
    window_slide_step = window_slide_step,
    max_edge_length = max_edge_length
  ))
  data.table::setkey(edges, from, to)

  return(edges)
}

# =====================
# Internal helper functions
# =====================

#' Default nearest-neighbor rank for local spacing summaries
#' @keywords internal
.default_nn_rank <- function() 10L

#' Per-cell distance to the k-th nearest neighbor (excluding self)
#' @param coords numeric matrix with X/Y columns
#' @param k_nn neighbor rank (default: 10)
#' @return numeric vector of k_nn-th NN distances per cell
#' @keywords internal
.nn_dist_at_rank <- function(coords, k_nn = .default_nn_rank()) {
  n <- nrow(coords)
  if (n < 2L) return(rep(NA_real_, n))
  k_nn <- as.integer(k_nn)
  if (k_nn < 1L) stop("k_nn must be at least 1")
  k_query <- min(k_nn + 1L, n)
  nn <- RANN::nn2(coords, coords, k = k_query)
  nn$nn.dists[, k_query]
}

#' Default search radius from mean NN distance
#' @param mean_nn_dist mean k-th nearest-neighbor distance
#' @param multiplier scale factor (default: 2.0)
#' @return numeric search radius
#' @keywords internal
.default_search_radius <- function(mean_nn_dist, multiplier = 2.0) {
  mean_nn_dist * multiplier
}

#' Round to the nearest even integer (minimum 2)
#' @keywords internal
.round_to_even <- function(x) {
  r <- 2 * round(x / 2)
  max(2L, as.integer(r))
}

#' Number of grid windows covering the spatial extent (Stereopy-style binning)
#' @keywords internal
.n_grid_windows <- function(x, y, win_size) {
  xmin <- min(x, na.rm = TRUE)
  xmax <- max(x, na.rm = TRUE)
  ymin <- min(y, na.rm = TRUE)
  ymax <- max(y, na.rm = TRUE)
  nx <- floor(xmax / win_size) - floor(xmin / win_size) + 1L
  ny <- floor(ymax / win_size) - floor(ymin / win_size) + 1L
  max(1L, nx) * max(1L, ny)
}

#' Mean cells per grid window over the full tissue grid
#' @keywords internal
.mean_cells_per_window_grid <- function(x, y, win_size) {
  length(x) / .n_grid_windows(x, y, win_size)
}

#' Optimal window size via Stereopy CCD iterative adjustment
#' @keywords internal
.calc_stereopy_window_size <- function(x,
                                       y,
                                       min_cells = 30L,
                                       max_cells = 50L,
                                       range_divisor = 100L,
                                       max_iter = 100L) {
  x_range <- diff(range(x, na.rm = TRUE))
  y_range <- diff(range(y, na.rm = TRUE))
  if (!is.finite(x_range) || !is.finite(y_range) || x_range <= 0 || y_range <= 0) {
    return(list(window_size = 2L, sliding_step = 1L, mean_cells = NA_real_))
  }

  win_size <- .round_to_even(min(x_range, y_range) / range_divisor)
  mean_cells <- NA_real_
  iter <- 0L

  while (iter < max_iter) {
    mean_cells <- .mean_cells_per_window_grid(x, y, win_size)
    if (mean_cells >= min_cells && mean_cells <= max_cells) break
    if (mean_cells < min_cells) {
      win_size <- .round_to_even(win_size * 1.1)
    } else {
      win_size <- .round_to_even(win_size * 0.9)
    }
    iter <- iter + 1L
  }

  if (iter >= max_iter) {
    warning(
      "Optimal window size not obtained in ", max_iter,
      " iterations (mean cells per window: ", round(mean_cells, 2), ")."
    )
  }

  list(
    window_size = win_size,
    sliding_step = max(1L, win_size %/% 2L),
    mean_cells = mean_cells
  )
}

#' Resolve X/Y coordinate column names
#' @param df data.frame/data.table
#' @param x_col optional explicit X column name
#' @param y_col optional explicit Y column name
#' @return named character vector c(x = ..., y = ...)
#' @keywords internal
.resolve_xy_cols <- function(df, x_col = NULL, y_col = NULL) {
  nm <- names(df)
  pick <- function(explicit, candidates) {
    if (!is.null(explicit)) {
      if (!explicit %in% nm) {
        stop("Coordinate column '", explicit, "' not found in data")
      }
      return(explicit)
    }
    hit <- candidates[candidates %in% nm]
    if (length(hit)) return(hit[1L])
    stop(
      "Could not find coordinate columns. Specify x_col/y_col explicitly, ",
      "or include one of: ", paste(candidates, collapse = ", ")
    )
  }
  c(
    x = pick(x_col, c("X", "x")),
    y = pick(y_col, c("Y", "y"))
  )
}

#' Extract numeric X/Y coordinate matrix from a data frame
#' @param df data.frame/data.table
#' @param x_col optional explicit X column name
#' @param y_col optional explicit Y column name
#' @return numeric matrix with columns X and Y
#' @keywords internal
.extract_xy_matrix <- function(df, x_col = NULL, y_col = NULL) {
  cols <- .resolve_xy_cols(df, x_col, y_col)
  mat <- cbind(
    as.numeric(as.character(df[[cols["x"]]])),
    as.numeric(as.character(df[[cols["y"]]]))
  )
  colnames(mat) <- c("X", "Y")
  mat
}

#' Categorical Moran's I via proportion-weighted binary Moran's I per type
#' @param coords numeric matrix of X/Y coordinates
#' @param types character/factor vector of category labels
#' @param k number of nearest neighbors (excluding self)
#' @return scalar Moran's I aggregated across categories, or NA
#' @keywords internal
.cat_morans_I_knn <- function(coords, types, k) {
  n <- nrow(coords)
  if (n < 2L || length(unique(types)) < 2L) return(NA_real_)

  nnI <- RANN::nn2(coords, coords, k = k + 1L)
  nbr_idx <- nnI$nn.idx[, -1L, drop = FALSE]
  Wd <- 1 / (nnI$nn.dists[, -1L, drop = FALSE] + .Machine$double.eps)
  W <- Wd / rowSums(Wd)

  p_tab <- table(types) / n
  I_acc <- 0
  w_acc <- 0
  for (t in names(p_tab)) {
    z <- as.numeric(types == t)
    zc <- z - mean(z)
    denom <- sum(zc^2)
    if (denom < .Machine$double.eps) next
    Zlag <- rowSums(W * matrix(zc[nbr_idx], nrow = n))
    I_acc <- I_acc + as.numeric(p_tab[t]) * sum(zc * Zlag) / denom
    w_acc <- w_acc + as.numeric(p_tab[t])
  }
  if (w_acc == 0) return(NA_real_)
  I_acc / w_acc
}

#' Calculate spatial metrics for method selection
#' @param df data.table with spatial data
#' @param celltype_col cell type column name
#' @param x_col optional X coordinate column name (auto-detects X/x)
#' @param y_col optional Y coordinate column name (auto-detects Y/y)
#' @return list of spatial metrics
#' @keywords internal
.calc_spatial_metrics <- function(df, celltype_col = NULL, x_col = NULL, y_col = NULL) {
  if (!data.table::is.data.table(df)) {
    df <- data.table::as.data.table(df)
  }
  cols <- .resolve_xy_cols(df, x_col, y_col)
  x_vals <- as.numeric(as.character(df[[cols["x"]]]))
  y_vals <- as.numeric(as.character(df[[cols["y"]]]))
  coords <- cbind(x_vals, y_vals)
  n <- nrow(coords)
  k_nn <- .default_nn_rank()

  d_knn <- .nn_dist_at_rank(coords, k_nn = k_nn)
  mean_nn <- mean(d_knn)
  cv_nn   <- stats::sd(d_knn) / mean_nn

  # Clark-Evans R uses 1-NN distance (theoretical CSR expectation is for 1-NN)
  d1 <- RANN::nn2(coords, coords, k = 2)$nn.dists[, 2]

  # Clark-Evans R index under CSR: E[NN] = 0.5 / sqrt(lambda)
  rngX <- range(x_vals, na.rm = TRUE); rngY <- range(y_vals, na.rm = TRUE)
  area <- (diff(rngX) + 1e-8) * (diff(rngY) + 1e-8)
  lambda <- n / area
  E_nn <- 0.5 / sqrt(lambda)
  clark_evans_R <- mean(d1) / E_nn

  # Heterogeneity (Gini-Simpson)
  het_gs <- NA_real_
  if (!is.null(celltype_col) && celltype_col %in% names(df)) {
    p <- table(df[[celltype_col]]) / n
    het_gs <- 1 - sum((p)^2)
  }

  # Categorical Moran's I: binary Moran's I per cell type, weighted by type proportion
  cat_morans_I <- NA_real_
  if (!is.null(celltype_col) && celltype_col %in% names(df)) {
    kI <- max(5, min(20, round(sqrt(n))))
    cat_morans_I <- .cat_morans_I_knn(coords, df[[celltype_col]], k = kI)
  }

  list(
    mean_nn_dist = mean_nn,
    mean_nn_rank = k_nn,
    cv_nn_dist = cv_nn,
    clark_evans_R = clark_evans_R,
    cat_morans_I = cat_morans_I,
    het_gini_simpson = het_gs
  )
}

#' Select optimal method based on spatial metrics
#'
#' Decision tree (when called from `method = "auto"`):
#' \enumerate{
#'   \item n > 50000 -> window
#'   \item cat_morans_I >= 0.4 -> radius (strong type clustering)
#'   \item cv_nn_dist < 0.3 -> knn (CSR-like)
#'     or radius (mild deviation from CSR within the same band)
#'   \item otherwise -> delaunay (fallback)
#' }
#'
#' @param m spatial metrics
#' @param n number of cells
#' @return list with `method` (character) and `params` (list of recommended parameters)
#' @keywords internal
.select_method <- function(m, n) {
  params <- list()

  # 1) Large datasets -> windowed kNN
  if (n > 50000L) {
    params$k <- 6L
    return(list(method = "window", params = params))
  }

  morans <- m$cat_morans_I

  # 2) Strong spatial clustering of cell types -> radius
  if (!is.na(morans) && morans >= 0.4) {
    params$radius <- .default_search_radius(m$mean_nn_dist)
    return(list(method = "radius", params = params))
  }

  if (is.na(morans)) {
    warning("Categorical Moran's I unavailable; treating as weak clustering (< 0.4).")
  }

  # 3) Uniform density + CSR-like layout -> knn or radius
  cv_ok <- !is.na(m$cv_nn_dist) && m$cv_nn_dist < 0.3
  ce_r <- m$clark_evans_R
  ce_ok <- !is.na(ce_r) && ce_r >= 0.85 && ce_r <= 1.15

  if (cv_ok && ce_ok) {
    if (ce_r >= 0.95 && ce_r <= 1.05) {
      return(list(method = "knn", params = params))
    }
    params$radius <- .default_search_radius(m$mean_nn_dist)
    return(list(method = "radius", params = params))
  }

  # 4) Fallback -> Delaunay (parameter-light, robust on irregular layouts)
  params$max_edge_length <- m$mean_nn_dist * 2
  list(method = "delaunay", params = params)
}

#' Build edges using kNN method
#' @param df data.table with spatial data
#' @param k number of neighbors
#' @return data.table with edges
#' @keywords internal
.edges_knn <- function(df, k) {
  coords <- as.matrix(df[, .(X, Y)])
  nn <- RANN::nn2(coords, coords, k = k + 1)
  from_idx <- rep.int(seq_len(nrow(df)), k)
  to_idx   <- as.vector(t(nn$nn.idx[, -1, drop = FALSE]))
  d        <- as.vector(t(nn$nn.dists[, -1, drop = FALSE]))

  ed <- data.table::data.table(
    from = df$Cell_ID[from_idx],
    to   = df$Cell_ID[to_idx],
    dist = d
  )
  # Make undirected unique (keep shortest duplicate)
  ed[, key := ifelse(from < to, paste(from,to), paste(to,from))]
  ed <- ed[order(key, dist)][!duplicated(key)][, key := NULL][]
  ed
}

#' Build edges using radius method
#' @param df data.table with spatial data
#' @param radius search radius
#' @param k_cap maximum neighbors to consider
#' @return data.table with edges
#' @keywords internal
.edges_radius <- function(df, radius, k_cap = 50) {
  coords <- as.matrix(df[, .(X, Y)])
  # Expected neighbors ~ lambda * pi r^2; cap to k_cap for speed
  rngX <- range(df$X); rngY <- range(df$Y)
  area <- (diff(rngX) + 1e-8) * (diff(rngY) + 1e-8)
  lambda <- nrow(df) / area
  exp_deg <- max(5, min(k_cap, round(lambda * pi * radius^2)))

  nn <- RANN::nn2(coords, coords, k = exp_deg + 1)
  from_idx <- rep.int(seq_len(nrow(df)), exp_deg)
  to_idx   <- as.vector(t(nn$nn.idx[, -1, drop = FALSE]))
  d        <- as.vector(t(nn$nn.dists[, -1, drop = FALSE]))

  keep <- d <= radius
  ed <- data.table::data.table(
    from = df$Cell_ID[from_idx][keep],
    to   = df$Cell_ID[to_idx][keep],
    dist = d[keep]
  )
  ed[, key := ifelse(from < to, paste(from,to), paste(to,from))]
  ed <- ed[order(key, dist)][!duplicated(key)][, key := NULL][]
  ed
}

#' Build edges using Delaunay triangulation
#' @param df data.table with spatial data
#' @param max_edge_length optional maximum edge length filter
#' @return data.table with edges
#' @keywords internal
.edges_delaunay <- function(df, max_edge_length = NULL) {
  if (!requireNamespace("deldir", quietly = TRUE)) {
    warning("Package 'deldir' not installed; falling back to kNN (k=6)")
    return(.edges_knn(df, k = 6))
  }
  dd <- deldir::deldir(df$X, df$Y, suppressMsge = TRUE)
  te <- deldir::triang.list(dd)
  # Collect triangle edges
  ed_list <- lapply(te, function(tri) {
    idx <- tri$ptNum
    matrix(c(idx[1],idx[2], idx[2],idx[3], idx[3],idx[1]), ncol=2, byrow=TRUE)
  })
  ed_idx <- do.call(rbind, ed_list)
  ed_idx <- unique(t(apply(ed_idx, 1, function(x) sort(x))))
  d <- sqrt((df$X[ed_idx[,1]] - df$X[ed_idx[,2]])^2 + (df$Y[ed_idx[,1]] - df$Y[ed_idx[,2]])^2)
  ed <- data.table::data.table(
    from = df$Cell_ID[ed_idx[,1]],
    to   = df$Cell_ID[ed_idx[,2]],
    dist = d
  )
  if (!is.null(max_edge_length)) {
    ed <- ed[dist <= max_edge_length]
  }
  ed
}

#' Build edges using window method
#' @param df data.table with spatial data
#' @param k number of neighbors
#' @param tile window size
#' @return data.table with edges
#' @keywords internal
.edges_window <- function(df, k, tile) {
  DT <- df[, .(Cell_ID, X, Y)]
  # Tile assignment
  DT[, gx := floor(X / tile)]
  DT[, gy := floor(Y / tile)]
  # Neighborhood map (self + 8 neighbors)
  neigh <- data.table::CJ(dx = -1:1, dy = -1:1)

  res <- vector("list", 0L)
  # Iterate unique tiles; keep memory local
  tiles <- unique(DT[, .(gx, gy)])
  for (i in seq_len(nrow(tiles))) {
    tx <- tiles$gx[i]; ty <- tiles$gy[i]
    # Gather cells in tile neighborhood
    nb <- data.table::CJ(gx = tx + neigh$dx,
                         gy = ty + neigh$dy,
                         unique = TRUE)
    sub <- DT[nb, on = .(gx, gy), nomatch = 0L]
    if (nrow(sub) < 2) next
    # Run small-k kNN within local block
    coords <- as.matrix(sub[, .(X, Y)])
    kk <- min(k + 1, max(2, nrow(sub)))
    nn <- RANN::nn2(coords, coords, k = kk)
    from_idx <- rep.int(seq_len(nrow(sub)), kk - 1)
    to_idx   <- as.vector(t(nn$nn.idx[, -1, drop = FALSE]))
    d        <- as.vector(t(nn$nn.dists[, -1, drop = FALSE]))
    ed <- data.table::data.table(
      from = sub$Cell_ID[from_idx],
      to   = sub$Cell_ID[to_idx],
      dist = d
    )
    res[[length(res) + 1L]] <- ed
  }
  ed <- data.table::rbindlist(res, use.names = TRUE, fill = TRUE)
  if (nrow(ed) == 0L) return(ed)
  ed[, key := ifelse(from < to, paste(from,to), paste(to,from))]
  ed <- ed[order(key, dist)][!duplicated(key)][, key := NULL][]
  ed
}

#' Calculate neighborhood composition features
#'
#' @param df Spatial data with cell types
#' @param edges Data.frame of edges from build_spatial_network()
#' @param cell_id_col Column name for cell IDs (default: "Cell_ID")
#' @param celltype_col Column name for cell types (default: "celltype")
#' @return Enhanced data.table with neighborhood type proportions
#' @export
calculate_neighborhood_features <- function(df, edges,
                                            cell_id_col = "Cell_ID",
                                            celltype_col = "celltype") {
  if (!data.table::is.data.table(df)) {
    df <- data.table::as.data.table(df)
  }

  # Verify required columns exist
  required_cols <- c(cell_id_col, celltype_col, "X", "Y")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Ensure edge data is data.table
  if (!data.table::is.data.table(edges)) {
    edges <- data.table::as.data.table(edges)
  }

  # Verify edge data contains correct columns
  if (!all(c("from", "to") %in% colnames(edges))) {
    stop("edges must contain 'from' and 'to' columns")
  }

  # Create ID to row number mapping
  df[, row_id := .I]
  id_to_row <- setNames(df$row_id, df[[cell_id_col]])

  # Build adjacency mapping (using actual IDs not row numbers)
  edge_map <- list()
  unique_ids <- unique(c(edges$from, edges$to))

  for (id in unique_ids) {
    neighbors_from <- edges[from == id, to]
    neighbors_to <- edges[to == id, from]
    edge_map[[as.character(id)]] <- unique(c(neighbors_from, neighbors_to))
  }

  # Get all cell types
  all_types <- unique(df[[celltype_col]])
  type_cols <- paste0("Nbr_", all_types)

  # Initialize feature matrix (all zeros)
  nbr_features <- matrix(0, nrow = nrow(df), ncol = length(all_types) + 2)
  colnames(nbr_features) <- c(type_cols, "Centroid_X", "Centroid_Y")

  # Process each cell
  for (i in 1:nrow(df)) {
    cell_id <- df[[cell_id_col]][i]
    neighbors <- edge_map[[as.character(cell_id)]]

    if (is.null(neighbors) || length(neighbors) == 0) next

    # Get neighbor row numbers
    neighbor_rows <- id_to_row[as.character(neighbors)]
    neighbor_rows <- neighbor_rows[!is.na(neighbor_rows)]

    if (length(neighbor_rows) == 0) next

    nbr_data <- df[neighbor_rows]

    # Calculate neighbor type proportions
    type_counts <- table(factor(nbr_data[[celltype_col]], levels = all_types))
    nbr_features[i, type_cols] <- type_counts / sum(type_counts)

    # Calculate centroid
    centroid_x <- mean(nbr_data$X, na.rm = TRUE)
    centroid_y <- mean(nbr_data$Y, na.rm = TRUE)
    nbr_features[i, c("Centroid_X", "Centroid_Y")] <- c(centroid_x, centroid_y)
  }

  # Add new features to data
  new_features <- data.table::as.data.table(nbr_features)
  df <- cbind(df, new_features)

  # Add relative position features
  df[, `:=`(
    Rel_X = Centroid_X - X,
    Rel_Y = Centroid_Y - Y
  )]

  # Handle missing values (replace with 0)
  nbr_cols <- c(type_cols, "Centroid_X", "Centroid_Y", "Rel_X", "Rel_Y")
  for (col in nbr_cols) {
    df[is.na(get(col)), (col) := 0]
  }

  # Remove temporary row number column
  df[, row_id := NULL]

  return(df)
}

#' Cluster neighborhoods using combined spatial and compositional features
#'
#' @param feature_df Data.table with neighborhood features
#' @param spatial_edges Spatial edges from build_spatial_network()
#' @param method Clustering method ("kmeans", "hdbscan", or "louvain")
#' @param k Number of clusters (for kmeans)
#' @param use_pca Whether to use PCA for dimensionality reduction
#' @param var_threshold Variance threshold for PCA components
#' @param n_components Explicit number of PCA components (overrides var_threshold)
#' @param min_cluster_size Minimum points per cluster (for hdbscan)
#' @param cluster_colname Name for the output cluster column
#' @return Data.table with cluster assignments in cluster_colname
#' @export
cluster_neighborhoods <- function(feature_df,
                                  spatial_edges = NULL,
                                  method = "kmeans",
                                  k = 10,
                                  use_pca = TRUE,
                                  var_threshold = 0.9,
                                  n_components = NULL,
                                  min_cluster_size = 5,
                                  cluster_colname = "Neighborhood_Cluster") {

  # Basic checks
  if (!data.table::is.data.table(feature_df)) {
    feature_df <- data.table::as.data.table(feature_df)
  }

  if (nrow(feature_df) < 2) stop("Feature dataframe has less than 2 rows")

  # Auto-detect numeric feature columns
  feat_cols <- names(feature_df)[sapply(feature_df, is.numeric)]
  feat_cols <- setdiff(feat_cols, c("Cell_ID", "X", "Y", cluster_colname))
  if (length(feat_cols) == 0) stop("No valid numeric features found.")

  # Build feature matrix
  features <- as.matrix(feature_df[, ..feat_cols])

  # 1. Remove zero-variance features
  col_var <- apply(features, 2, stats::var, na.rm = TRUE)
  zero_var_cols <- which(col_var <= .Machine$double.eps)
  if (length(zero_var_cols) > 0) {
    message("Removing ", length(zero_var_cols), " zero-variance features.")
    features <- features[, -zero_var_cols, drop = FALSE]
  }
  if (ncol(features) == 0) stop("All features have zero variance.")

  # 2. Impute missing values
  na_count <- sum(is.na(features))
  if (na_count > 0) {
    message("Imputing ", na_count, " missing values with column means.")
    features <- apply(features, 2, function(x) {
      x[is.na(x)] <- mean(x, na.rm = TRUE)
      x
    })
  }

  # 3. Scale
  features <- scale(features)

  # 4. Remove columns that became NaN (zero std)
  nan_cols <- which(apply(features, 2, function(x) any(is.nan(x))))
  if (length(nan_cols) > 0) {
    message("Removing ", length(nan_cols), " NaN columns after scaling.")
    features <- features[, -nan_cols, drop = FALSE]
  }
  if (ncol(features) == 0) stop("No valid features after preprocessing.")

  # PCA
  if (use_pca) {
    if (is.null(n_components)) {
      pca <- stats::prcomp(features, scale. = FALSE)   # Already scaled
      var_exp <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
      n_components <- which(var_exp >= var_threshold)[1]
      message("Using ", n_components, " PCA components (",
              round(var_threshold * 100, 1), "% variance)")
    }
    if (nrow(features) > 1000 && requireNamespace("irlba", quietly = TRUE)) {
      pca <- irlba::prcomp_irlba(features, n = n_components, scale. = FALSE)
    } else {
      pca <- stats::prcomp(features, scale. = FALSE, rank. = n_components)
    }
    features <- pca$x
  }

  # Clustering
  clusters <- switch(method,
                     "kmeans" = {
                       set.seed(123)
                       km <- stats::kmeans(features, centers = k, nstart = 25, iter.max = 100)
                       km$cluster
                     },
                     "hdbscan" = {
                       if (!requireNamespace("dbscan", quietly = TRUE))
                         stop("dbscan package required for HDBSCAN")
                       minPts <- max(min_cluster_size,
                                     min(50, nrow(features) %/% 100))
                       hdb <- dbscan::hdbscan(features, minPts = minPts)
                       hdb$cluster
                     },
                     "louvain" = {
                       if (is.null(spatial_edges))
                         stop("spatial_edges required for Louvain clustering")
                       if (!requireNamespace("igraph", quietly = TRUE))
                         stop("igraph package required for Louvain clustering")

                       # Build graph
                       graph <- igraph::graph_from_data_frame(spatial_edges,
                                                              directed = FALSE)
                       v_in_graph <- igraph::V(graph)$name   # Character vertex IDs

                       # Align feature_df to graph vertices
                       feature_df_aligned <- feature_df[match(v_in_graph,
                                                              feature_df$Cell_ID), ]
                       if (anyNA(feature_df_aligned$Cell_ID))
                         stop("Some graph vertices are missing from feature_df")

                       # Add vertex attributes
                       for (col in setdiff(names(feature_df), "Cell_ID")) {
                         igraph::vertex_attr(graph, name = col) <-
                           feature_df_aligned[[col]]
                       }

                       # Louvain clustering
                       lv <- igraph::cluster_louvain(graph)
                       clusters <- lv$membership

                       # Add cluster labels to the *aligned* subset
                       feature_df_aligned[, (cluster_colname) := paste0("CN", clusters)]

                       # Return the aligned subset only
                       return(feature_df_aligned)
                     },
                     stop("Unsupported clustering method: ", method))

  # Non-louvain methods: add labels to original feature_df
  cluster_labels <- paste0("CN", clusters)
  if (cluster_colname %in% names(feature_df))
    warning("Overwriting existing column: ", cluster_colname)
  feature_df[, (cluster_colname) := cluster_labels]

  # Attach metadata
  attr(feature_df, "cluster_info") <- list(
    method = method,
    n_clusters = length(unique(clusters)),
    features_used = feat_cols,
    pca_used = use_pca,
    pca_components = if (use_pca) ncol(features) else NULL
  )

  return(feature_df)
}

#' Calculate neighborhood purity using flexible neighbor definitions
#'
#' @param df A data.frame or data.table containing spatial coordinates and cell type labels.
#' @param x_col Character, name of the X-coordinate column.
#' @param y_col Character, name of the Y-coordinate column.
#' @param celltype_col Character, name of the cell type column.
#' @param method Character, neighbor definition: "window", "radius", "knn", or "delaunay".
#' @param k Integer, number of nearest neighbors for method = "knn" (ignored otherwise).
#' @param radius Numeric, distance threshold for methods "radius" and "window" (ignored otherwise).
#' @param min_cells Integer, minimum number of neighbors required to compute purity.
#' @param verbose Logical, print progress messages.
#' @return A data.table with an added column `Neighborhood_Purity`.
#' @export
calculate_neighborhood_purity <- function(df,
                                          x_col = "X",
                                          y_col = "Y",
                                          celltype_col = "celltype",
                                          method = c("window", "radius", "knn", "delaunay"),
                                          k = 30,
                                          radius = 30,
                                          min_cells = 5,
                                          verbose = TRUE) {

  ## ---- 0. Ensure data.table with zero-copy --------------------------
  if (!data.table::is.data.table(df)) {
    df <- data.table::as.data.table(df)
  }

  ## ---- 1. Check required columns ------------------------------------
  required_cols <- c(x_col, y_col, celltype_col)
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols))
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))

  method <- match.arg(method)

  if (verbose) {
    message("Method: ", method)
    message("Total cells: ", nrow(df))
  }

  ## ---- 2. Extract vectors for speed ---------------------------------
  x <- df[[x_col]]
  y <- df[[y_col]]
  celltype <- df[[celltype_col]]
  purity <- rep(NA_real_, length(x))

  ## ---- 3. Branch by method ------------------------------------------
  if (method == "window") {
    half <- radius / 2
    for (i in seq_along(x)) {
      if (verbose && i %% 1000 == 0) message("Processed ", i, "/", length(x))
      mask <- (x >= x[i] - half) & (x <= x[i] + half) &
        (y >= y[i] - half) & (y <= y[i] + half)
      nbs <- which(mask)
      nbs <- setdiff(nbs, i)              # Exclude self
      if (length(nbs) < min_cells) next
      tab <- table(celltype[nbs])
      purity[i] <- max(tab) / sum(tab)
    }
  }

  if (method == "radius") {
    for (i in seq_along(x)) {
      if (verbose && i %% 1000 == 0) message("Processed ", i, "/", length(x))
      dx <- x - x[i]
      dy <- y - y[i]
      dist <- sqrt(dx^2 + dy^2)
      nbs <- which(dist <= radius & dist > 0)
      if (length(nbs) < min_cells) next
      tab <- table(celltype[nbs])
      purity[i] <- max(tab) / sum(tab)
    }
  }

  if (method == "knn") {
    if (!requireNamespace("RANN", quietly = TRUE))
      stop("Package 'RANN' is required for knn.")
    coords <- cbind(x, y)
    nn_idx <- RANN::nn2(coords, coords, k + 1)$nn.idx[, -1]  # Drop self
    for (i in seq_len(nrow(nn_idx))) {
      if (verbose && i %% 1000 == 0) message("Processed ", i, "/", nrow(nn_idx))
      nbs <- nn_idx[i, ]
      nbs <- nbs[nbs > 0]                    # Safety
      if (length(nbs) < min_cells) next
      tab <- table(celltype[nbs])
      purity[i] <- max(tab) / sum(tab)
    }
  }

  if (method == "delaunay") {
    if (!requireNamespace("deldir", quietly = TRUE))
      stop("Package 'deldir' is required for delaunay.")
    del <- deldir::deldir(x, y, suppressMsge = TRUE)
    edges <- rbind(
      data.frame(from = del$delsgs$ind1, to = del$delsgs$ind2),
      data.frame(from = del$delsgs$ind2, to = del$delsgs$ind1)
    )
    adj <- split(edges$to, edges$from)
    for (i in seq_along(x)) {
      if (verbose && i %% 1000 == 0) message("Processed ", i, "/", length(x))
      nbs <- unlist(adj[as.character(i)], use.names = FALSE)
      if (length(nbs) < min_cells) next
      tab <- table(celltype[nbs])
      purity[i] <- max(tab) / sum(tab)
    }
  }

  ## ---- 4. Attach result and return ----------------------------------
  df[, Neighborhood_Purity := purity]

  if (verbose) {
    message("Purity calculation completed")
    message("Summary of Neighborhood_Purity:")
    print(summary(df$Neighborhood_Purity))
    message("Cells with NA purity (insufficient neighbors): ", sum(is.na(df$Neighborhood_Purity)))
  }

  return(df)
}

#' Analyze spatial interactions between cell types
#'
#' @param df Spatial data with cell types
#' @param edges Spatial network edges (data.frame with "from" and "to")
#' @param celltype_col Column name for cell types in metadata.
#' @return List with interaction matrix and network graph
#' @export
analyze_spatial_interactions <- function(df, edges,
                                         celltype_col = "celltype") {
  # Basic checks
  if (!celltype_col %in% names(df)) stop("'", celltype_col, "' column not found in df")
  if (!"Cell_ID" %in% names(df)) stop("'Cell_ID' column not found in df")
  if (is.null(edges) || !is.data.frame(edges) || nrow(edges) == 0)
    stop("'edges' is empty or not a valid data.frame")
  if (!all(c("from", "to") %in% names(edges)))
    stop("'edges' must contain 'from' and 'to' columns")

  # Ensure IDs are character type
  df$Cell_ID <- as.character(df$Cell_ID)
  edges$from <- as.character(edges$from)
  edges$to <- as.character(edges$to)

  # Build mapping table
  cell_type_map <- setNames(df[[celltype_col]], df$Cell_ID)

  # Get cell types of interacting pairs
  edges$type1 <- cell_type_map[edges$from]
  edges$type2 <- cell_type_map[edges$to]

  # Remove NA or unknown types
  valid_edges <- edges[!is.na(edges$type1) & !is.na(edges$type2), ]

  # Build interaction frequency table
  inter_table <- table(pmin(valid_edges$type1, valid_edges$type2),
                       pmax(valid_edges$type1, valid_edges$type2))

  # Convert to symmetric matrix
  all_types <- sort(unique(df[[celltype_col]]))
  int_mat <- matrix(0, nrow = length(all_types), ncol = length(all_types),
                    dimnames = list(all_types, all_types))
  for (i in seq_len(nrow(inter_table))) {
    type1 <- rownames(inter_table)[i]
    for (j in seq_len(ncol(inter_table))) {
      type2 <- colnames(inter_table)[j]
      count <- inter_table[i, j]
      int_mat[type1, type2] <- int_mat[type1, type2] + count
      int_mat[type2, type1] <- int_mat[type2, type1] + count
    }
  }

  # Build graph
  g <- igraph::graph_from_adjacency_matrix(int_mat, mode = "undirected", weighted = TRUE)

  return(list(interaction_matrix = int_mat, network = g))
}


#' Generate color palette for visualizations
#'
#' @param n Number of colors needed
#' @return Vector of color codes
#' @export
#'
#' @description
#' Creates optimized color palettes:
#' 1. Small sets: ColorBrewer Set1
#' 2. Medium sets: Set3
#' 3. Large sets: Hue palette
#' Ensures distinct colors for all categories
get_color_palette <- function(n) {
  if (n <= 2) {
    return(c("#1f77b4", "#ff7f0e")[1:n])  # Minimal palette
  } else if (n <= 9) {
    return(RColorBrewer::brewer.pal(n, "Set1"))  # ColorBrewer qualitative
  } else if (n <= 12) {
    return(RColorBrewer::brewer.pal(n, "Set3"))  # Extended qualitative
  } else {
    return(scales::hue_pal()(n))  # Large palette
  }
}

#' Visualize spatial cell type distribution
#'
#' @param df Spatial data
#' @param x_col X coordinate column name (default: "X")
#' @param y_col Y coordinate column name (default: "Y")
#' @param celltype_col Cell type column name (default: "celltype")
#' @param point_size Point size (default: 1.5)
#' @param point_alpha Point transparency (default: 0.6)
#' @param point_shape Point shape (default: 16)
#' @param legend_point_size Legend point size (default: 3)
#' @param title Plot title
#' @param legend.position Legend position (default: "right")
#' @param save_path Output file path (optional)
#' @param width Plot width in inches (default: 10)
#' @param height Plot height in inches (default: 8)
#' @return ggplot object and saves plot to file if save_path provided
#' @export
visualize_spatial_distribution <- function(df,
                                           x_col = "X",
                                           y_col = "Y",
                                           celltype_col = "celltype",
                                           point_size = 1.5,
                                           point_alpha = 0.6,
                                           point_shape = 16,
                                           legend_point_size = 3,
                                           title = NULL,
                                           legend.position = "right",
                                           save_path = NULL,
                                           width = 10,
                                           height = 8) {

  set.seed(123)

  # Create local copies of columns
  df$plot_X <- df[[x_col]]
  df$plot_Y <- df[[y_col]]
  df$plot_celltype <- df[[celltype_col]]

  # Generate color palette
  n_types <- length(unique(df$plot_celltype))
  palette <- get_color_palette(n_types)

  # Create spatial plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = plot_X, y = plot_Y, color = plot_celltype)) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha, shape = point_shape) +
    ggplot2::scale_color_manual(
      values = palette,
      guide = ggplot2::guide_legend(
        override.aes = list(
          size = legend_point_size,
          alpha = 1,
          shape = point_shape
        )
      )
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = if (!is.null(title)) title else "Spatial Distribution of Cell Types",
      x = "X Coordinate",
      y = "Y Coordinate",
      color = celltype_col
    ) +
    ggplot2::coord_fixed() +
    ggplot2::theme(legend.position = legend.position)

  # Save plot if path provided
  if (!is.null(save_path)) {
    ggplot2::ggsave(save_path, p, width = width, height = height)
    message("Plot saved to: ", save_path)
  }

  return(p)
}

#' Visualize distance heatmap between cell types
#'
#' @param dist_result Distance matrix result from calculate_celltype_distances()
#' @param save_path Output file path (optional)
#' @param width Plot width in inches (default: 12)
#' @param height Plot height in inches (default: 10)
#' @return ggplot object and saves plot to file if save_path provided
#' @export
visualize_distance_heatmap <- function(dist_result,
                                       save_path = NULL,
                                       width = 12,
                                       height = 10) {

  # Prepare data
  dist_matrix <- dist_result$distance_matrix
  dist_df <- reshape2::melt(dist_matrix, varnames = c("Source", "Target"), value.name = "Distance")

  # Calculate optimal text color
  dist_df$TextColor <- ifelse(
    dist_df$Distance < stats::quantile(dist_df$Distance, 0.5, na.rm = TRUE),
    "black",
    "white"
  )

  # Create heatmap
  p <- ggplot2::ggplot(dist_df, ggplot2::aes(x = Source, y = Target, fill = Distance)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(
      ggplot2::aes(label = round(Distance, 1), color = TextColor),
      size = 3.5,
      fontface = "bold",
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_gradientn(
      name = "Mean Distance",
      colours = c("#FFFFE0", "#FFD700", "#FF8C00", "#FF4500", "#8B0000", "#4B0082", "#2E0854"),
      values = scales::rescale(c(
        min(dist_df$Distance, na.rm = TRUE),
        stats::quantile(dist_df$Distance, 0.2, na.rm = TRUE),
        stats::quantile(dist_df$Distance, 0.4, na.rm = TRUE),
        stats::quantile(dist_df$Distance, 0.5, na.rm = TRUE),
        stats::quantile(dist_df$Distance, 0.7, na.rm = TRUE),
        stats::quantile(dist_df$Distance, 0.9, na.rm = TRUE),
        max(dist_df$Distance, na.rm = TRUE)
      ))
    ) +
    ggplot2::scale_color_identity() +
    ggplot2::labs(
      title = "Mean Distance Between Cell Types",
      x = "Source Cell Type",
      y = "Target Cell Type"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 10, face = "bold"),
      axis.title = ggplot2::element_text(size = 12, face = "bold"),
      plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5),
      legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold"),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::geom_tile(
      data = dist_df[dist_df$Source == dist_df$Target, ],
      ggplot2::aes(fill = Distance),
      color = "black",
      linewidth = 1
    )

  # Save plot if path provided
  if (!is.null(save_path)) {
    ggplot2::ggsave(save_path, p, width = width, height = height)
    message("Heatmap saved to: ", save_path)
  }

  return(p)
}

#' Visualize distance relationships using parallel coordinates plot
#'
#' @param dist_result Distance matrix result from calculate_celltype_distances()
#' @param save_path Output file path (optional)
#' @param width Plot width in inches (default: 14)
#' @param height Plot height in inches (default: 8)
#' @return ggplot object and saves plot to file if save_path provided
#' @export
visualize_distance_parallel <- function(dist_result,
                                        save_path = NULL,
                                        width = 14,
                                        height = 8) {

  # Prepare data
  dist_matrix <- dist_result$distance_matrix
  cell_types <- rownames(dist_matrix)

  # Convert to long format
  plot_data <- data.frame(dist_matrix)
  plot_data$Source <- cell_types
  plot_data <- reshape2::melt(plot_data, id.vars = "Source",
                              variable.name = "Target",
                              value.name = "Distance")

  # Create color scheme
  n_types <- length(cell_types)
  colors <- viridis::viridis(n_types, option = "D", direction = -1)

  # Create parallel coordinates plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Target, y = Distance, group = Source, color = Source)) +
    ggplot2::geom_line(size = 1.2, alpha = 0.9) +
    ggplot2::geom_point(size = 3, alpha = 0.9) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(
      title = "Cell Type Distance Relationships",
      x = "Target Cell Type",
      y = "Distance to Target"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5,
                                         margin = ggplot2::margin(b = 15)),
      axis.title = ggplot2::element_text(size = 12, face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 9),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold", size = 10),
      legend.text = ggplot2::element_text(size = 9),
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(title = "Source Cell Type",
                                                  title.position = "top"))

  # Save plot if path provided
  if (!is.null(save_path)) {
    ggplot2::ggsave(save_path, p, width = width, height = height)
    message("Parallel coordinates plot saved to: ", save_path)
  }

  return(p)
}

#' Visualize neighborhood purity
#'
#' @param df Spatial data with Neighborhood_Purity column
#' @param x_col X coordinate column name (default: "X")
#' @param y_col Y coordinate column name (default: "Y")
#' @param max_points Maximum points to plot (default: 10000)
#' @param point_size Point size (default: 1.5)
#' @param point_alpha Point transparency (default: 0.8)
#' @param point_shape Point shape (default: 16)
#' @param title Plot title
#' @param legend.position Legend position (default: "right")
#' @param save_path Output file path (optional)
#' @param width Plot width in inches (default: 10)
#' @param height Plot height in inches (default: 8)
#' @return ggplot object and saves plot to file if save_path provided
#' @export
visualize_neighborhood_purity <- function(df,
                                          x_col = "X",
                                          y_col = "Y",
                                          max_points = 10000,
                                          point_size = 1.5,
                                          point_alpha = 0.8,
                                          point_shape = 16,
                                          title = NULL,
                                          legend.position = "right",
                                          save_path = NULL,
                                          width = 10,
                                          height = 8) {

  # Sample large datasets
  if (nrow(df) > max_points) {
    set.seed(123)
    df <- df[sample(nrow(df), max_points), ]
  }

  # Validate required column
  if (!"Neighborhood_Purity" %in% names(df)) {
    stop("Neighborhood_Purity column not found")
  }

  # Create local copies of columns
  df$plot_X <- df[[x_col]]
  df$plot_Y <- df[[y_col]]

  # Create purity visualization
  p <- ggplot2::ggplot(df, ggplot2::aes(x = plot_X, y = plot_Y, color = Neighborhood_Purity)) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha, shape = point_shape) +
    ggplot2::scale_color_viridis_c(option = "magma", direction = -1) +
    ggplot2::theme_minimal() +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title = if (!is.null(title)) title else "Neighborhood Purity",
      x = "X Coordinate",
      y = "Y Coordinate"
    ) +
    ggplot2::theme(legend.position = legend.position)

  # Save plot if path provided
  if (!is.null(save_path)) {
    ggplot2::ggsave(save_path, p, width = width, height = height)
    message("Neighborhood purity plot saved to: ", save_path)
  }

  return(p)
}

#' Calculate cluster composition metrics
#'
#' @param df Spatial data
#' @param cluster_col Cluster column name (default: "Neighborhood_Cluster")
#' @param celltype_col Cell type column name (default: "celltype")
#' @return Data frame with cluster composition statistics
#' @export
calculate_cluster_composition <- function(df,
                                          cluster_col = "Neighborhood_Cluster",
                                          celltype_col = "celltype") {
  # Calculate counts and proportions
  composition <- df %>%
    dplyr::group_by(!!rlang::sym(cluster_col), !!rlang::sym(celltype_col)) %>%
    dplyr::summarise(count = dplyr::n(), .groups = 'drop') %>%
    dplyr::group_by(!!rlang::sym(cluster_col)) %>%
    dplyr::mutate(total = sum(.data$count),
                  proportion = .data$count / .data$total) %>%
    dplyr::ungroup()

  return(composition)
}

#' Visualize cluster composition as heatmap
#'
#' @param composition_df Composition data from calculate_cluster_composition()
#' @param cluster_col Cluster column name (default: "Neighborhood_Cluster")
#' @param celltype_col Cell type column name (default: "celltype")
#' @param value_col Value column to visualize (default: "proportion")
#' @param save_path Output file path (optional)
#' @param width Plot width in inches (default: 12)
#' @param height Plot height in inches (default: 10)
#' @return ComplexHeatmap object and saves plot to file if save_path provided
#' @export
plot_composition_heatmap <- function(composition_df,
                                     cluster_col = "Neighborhood_Cluster",
                                     celltype_col = "celltype",
                                     value_col = "proportion",
                                     save_path = NULL,
                                     width = 12,
                                     height = 10) {
  # Check required columns
  required_cols <- c(cluster_col, celltype_col, value_col, "count")
  missing_cols <- setdiff(required_cols, colnames(composition_df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Convert cluster column to character
  composition_df[[cluster_col]] <- as.character(composition_df[[cluster_col]])

  # Remove NA values
  composition_df <- composition_df %>%
    dplyr::filter(!is.na(!!rlang::sym(cluster_col)), !is.na(!!rlang::sym(celltype_col)))

  # Create composition matrix
  comp_matrix <- composition_df %>%
    tidyr::pivot_wider(
      id_cols = !!rlang::sym(cluster_col),
      names_from = !!rlang::sym(celltype_col),
      values_from = !!rlang::sym(value_col),
      values_fill = 0
    ) %>%
    as.data.frame()

  # Set row names
  rownames(comp_matrix) <- comp_matrix[[cluster_col]]
  comp_matrix <- comp_matrix[, !(names(comp_matrix) %in% cluster_col), drop = FALSE]
  comp_matrix <- as.matrix(comp_matrix)

  # Ensure values are within [0,1]
  comp_matrix[comp_matrix < 0] <- 0
  comp_matrix[comp_matrix > 1] <- 1

  # Apply Z-score normalization across rows (clusters)
  z_matrix <- t(apply(comp_matrix, 1, scale))
  colnames(z_matrix) <- colnames(comp_matrix)
  rownames(z_matrix) <- rownames(comp_matrix)

  # Define color function for Z-scores
  z_range <- stats::quantile(z_matrix, c(0.05, 0.95), na.rm = TRUE)
  max_abs <- max(abs(z_range))
  col_limits <- c(-max_abs, max_abs)
  col_fun <- circlize::colorRamp2(
    c(col_limits[1], 0, col_limits[2]),
    c("#4575B4", "white", "#D73027")
  )

  # Calculate cluster sizes
  cluster_counts <- composition_df %>%
    dplyr::group_by(!!rlang::sym(cluster_col)) %>%
    dplyr::summarise(total = sum(.data$count)) %>%
    tibble::deframe() %>%
    .[rownames(z_matrix)]

  # Create row annotations
  row_ha <- ComplexHeatmap::rowAnnotation(
    "Cell Count" = ComplexHeatmap::anno_barplot(
      cluster_counts,
      bar_width = 0.8,
      gp = grid::gpar(fill = "grey70", col = NA),
      width = grid::unit(2, "cm")
    )
  )

  # Generate heatmap
  hm <- ComplexHeatmap::Heatmap(
    z_matrix,
    name = "Z-score",
    col = col_fun,
    row_title = "Neighborhood Cluster",
    column_title = "Cell Type Composition",
    show_row_names = TRUE,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "ward.D2",
    row_names_gp = grid::gpar(fontsize = 10),
    column_names_gp = grid::gpar(fontsize = 10),
    heatmap_legend_param = list(title_position = "topcenter"),
    right_annotation = row_ha,
    row_dend_width = grid::unit(1.5, "cm"),
    column_dend_height = grid::unit(1.5, "cm"),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid::grid.text(sprintf("%.2f", comp_matrix[i, j]),
                      x, y, gp = grid::gpar(fontsize = 8, col = "black"))
    }
  )

  # Save heatmap if path provided
  if (!is.null(save_path)) {
    grDevices::pdf(save_path, width = width, height = height)
    ComplexHeatmap::draw(hm)
    grDevices::dev.off()
    message("Composition heatmap saved to: ", save_path)
  }

  return(hm)
}

#' Visualize cluster composition as bar plot
#'
#' @param composition_df Composition data from calculate_cluster_composition()
#' @param cluster_col Cluster column name (default: "Neighborhood_Cluster")
#' @param celltype_col Cell type column name (default: "celltype")
#' @param value_col Value column to visualize (default: "proportion")
#' @param save_path Output file path (optional)
#' @param width Plot width in inches (default: 12)
#' @param height Plot height in inches (default: 8)
#' @return ggplot object and saves plot to file if save_path provided
#' @export
plot_composition_barplot <- function(composition_df,
                                     cluster_col = "Neighborhood_Cluster",
                                     celltype_col = "celltype",
                                     value_col = "proportion",
                                     save_path = NULL,
                                     width = 12,
                                     height = 8) {

  # Generate color palette
  celltypes <- unique(composition_df[[celltype_col]])
  palette <- get_color_palette(length(celltypes))
  names(palette) <- celltypes

  # Order cell types by frequency
  celltype_freq <- composition_df %>%
    dplyr::group_by(!!rlang::sym(celltype_col)) %>%
    dplyr::summarise(total = sum(.data$count)) %>%
    dplyr::arrange(dplyr::desc(.data$total)) %>%
    dplyr::pull(!!rlang::sym(celltype_col))

  composition_df[[celltype_col]] <- factor(
    composition_df[[celltype_col]],
    levels = rev(celltype_freq)
  )

  # Create bar plot
  p <- ggplot2::ggplot(composition_df,
                       ggplot2::aes(x = factor(!!rlang::sym(cluster_col)),
                                    y = !!rlang::sym(value_col),
                                    fill = !!rlang::sym(celltype_col))) +
    ggplot2::geom_bar(stat = "identity", position = "fill", width = 0.8) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(),
                                expand = c(0, 0)) +
    ggplot2::labs(
      x = "Neighborhood Cluster",
      y = "Cell Type Proportion",
      fill = "Cell Type"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = ggplot2::element_text(size = 10),
      axis.title = ggplot2::element_text(size = 12, face = "bold"),
      legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold"),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(10, 10, 10, 20)
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE))

  # Save plot if path provided
  if (!is.null(save_path)) {
    ggplot2::ggsave(save_path, p, width = width, height = height)
    message("Composition bar plot saved to: ", save_path)
  }

  return(p)
}

#' Visualize cell-cell interaction matrix
#'
#' @param interaction_matrix Interaction matrix from analyze_spatial_interactions()
#' @param transform Apply log2 transformation (default: TRUE)
#' @param color_palette Color palette function (default: viridis::inferno)
#' @param save_path Output file path (optional)
#' @param width Plot width in inches (default: 10)
#' @param height Plot height in inches (default: 8)
#' @return pheatmap object and saves plot to file if save_path provided
#' @export
visualize_interaction_heatmap <- function(interaction_matrix,
                                          transform = TRUE,
                                          color_palette = viridis::inferno,
                                          save_path = NULL,
                                          width = 10,
                                          height = 8) {

  mat <- interaction_matrix

  # Optional transformation
  if (isTRUE(transform)) {
    mat <- log2(mat + 1)
  }

  # Create color palette
  col_fun <- color_palette(256)

  # Create heatmap
  hm <- pheatmap::pheatmap(
    mat = mat,
    main = "Cell-Cell Interaction Frequency",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    fontsize_row = 9,
    fontsize_col = 9,
    display_numbers = TRUE,
    number_format = "%.0f",
    color = col_fun,
    border_color = NA,
    silent = is.null(save_path)
  )

  # Save heatmap if path provided
  if (!is.null(save_path)) {
    grDevices::pdf(save_path, width = width, height = height)
    grid::grid.newpage()
    grid::grid.draw(hm$gtable)
    grDevices::dev.off()
    message("Interaction heatmap saved to: ", save_path)
  }

  return(hm)
}

#' Visualize spatial network with flexible edge selection
#'
#' @param df Spatial data
#' @param edges Spatial network edges
#' @param celltype_col Cell type column name (default: "celltype")
#' @param x_col X coordinate column name (default: "X")
#' @param y_col Y coordinate column name (default: "Y")
#' @param edge_mode Edge selection mode: "all", "top", or "random" (default: "all")
#' @param top_n When edge_mode = "top", keep this many strongest edges (default: 1000)
#' @param max_edges When edge_mode = "random", sample this many edges (default: 1000)
#' @param point_size Point size (default: 1.5)
#' @param point_alpha Point transparency (default: 0.7)
#' @param edge_size_range Edge size range (default: c(0.3, 2.0))
#' @param edge_alpha_range Edge alpha range (default: c(0.3, 0.9))
#' @param edge_color Edge color (default: "darkred")
#' @param show_points Whether to show points (default: TRUE)
#' @param legend_point_size Legend point size (default: 3)
#' @param save_path Output file path (optional)
#' @param width Plot width in inches (default: 12)
#' @param height Plot height in inches (default: 10)
#' @return ggplot object and saves plot to file if save_path provided
#' @export
visualize_spatial_network <- function(df, edges,
                                      celltype_col = "celltype",
                                      x_col = "X",
                                      y_col = "Y",
                                      edge_mode = c("all", "top", "random"),
                                      top_n = 1000,
                                      max_edges = 1000,
                                      point_size = 1.5,
                                      point_alpha = 0.7,
                                      edge_size_range = c(0.3, 2.0),
                                      edge_alpha_range = c(0.3, 0.9),
                                      edge_color = "darkred",
                                      show_points = TRUE,
                                      legend_point_size = 3,
                                      save_path = NULL,
                                      width = 12,
                                      height = 10) {

  edge_mode <- match.arg(edge_mode)

  # Convert to data.table if needed
  if (!data.table::is.data.table(df)) {
    df <- data.table::as.data.table(df)
  }
  if (!data.table::is.data.table(edges)) {
    edges <- data.table::as.data.table(edges)
  }

  # Prepare coordinates
  edges$from <- as.character(edges$from)
  edges$to <- as.character(edges$to)
  df$Cell_ID <- as.character(df$Cell_ID)

  df_coords <- df[, .(Cell_ID, get(x_col), get(y_col))]
  data.table::setnames(df_coords, c("Cell_ID", "X", "Y"))

  edges_with_coords <- merge(edges, df_coords, by.x = "from", by.y = "Cell_ID")
  data.table::setnames(edges_with_coords, c("X", "Y"), c("from_x", "from_y"))

  edges_with_coords <- merge(edges_with_coords, df_coords, by.x = "to", by.y = "Cell_ID")
  data.table::setnames(edges_with_coords, c("X", "Y"), c("to_x", "to_y"))

  # Calculate edge strength
  if ("dist" %in% names(edges_with_coords)) {
    edges_with_coords[, strength := 1 / dist]
  } else {
    edges_with_coords[, dist := sqrt((to_x - from_x)^2 + (to_y - from_y)^2)]
    edges_with_coords[, strength := 1 / dist]
  }

  # Filter edges based on mode
  if (edge_mode == "top") {
    if (nrow(edges_with_coords) > top_n) {
      edges_with_coords <- edges_with_coords[order(-strength)][1:top_n]
      message("Keep top ", top_n, " edges by strength")
    }
  } else if (edge_mode == "random") {
    if (nrow(edges_with_coords) > max_edges) {
      set.seed(123)
      edges_with_coords <- edges_with_coords[sample(.N, max_edges)]
      message("Randomly sample ", max_edges, " edges")
    }
  } else {
    message("Show all ", nrow(edges_with_coords), " edges")
  }

  # Create base plot
  p <- ggplot2::ggplot() +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Spatial Interaction Network") +
    ggplot2::coord_fixed() +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "white", color = "grey90"),
      panel.grid.major = ggplot2::element_line(color = "grey95", linewidth = 0.2),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14)
    )

  # Add points if requested
  if (show_points) {
    cell_types <- unique(df[[celltype_col]])
    n_types <- length(cell_types)
    palette <- get_color_palette(n_types)
    names(palette) <- cell_types

    p <- p + ggplot2::geom_point(
      data = df,
      ggplot2::aes_string(x = x_col, y = y_col, color = celltype_col),
      size = point_size,
      alpha = point_alpha
    ) +
      ggplot2::scale_color_manual(
        name = "Cell Type",
        values = palette,
        guide = ggplot2::guide_legend(
          override.aes = list(
            size = legend_point_size,
            alpha = 1
          )
        )
      )
  }

  # Add edges
  p <- p + ggplot2::geom_segment(
    data = edges_with_coords,
    ggplot2::aes(x = from_x, y = from_y,
                 xend = to_x, yend = to_y,
                 size = strength, alpha = strength),
    color = edge_color,
    show.legend = FALSE
  ) +
    ggplot2::scale_size_continuous(range = edge_size_range) +
    ggplot2::scale_alpha_continuous(range = edge_alpha_range)

  # Save plot if path provided
  if (!is.null(save_path)) {
    ggplot2::ggsave(save_path, p, width = width, height = height)
    message("Spatial network plot saved to: ", save_path)
  }

  return(p)
}

#' Visualize spatial interaction network using graph layout
#'
#' @param network An igraph network object from analyze_spatial_interactions()
#' @param node_size_range Range of node sizes (default: c(5, 15))
#' @param edge_size_range Range of edge sizes (default: c(0.5, 3))
#' @param label_size Text label size (default: 4)
#' @param show_labels Whether to show node labels (default: TRUE)
#' @param max_nodes Maximum number of nodes to display (default: 50)
#' @param save_path Output file path (optional)
#' @param width Plot width in inches (default: 12)
#' @param height Plot height in inches (default: 10)
#' @param layout Network layout algorithm (default: "fr")
#' @return ggraph object and saves plot to file if save_path provided
#' @export
visualize_interaction_network <- function(network,
                                          node_size_range = c(5, 15),
                                          edge_size_range = c(0.5, 3),
                                          label_size = 4,
                                          show_labels = TRUE,
                                          max_nodes = 50,
                                          save_path = NULL,
                                          width = 12,
                                          height = 10,
                                          layout = "fr") {

  # Check if network is valid
  if (igraph::vcount(network) == 0) {
    warning("Empty network provided")
    empty_plot <- ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = "Empty Network") +
      ggplot2::theme_void()

    if (!is.null(save_path)) {
      ggplot2::ggsave(save_path, empty_plot, width = width, height = height)
      message("Empty network plot saved to: ", save_path)
    }
    return(empty_plot)
  }

  # Calculate node metrics
  igraph::V(network)$degree <- igraph::degree(network)
  igraph::V(network)$strength <- igraph::strength(network)

  # Get cell types
  cell_types <- unique(igraph::V(network)$name)
  n_types <- length(cell_types)

  # Create scalable color scheme
  if (n_types <= 9) {
    color_scale <- ggplot2::scale_color_brewer(
      palette = "Set1",
      name = "Cell Type"
    )
  } else if (n_types <= 12) {
    color_scale <- ggplot2::scale_color_brewer(
      palette = "Paired",
      name = "Cell Type"
    )
  } else {
    color_scale <- ggplot2::scale_color_hue(
      h = c(0, 360) + 15,
      c = 100,
      l = 65,
      h.start = 0,
      direction = 1,
      name = "Cell Type"
    )
  }

  # Simplify large networks
  if (igraph::vcount(network) > max_nodes) {
    node_importance <- igraph::strength(network)
    top_nodes <- order(node_importance, decreasing = TRUE)[1:max_nodes]
    network <- igraph::induced_subgraph(network, top_nodes)
    message("Network too large, filtered to top ", max_nodes, " important nodes")
  }

  # Create network plot
  plot <- ggraph::ggraph(network, layout = layout) +
    # 1. Edges
    ggraph::geom_edge_link(
      ggplot2::aes(width = weight, alpha = weight),
      color = "grey50",
      show.legend = TRUE
    ) +
    ggraph::scale_edge_width_continuous(
      range = edge_size_range,
      name = "Interaction\nStrength"
    ) +
    ggraph::scale_edge_alpha_continuous(
      range = c(0.3, 0.8),
      guide = "none"
    ) +
    # 2. Nodes
    ggraph::geom_node_point(
      ggplot2::aes(size = strength, color = name),
      alpha = 0.85
    ) +
    ggplot2::scale_size_continuous(
      range = node_size_range,
      name = "Interaction\nStrength"
    ) +
    color_scale +
    # 3. Labels
    {
      if (show_labels) {
        ggraph::geom_node_text(
          ggplot2::aes(label = name),
          size = label_size,
          repel = TRUE,
          bg.color = "white",
          bg.r = 0.15,
          max.overlaps = 100
        )
      }
    } +
    # 4. Titles and theme
    ggplot2::labs(title = "Spatial Cell Interaction Network") +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
    )

  # Save plot if path provided
  if (!is.null(save_path)) {
    ggplot2::ggsave(save_path, plot, width = width, height = height)
    message("Interaction network plot saved to: ", save_path)
  }

  return(plot)
}

#' Visualize spatial data using Voronoi diagrams
#'
#' @param df Spatial data with coordinates and annotations
#' @param x_col X coordinate column name (default: "X")
#' @param y_col Y coordinate column name (default: "Y")
#' @param coloring Coloring method: "celltype" or "neighborhood" (default: "celltype")
#' @param highlight_cluster Specific cluster to highlight (optional)
#' @param celltype_col Cell type column name (default: "annotation")
#' @param neighborhood_col Neighborhood cluster column name (default: "Neighborhood_Cluster")
#' @param background_color Background color for non-highlighted cells (default: "gray90")
#' @param highlight_alpha Alpha for highlighted cells (default: 0.9)
#' @param background_alpha Alpha for background cells (default: 0.3)
#' @param celltype_palette Custom color palette for cell types (optional)
#' @param show_composition Show cell type composition in highlight (default: TRUE)
#' @param save_path Output file path (optional)
#' @param width Plot width in inches (default: 12)
#' @param height Plot height in inches (default: 10)
#' @return ggplot object and saves plot to file if save_path provided
#' @export
visualize_voronoi <- function(df,
                              x_col = "X",
                              y_col = "Y",
                              coloring = c("celltype", "neighborhood"),
                              highlight_cluster = NULL,
                              celltype_col = "annotation",
                              neighborhood_col = "Neighborhood_Cluster",
                              background_color = "gray90",
                              highlight_alpha = 0.9,
                              background_alpha = 0.3,
                              celltype_palette = NULL,
                              show_composition = TRUE,
                              save_path = NULL,
                              width = 12,
                              height = 10) {

  # Check required packages
  required_packages <- c("deldir", "ggplot2", "purrr", "dplyr")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Please install package: ", pkg, " - install.packages('", pkg, "')")
    }
  }

  # Validate parameters
  coloring <- match.arg(coloring)

  # Check if required columns exist
  if (!x_col %in% names(df)) stop("Column '", x_col, "' not found in data frame")
  if (!y_col %in% names(df)) stop("Column '", y_col, "' not found in data frame")
  if (!celltype_col %in% names(df)) stop("Column '", celltype_col, "' not found in data frame")
  if (!neighborhood_col %in% names(df)) stop("Column '", neighborhood_col, "' not found in data frame")

  # Create local column copies
  df <- df %>%
    dplyr::mutate(
      plot_X = .data[[x_col]],
      plot_Y = .data[[y_col]],
      plot_celltype = as.factor(.data[[celltype_col]]),
      plot_neighborhood = as.factor(.data[[neighborhood_col]])
    )

  # Create Voronoi tessellation
  voronoi <- deldir::deldir(df$plot_X, df$plot_Y)
  tiles <- deldir::tile.list(voronoi)

  # Create polygon data
  polygons <- purrr::map_df(seq_along(tiles), function(i) {
    tile <- tiles[[i]]
    data.frame(
      x = tile$x,
      y = tile$y,
      cell_id = i,
      celltype = df$plot_celltype[i],
      neighborhood = df$plot_neighborhood[i]
    )
  })

  # Set coloring method
  if (coloring == "celltype") {
    polygons$color_group <- polygons$celltype
    legend_title <- "Cell Type"

    # Use custom or default palette
    if (is.null(celltype_palette)) {
      n_celltypes <- length(unique(polygons$color_group))
      color_palette <- get_color_palette(n_celltypes)
    } else {
      color_palette <- celltype_palette
    }
  } else {
    polygons$color_group <- polygons$neighborhood
    legend_title <- "Neighborhood Cluster"

    # Generate new color scheme for neighborhood types
    n_clusters <- length(unique(polygons$color_group))
    color_palette <- scales::hue_pal()(n_clusters)
  }

  # Handle highlighting logic
  if (!is.null(highlight_cluster) && coloring == "neighborhood") {
    # Add highlight identifier column
    polygons$highlight <- ifelse(polygons$neighborhood == highlight_cluster,
                                 "Highlighted", "Background")

    # Get cell type color scheme
    if (is.null(celltype_palette)) {
      n_celltypes <- length(unique(polygons$celltype))
      celltype_palette <- get_color_palette(n_celltypes)
      names(celltype_palette) <- levels(polygons$celltype)
    }

    # Create highlight color scheme
    polygons$fill_color <- ifelse(
      polygons$highlight == "Highlighted",
      as.character(polygons$celltype),
      background_color
    )

    polygons$border_color <- ifelse(polygons$highlight == "Highlighted", "white", "grey75")

    # Create plot object
    p <- ggplot2::ggplot(polygons, ggplot2::aes(x = x, y = y, group = cell_id)) +
      ggplot2::geom_polygon(
        ggplot2::aes(fill = fill_color, alpha = highlight, color = border_color),
        size = 0.2
      ) +
      ggplot2::scale_alpha_manual(
        values = c("Highlighted" = highlight_alpha, "Background" = background_alpha),
        guide = "none"
      ) +
      ggplot2::scale_color_identity() +
      ggplot2::theme_void() +
      ggplot2::coord_fixed() +
      ggplot2::labs(
        title = paste("Voronoi Diagram - Highlighting", highlight_cluster),
        subtitle = if (show_composition) {
          paste("Cell composition in", highlight_cluster)
        } else {
          NULL
        }
      )

    # Add cell type legend (only for highlighted area)
    if (show_composition) {
      # Get cell types in highlighted area
      highlight_celltypes <- unique(polygons$celltype[polygons$highlight == "Highlighted"])
      highlight_palette <- celltype_palette[names(celltype_palette) %in% as.character(highlight_celltypes)]

      p <- p +
        ggplot2::scale_fill_manual(
          values = c(highlight_palette, "Background" = background_color),
          breaks = names(highlight_palette),
          name = "Cell Types"
        ) +
        ggplot2::guides(
          fill = ggplot2::guide_legend(
            override.aes = list(alpha = 1, color = NA),
            ncol = min(3, length(highlight_palette))
          )
        )
    } else {
      p <- p + ggplot2::scale_fill_identity()
    }
  } else {
    # Standard plot without highlighting
    if (coloring == "celltype") {
      # Ensure cell type colors match original function
      if (is.null(celltype_palette)) {
        color_palette <- get_color_palette(length(unique(polygons$color_group)))
      } else {
        color_palette <- celltype_palette
      }
    }

    polygons$border_color <- "white"

    p <- ggplot2::ggplot(polygons, ggplot2::aes(x = x, y = y, group = cell_id, fill = color_group)) +
      ggplot2::geom_polygon(
        ggplot2::aes(color = border_color),
        alpha = 0.8, size = 0.2
      ) +
      ggplot2::scale_fill_manual(values = color_palette, name = legend_title) +
      ggplot2::scale_color_identity() +
      ggplot2::theme_void() +
      ggplot2::coord_fixed() +
      ggplot2::labs(title = paste("Voronoi Diagram -", legend_title))
  }

  # Calculate plot boundaries with margin
  x_range <- range(df$plot_X, na.rm = TRUE)
  y_range <- range(df$plot_Y, na.rm = TRUE)

  # Add 5% margin
  x_margin <- (x_range[2] - x_range[1]) * 0.05
  y_margin <- (y_range[2] - y_range[1]) * 0.05

  p <- p +
    ggplot2::xlim(x_range[1] - x_margin, x_range[2] + x_margin) +
    ggplot2::ylim(y_range[1] - y_margin, y_range[2] + y_margin)

  # Save plot if path provided
  if (!is.null(save_path)) {
    ggplot2::ggsave(save_path, p, width = width, height = height)
    message("Voronoi diagram saved to: ", save_path)
  }

  return(p)
}
