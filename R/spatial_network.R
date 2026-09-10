# R/spatial_network.R

#' Prepare spatial data for network analysis
#'
#' @param df Data frame containing spatial data
#' @param cell_id_col Column name for cell IDs (default: NULL, auto-detect)
#' @param x_col Column name for X coordinates (default: "X")
#' @param y_col Column name for Y coordinates (default: "Y")
#' @param celltype_col Column name for cell types (default: "celltype")
#' @return data.table with standardized structure for spatial analysis
#' @examples
#' df <- Sphinx:::.sphinx_example_df(30)
#' df <- prepare_data(df)
#' head(df[, .(Cell_ID, X, Y, celltype)])
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
#' @examples
#' df <- prepare_data(Sphinx:::.sphinx_example_df(50))
#' params <- calculate_optimal_radius(df)
#' names(params)
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

#' Calculate an adaptive window size from cell density
#'
#' Iteratively adjusts window side length so that the mean number of cells per
#' window falls between `min_cells` and `max_cells`.
#'
#' @param df Spatial data with coordinates
#' @param x_col Column name for X coordinates (default: "X")
#' @param y_col Column name for Y coordinates (default: "Y")
#' @param min_cells Target minimum mean cells per window (default: 30)
#' @param max_cells Target maximum mean cells per window (default: 50)
#' @param range_divisor Divisor for initial window size from spatial range (default: 100)
#' @param max_iter Maximum adjustment iterations (default: 100)
#' @return List with window_size, sliding_step, mean_cells, x_range, y_range
#' @examples
#' df <- prepare_data(Sphinx:::.sphinx_example_df(50))
#' ws <- calculate_optimal_window_size(df)
#' ws$window_size
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
  res <- .calc_adaptive_window_size(
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
#' @param celltype_col Column name for cell types (default: "celltype")
#' @return List with distance matrix and distributions
#' @examples
#' df <- prepare_data(Sphinx:::.sphinx_example_df(40))
#' dist_res <- calculate_celltype_distances(df, celltype_col = "celltype")
#' dim(dist_res$distance_matrix)
#' @export
calculate_celltype_distances <- function(df,
                                         x_col = "X",
                                         y_col = "Y",
                                         celltype_col = "celltype") {
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
#' selects the method from within-sample spatial metrics using an empirical
#' decision tree (multi-cohort network evaluation; see `.select_method()`):
#' baseline knn(k=10, mutual); radius when types are strongly clustered and
#' density is fairly uniform; delaunay for continuous / large / heterogeneous
#' layouts (optional length gate); window only for ordered tiled-like FOVs.
#'
#' By default there is **no** distance gate for knn / delaunay / window (same
#' idea as leaving `max_dist` unset), except that `method = "auto"` may
#' recommend `max_edge_length = "auto"` for delaunay. Optionally set
#' `max_dist` (absolute), a numeric `max_edge_length`, or
#' `max_edge_length = "auto"` (upper quantile of per-cell k-th NN distances)
#' to drop long edges. Each cell can also be limited to its `n_neighbors` /
#' `max_degree` nearest contacts.
#'
#' @param df data.frame/data.table with columns: X, Y, Cell_ID; optional cell type column
#' @param method one of c("auto","radius","knn","delaunay","window")
#' @param k **deprecated** neighbor count. Use `n_neighbors`. Clustering `k`
#'   belongs to `cluster_neighborhoods()`, not this function.
#' @param radius numeric, search radius for radius method (auto if NULL)
#' @param window_size numeric, sliding-window side length for window method
#'   (auto if NULL). Stereopy CCD: target 30-50 cells per window.
#' @param window_slide_step numeric, stride of the sliding window. `NULL`
#'   (default) uses half of `window_size`, as in Stereopy CCD.
#' @param max_edge_length distance gate. `NULL` / `Inf` (default): no gate for
#'   knn/delaunay/window. Numeric: keep edges with `dist <=` this value.
#'   `"auto"`: use `edge_length_quantile` of per-cell k-th NN distances
#'   (scaled by `edge_length_mult`). For `method = "radius"`, the search
#'   radius itself is the gate.
#' @param max_dist numeric absolute distance cap in the same units as X/Y.
#'   If set, overrides `max_edge_length`.
#' @param edge_length_mult numeric multiplier for `max_edge_length = "auto"`
#'   (default 1.0). Ignored unless auto gating is requested.
#' @param edge_length_quantile quantile of per-cell k-th NN distances when
#'   `max_edge_length = "auto"` (default 0.95).
#' @param n_neighbors integer number of spatial neighbors (knn/window search
#'   rank and post-filter degree cap). Default `10` for knn/window (aligned
#'   with `.default_nn_rank()` and the auto baseline). For
#'   `method = "radius"` or `"delaunay"`, omitted/`NULL` means no degree cap
#'   (`Inf`: keep the ball graph / triangulation). Alias of `max_degree`.
#'   This is **not** the clustering `k` used by `cluster_neighborhoods()`.
#' @param max_degree integer. Alias of `n_neighbors`. Used only when
#'   `n_neighbors` is not supplied. If both are set, `n_neighbors` wins.
#'   Set `Inf` to disable degree capping.
#' @param require_mutual logical. If TRUE, keep an edge only when both cells list
#'   each other among their `n_neighbors` nearest contacts. Default `TRUE` for
#'   knn/window (auto baseline). For `method = "radius"` defaults to `FALSE`
#'   unless set explicitly. Delaunay edges are already mutual; the flag only
#'   matters if a degree cap is also set.
#' @param max_edges deprecated; ignored. Edges are filtered by distance and
#'   per-cell degree only (no global random downsampling).
#' @param n_cores integer parallel workers for knn / radius construction. NULL
#'   uses `SPHINX_N_CORES` if set, otherwise up to 8 detected cores (fork on Unix).
#'   Also passed as `num.threads` to kNN search when using a single chunk.
#' @param celltype_col character, column with cell-type labels (optional)
#' @param verbose logical, print decisions and key metrics
#' @return data.table with columns: from, to, dist (and optional context columns)
#' @examples
#' df <- prepare_data(Sphinx:::.sphinx_example_df(40))
#' edges <- build_spatial_network(df, method = "knn", n_neighbors = 5, verbose = FALSE)
#' # Optional absolute distance gate (same units as X/Y):
#' edges2 <- build_spatial_network(
#'   df, method = "knn", n_neighbors = 10, max_dist = 40, verbose = FALSE
#' )
#' # Optional adaptive gate (quantile of k-th NN distances):
#' edges3 <- build_spatial_network(
#'   df, method = "knn", n_neighbors = 10, max_edge_length = "auto", verbose = FALSE
#' )
#' head(edges)
#' @export
build_spatial_network <- function(
    df,
    method = "auto",
    k = NULL,
    radius = NULL,
    window_size = NULL,
    window_slide_step = NULL,
    max_edge_length = NULL,
    max_dist = NULL,
    edge_length_mult = 1.0,
    edge_length_quantile = 0.95,
    n_neighbors = 10L,
    max_degree = NULL,
    require_mutual = TRUE,
    max_edges = NULL,
    n_cores = NULL,
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
  contact_1nn <- metrics$median_1nn_dist
  if (!is.finite(contact_1nn) || contact_1nn <= 0) {
    contact_1nn <- metrics$mean_nn_dist
  }

  # n_neighbors = spatial neighbor count (knn/window rank + degree cap).
  # k is clustering-only (cluster_neighborhoods). Old calls that passed k as
  # neighbor count are still accepted with a warning.
  # Default n_neighbors=10 must not overwrite an explicit max_degree.
  user_set_n_neighbors <- !missing(n_neighbors) && !is.null(n_neighbors)
  user_set_max_degree <- !missing(max_degree) && !is.null(max_degree)
  user_set_degree <- user_set_n_neighbors || user_set_max_degree
  user_set_mutual <- !missing(require_mutual)
  user_set_mel <- !missing(max_edge_length)
  user_set_radius <- !missing(radius) && !is.null(radius)
  if (!is.null(k)) {
    warning(
      "`k` in build_spatial_network() is deprecated for neighbor count. ",
      "Use n_neighbors. Clustering k belongs to cluster_neighborhoods().",
      call. = FALSE
    )
    if (!user_set_degree) {
      n_neighbors <- as.integer(k)
      user_set_n_neighbors <- TRUE
      user_set_degree <- TRUE
    }
  }

  if (user_set_n_neighbors) {
    max_degree <- n_neighbors
  } else if (!user_set_max_degree) {
    max_degree <- if (!is.null(n_neighbors)) n_neighbors else .default_nn_rank()
  }
  if (is.finite(max_degree)) {
    max_degree <- as.integer(max_degree)
    if (max_degree < 1L) {
      stop("n_neighbors / max_degree must be >= 1, or Inf to disable.")
    }
  }
  # Neighbor rank used by knn / window builders (not clustering).
  nn_k <- if (is.finite(max_degree)) max_degree else .default_nn_rank()

  # ---- Auto method selection (returns method + recommended params) ----
  chosen <- method
  auto_params <- list()
  auto_reason <- NULL
  if (method == "auto") {
    sel <- .select_method(metrics, nrow(df))
    chosen <- sel$method
    auto_params <- sel$params
    auto_reason <- sel$reason
  }

  if (is.null(edge_length_mult) || !is.finite(edge_length_mult) || edge_length_mult <= 0) {
    edge_length_mult <- 1.0
  }
  if (is.null(edge_length_quantile) || !is.finite(edge_length_quantile) ||
      edge_length_quantile <= 0 || edge_length_quantile > 1) {
    edge_length_quantile <- 0.95
  }
  if (!is.null(auto_params$n_neighbors) && is.finite(auto_params$n_neighbors) &&
      !user_set_degree) {
    nn_k <- as.integer(auto_params$n_neighbors)
    max_degree <- nn_k
  }
  if (!is.null(auto_params$require_mutual) && !user_set_mutual) {
    require_mutual <- isTRUE(auto_params$require_mutual)
  }
  if (!is.null(auto_params$max_edge_length) && !user_set_mel) {
    max_edge_length <- auto_params$max_edge_length
  }
  if (!is.null(auto_params$radius) && !user_set_radius) {
    radius <- auto_params$radius
  }

  # Distance gate (optional). Default: none for knn/delaunay/window.
  # Priority: max_dist > max_edge_length ("auto" | numeric) > Inf.
  mel_auto <- FALSE
  if (!is.null(max_dist) && is.finite(max_dist)) {
    if (!is.null(max_edge_length) && !identical(max_edge_length, "auto") &&
        is.finite(suppressWarnings(as.numeric(max_edge_length)[1L])) &&
        !isTRUE(all.equal(as.numeric(max_dist), as.numeric(max_edge_length)))) {
      warning(
        "Both max_dist and max_edge_length were set; using max_dist=",
        max_dist, " (absolute)."
      )
    }
    max_edge_length <- as.numeric(max_dist)
  } else if (is.null(max_edge_length)) {
    if (identical(chosen, "radius")) {
      # Filled from search radius below.
      max_edge_length <- NA_real_
    } else {
      max_edge_length <- Inf
    }
  } else if (is.character(max_edge_length) &&
             tolower(max_edge_length[[1L]]) %in% c("auto", "quantile")) {
    mel_auto <- TRUE
    k_for_gate <- if (identical(chosen, "delaunay")) {
      6L
    } else if (identical(chosen, "radius")) {
      .default_nn_rank()
    } else {
      as.integer(nn_k)
    }
    max_edge_length <- .default_max_edge_length(
      as.matrix(df[, .(X, Y)]),
      k = k_for_gate,
      quantile = edge_length_quantile,
      multiplier = edge_length_mult
    )
  } else {
    max_edge_length <- suppressWarnings(as.numeric(max_edge_length)[1L])
    if (!is.finite(max_edge_length)) max_edge_length <- Inf
  }

  if (is.null(radius)) {
    radius <- if (!is.null(auto_params$radius)) {
      auto_params$radius
    } else if (identical(chosen, "radius")) {
      .default_search_radius(contact_1nn, multiplier = 3)
    } else {
      # Metadata only for non-radius methods when no distance gate is set.
      if (is.finite(max_edge_length)) {
        max_edge_length
      } else {
        .default_search_radius(contact_1nn, multiplier = 3)
      }
    }
  }
  if (chosen == "radius") {
    if (!is.finite(max_edge_length) || is.na(max_edge_length)) {
      # Radius graph is already a ball of radius `radius`.
      max_edge_length <- radius
    } else {
      # Optional tighter numeric / auto / max_dist gate on top of the ball.
      max_edge_length <- min(as.numeric(max_edge_length), as.numeric(radius))
    }
  }

  # Radius: all pairs with dist <= R. Delaunay: keep the triangulation.
  # Apply degree / mutual caps only if the user set them.
  if (identical(chosen, "radius") || identical(chosen, "delaunay")) {
    if (!user_set_degree) max_degree <- Inf
    if (!user_set_mutual && identical(chosen, "radius")) require_mutual <- FALSE
  }

  # window_size / window_slide_step: Stereopy CCD sliding window (d, s=d/2).
  if (is.null(window_size)) {
    if (!is.null(auto_params$window_size)) {
      window_size <- auto_params$window_size
      if (is.null(window_slide_step)) {
        window_slide_step <- auto_params$window_slide_step
      }
    } else if (identical(chosen, "window") || identical(method, "auto")) {
      opt_win <- .calc_adaptive_window_size(df$X, df$Y)
      window_size <- opt_win$window_size
      if (is.null(window_slide_step)) {
        window_slide_step <- opt_win$sliding_step
      }
    } else {
      xr <- diff(range(df$X, na.rm = TRUE))
      yr <- diff(range(df$Y, na.rm = TRUE))
      window_size <- max(10, min(xr, yr) / 40)
      if (is.null(window_slide_step)) window_slide_step <- window_size / 2
    }
  } else if (is.null(window_slide_step)) {
    window_slide_step <- as.numeric(window_size) / 2
  }
  if (is.null(window_slide_step) || !is.finite(window_slide_step) ||
      window_slide_step <= 0) {
    window_slide_step <- as.numeric(window_size) / 2
  }

  if (verbose) {
    if (identical(method, "auto") && !is.null(auto_reason) && nzchar(auto_reason)) {
      message("auto -> ", chosen, " | ", auto_reason)
    }
    mel_msg <- if (is.finite(max_edge_length)) {
      paste0(
        round(max_edge_length, 3),
        if (isTRUE(mel_auto)) " (auto)" else if (!is.null(max_dist) && is.finite(max_dist)) " (max_dist)" else ""
      )
    } else {
      "Inf (none)"
    }
    msg <- paste0(
      "Method=", chosen,
      " | n=", nrow(df),
      " | radius=", round(radius,3),
      " | window_size=", round(window_size,3),
      if (identical(chosen, "window")) {
        paste0(" | window_slide_step=", round(window_slide_step, 3))
      } else {
        ""
      },
      " | max_edge_length=", mel_msg,
      " | median1NN=", round(contact_1nn,3),
      " | meanNN10=", round(metrics$mean_nn_dist,3),
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
    knn      = .edges_knn(df, k = nn_k, max_edge_length = max_edge_length,
                          n_cores = n_cores),
    radius   = .edges_radius(df, radius = radius, n_cores = n_cores),
    delaunay = .edges_delaunay(df, max_edge_length = max_edge_length, k = nn_k),
    window   = .edges_window(
      df, k = nn_k, tile = window_size,
      sliding_step = window_slide_step,
      max_edge_length = max_edge_length
    ),
    stop("Unknown method: ", chosen)
  )

  if (!is.null(max_edges) && is.finite(max_edges)) {
    warning(
      "max_edges is deprecated and ignored; edges are filtered by ",
      "max_edge_length / max_dist + n_neighbors (no random downsampling)."
    )
  }

  edges <- .filter_edges_biological(
    edges,
    max_edge_length = max_edge_length,
    max_degree = max_degree,
    require_mutual = isTRUE(require_mutual),
    verbose = verbose
  )

  if (nrow(edges) == 0L) {
    warning("No edges produced; try increasing n_neighbors/radius/window_size/max_edge_length/max_dist.")
    return(edges)
  }

  # ---- Optional context (cheap) ----
  if (!is.null(celltype_col) && celltype_col %in% names(df)) {
    df_small <- df[, .(Cell_ID, celltype = get(celltype_col))]
    edges[df_small, on = .(from = Cell_ID), from_type := i.celltype]
    edges[df_small, on = .(to   = Cell_ID), to_type   := i.celltype]
    edges[, same_type := (from_type == to_type)]
  }

  data.table::setattr(edges, "metrics", metrics)
  data.table::setattr(edges, "method", chosen)
  data.table::setattr(edges, "parameters", list(
    n_neighbors = max_degree,
    radius = radius,
    window_size = window_size,
    window_slide_step = window_slide_step,
    max_edge_length = max_edge_length,
    max_dist = if (!is.null(max_dist) && is.finite(max_dist)) as.numeric(max_dist) else NA_real_,
    edge_length_mult = edge_length_mult,
    edge_length_quantile = edge_length_quantile,
    distance_gate = if (!is.null(max_dist) && is.finite(max_dist)) {
      "max_dist"
    } else if (isTRUE(mel_auto)) {
      "auto_quantile"
    } else if (is.finite(max_edge_length) && !identical(chosen, "radius")) {
      "max_edge_length"
    } else if (identical(chosen, "radius")) {
      "radius"
    } else {
      "none"
    },
    contact_1nn_median = contact_1nn,
    max_degree = max_degree,
    require_mutual = isTRUE(require_mutual),
    n_cores = .sphinx_default_cores(n_cores),
    auto_requested = identical(method, "auto"),
    auto_reason = auto_reason
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

#' Default search radius from local spacing (median 1-NN x multiplier)
#' @param mean_nn_dist local spacing distance
#' @param multiplier scale factor (default: 3 ~ a few cell spacings)
#' @return numeric search radius
#' @keywords internal
.default_search_radius <- function(mean_nn_dist, multiplier = 3) {
  mean_nn_dist * multiplier
}

#' Default maximum edge length from the k-th NN distance distribution
#'
#' Uses the upper quantile of each cell's distance to its k-th nearest neighbor.
#' This matches the natural length scale of a kNN / local triangulation graph:
#' typical neighborhoods are kept, while only sparse-region outliers are cut.
#' (A cap of ~1.5xmedian(1-NN) is far too tight once k >= 5-10.)
#'
#' @param coords numeric matrix with X/Y columns (or a single spacing scalar for
#'   backward-compatible callers - then treated as median 1-NN x 3)
#' @param k neighbor rank used for the distance distribution (default 10)
#' @param quantile upper quantile in (0, 1] (default 0.95)
#' @param multiplier optional scale on the quantile (default 1)
#' @return numeric max edge length
#' @keywords internal
.default_max_edge_length <- function(coords,
                                     k = .default_nn_rank(),
                                     quantile = 0.95,
                                     multiplier = 1) {
  if (is.null(multiplier) || !is.finite(multiplier) || multiplier <= 0) {
    multiplier <- 1
  }
  # Backward compatible: scalar spacing -> contact-scale fallback
  if (is.numeric(coords) && length(coords) == 1L) {
    return(as.numeric(coords) * 3 * multiplier)
  }
  k <- as.integer(k)
  if (!is.finite(k) || k < 1L) k <- .default_nn_rank()
  if (is.null(quantile) || !is.finite(quantile) || quantile <= 0 || quantile > 1) {
    quantile <- 0.95
  }
  dk <- .nn_dist_at_rank(coords, k_nn = k)
  thr <- as.numeric(stats::quantile(dk, probs = quantile, na.rm = TRUE, names = FALSE))
  if (!is.finite(thr) || thr <= 0) {
    thr <- stats::median(dk, na.rm = TRUE)
  }
  thr * multiplier
}

#' Filter edges by optional distance and degree rules
#'
#' Rules (no global random downsampling):
#' 1. Distance (optional): if `max_edge_length` is finite, keep edges with
#'    dist <= max_edge_length.
#' 2. Neighbor cap: each cell keeps at most `max_degree` / `n_neighbors`
#'    nearest contacts (default 10 for non-radius methods).
#' 3. Optional mutual filter: if `require_mutual`, keep an edge only when both
#'    endpoints nominate each other (reciprocal contact).
#'
#' @param edges data.table with from, to, dist
#' @param max_edge_length numeric distance cap (Inf to skip)
#' @param max_degree integer per-cell neighbor cap (Inf to skip); aka n_neighbors
#' @param require_mutual logical; mutual nearest neighbors if TRUE
#' @param verbose logical
#' @return filtered data.table with from, to, dist
#' @keywords internal
.filter_edges_biological <- function(edges,
                                     max_edge_length = Inf,
                                     max_degree = 10L,
                                     require_mutual = TRUE,
                                     verbose = FALSE) {
  if (is.null(edges) || !nrow(edges)) return(edges)
  if (!data.table::is.data.table(edges)) {
    edges <- data.table::as.data.table(edges)
  }
  n0 <- nrow(edges)

  if (is.finite(max_edge_length)) {
    edges <- edges[is.finite(dist) & dist <= max_edge_length]
  }
  n1 <- nrow(edges)
  if (!nrow(edges)) {
    if (verbose) {
      message("Biological edge filter: ", n0, " -> 0 (all edges beyond contact distance)")
    }
    return(edges[, .(from, to, dist)])
  }

  if (is.null(max_degree) || !is.finite(max_degree) || max_degree < 1) {
    if (verbose && n0 != n1) {
      message(
        "Distance filter: ", n0, " -> ", n1,
        " (max_edge_length=", round(max_edge_length, 3), ")"
      )
    }
    return(edges[, .(from, to, dist)])
  }

  max_degree <- as.integer(max_degree)
  d <- data.table::rbindlist(list(
    edges[, .(cell = as.character(from), nbr = as.character(to), dist)],
    edges[, .(cell = as.character(to), nbr = as.character(from), dist)]
  ), use.names = TRUE)
  data.table::setorder(d, cell, dist)
  d[, ord := seq_len(.N), by = cell]
  d <- d[ord <= max_degree]
  d[, key := paste(pmin(cell, nbr), pmax(cell, nbr), sep = "\t")]
  key_n <- d[, .N, by = key]
  if (isTRUE(require_mutual)) {
    keep <- key_n[N >= 2L, key]
  } else {
    keep <- key_n$key
  }
  if (!length(keep)) {
    return(edges[0L, .(from, to, dist)])
  }
  out <- d[key %in% keep, .(dist = min(dist)), by = key]
  split_key <- strsplit(out$key, "\t", fixed = TRUE)
  out[, `:=`(
    from = vapply(split_key, `[[`, character(1), 1L),
    to = vapply(split_key, `[[`, character(1), 2L),
    key = NULL
  )]

  if (verbose) {
    mel <- if (is.finite(max_edge_length)) round(max_edge_length, 3) else "Inf"
    message(
      "Biological edge filter: ", n0, " -> ", nrow(out),
      " (max_edge_length=", mel,
      ", n_neighbors=", max_degree,
      ", require_mutual=", isTRUE(require_mutual), ")"
    )
  }
  out[, .(from, to, dist)]
}

#' Round to the nearest even integer (minimum 2)
#' @keywords internal
.round_to_even <- function(x) {
  r <- 2 * round(x / 2)
  max(2L, as.integer(r))
}

#' Number of grid windows covering the spatial extent
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

#' Adaptive window size by iterative density targeting
#' @keywords internal
.calc_adaptive_window_size <- function(x,
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
  median_1nn <- stats::median(d1, na.rm = TRUE)

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
    # mean distance to the 10th NN - used for density CV, NOT contact length
    mean_nn_dist = mean_nn,
    mean_nn_rank = k_nn,
    # median 1-NN - contact / radius scale (not the knn distance gate)
    median_1nn_dist = as.numeric(median_1nn),
    cv_nn_dist = cv_nn,
    clark_evans_R = clark_evans_R,
    cat_morans_I = cat_morans_I,
    het_gini_simpson = het_gs
  )
}

#' Select optimal method based on spatial metrics
#'
#' Empirical decision tree for `method = "auto"`, calibrated on a multi-cohort
#' within-sample composite (PAS, SI_expr, speed) across ~180 FOVs. Rules are
#' applied in order; first match wins:
#' \enumerate{
#'   \item Tiled / blocked-like FOV (ordered Clark-Evans R > 1.15 and uniform
#'     density CV < 0.3, n < 5e4) -> `window` with knn(k=10, mutual).
#'   \item Strong type clustering (`cat_morans_I >= 0.4`) **and** fairly uniform
#'     density (`cv_nn_dist < 0.3`) -> `radius` at ~3x median 1-NN.
#'   \item Continuous epithelium / few holes: large n (`>= 2e5`), heterogeneous
#'     density (`cv >= 0.6`), or space-filling irregular layout -> `delaunay`
#'     with optional `max_edge_length = "auto"`.
#'   \item Otherwise baseline -> `knn` with `n_neighbors = 10`, `require_mutual = TRUE`.
#' }
#' Scoring parallelism (`n_jobs`) is an evaluation concern, not encoded here;
#' builders still use `n_cores`.
#'
#' @param m spatial metrics
#' @param n number of cells
#' @return list with `method`, `params`, and `reason`
#' @keywords internal
.select_method <- function(m, n) {
  params <- list()
  contact <- m$median_1nn_dist
  if (is.null(contact) || !is.finite(contact) || contact <= 0) {
    contact <- m$mean_nn_dist
  }

  cv <- m$cv_nn_dist
  ce_r <- m$clark_evans_R
  morans <- m$cat_morans_I
  uniform_density <- !is.na(cv) && is.finite(cv) && cv < 0.3
  hetero_density <- !is.na(cv) && is.finite(cv) && cv >= 0.6
  ordered_layout <- !is.na(ce_r) && is.finite(ce_r) && ce_r > 1.15
  space_filling_irregular <- !is.na(ce_r) && is.finite(ce_r) && ce_r < 0.95 &&
    !is.na(cv) && is.finite(cv) && cv >= 0.35 && n >= 10000L

  # 1) Tiled / blocked FOV (narrow proxy: ordered + uniform, not huge)
  if (ordered_layout && uniform_density && n < 50000L) {
    params$n_neighbors <- .default_nn_rank()
    params$require_mutual <- TRUE
    return(list(
      method = "window",
      params = params,
      reason = sprintf(
        "tiled/blocked-like (CE_R=%.3f>1.15, CV=%.3f<0.3, n=%s<5e4)",
        ce_r, cv, format(n, scientific = FALSE)
      )
    ))
  }

  # 2) Strong type clustering + fairly uniform density -> radius (~3x median 1-NN)
  if (!is.na(morans) && is.finite(morans) && morans >= 0.4 && uniform_density) {
    params$radius <- .default_search_radius(contact, multiplier = 3)
    return(list(
      method = "radius",
      params = params,
      reason = sprintf(
        "type clustering + uniform density (MoranI=%.3f>=0.4, CV=%.3f<0.3; R~3xmedian1NN)",
        morans, cv
      )
    ))
  }

  if (is.na(morans) || !is.finite(morans)) {
    warning(
      "Categorical Moran's I unavailable; skipping radius branch ",
      "(requires MoranI>=0.4 with uniform density)."
    )
  }

  # 3) Continuous epithelium / few holes -> delaunay (+ optional length gate)
  continuous <- (n >= 200000L) ||
    (hetero_density && n >= 5000L) ||
    space_filling_irregular
  if (continuous) {
    params$max_edge_length <- "auto"
    why <- if (n >= 200000L) {
      sprintf("large continuous FOV (n=%s>=2e5)", format(n, scientific = FALSE))
    } else if (hetero_density && n >= 5000L) {
      sprintf("heterogeneous density (CV=%.3f>=0.6, n=%s)", cv, format(n, scientific = FALSE))
    } else {
      sprintf(
        "space-filling irregular layout (CE_R=%.3f, CV=%.3f, n=%s)",
        ce_r, cv, format(n, scientific = FALSE)
      )
    }
    return(list(
      method = "delaunay",
      params = params,
      reason = paste0(why, "; max_edge_length=auto")
    ))
  }

  # 4) Baseline: knn(k=10, mutual=TRUE)
  params$n_neighbors <- .default_nn_rank()
  params$require_mutual <- TRUE
  list(
    method = "knn",
    params = params,
    reason = sprintf(
      "baseline knn(k=%d, mutual=TRUE)%s",
      as.integer(params$n_neighbors),
      if (!is.na(morans) && is.finite(morans) && morans >= 0.4 && !uniform_density) {
        sprintf(" (MoranI=%.3f but CV=%.3f not uniform -> skip radius)", morans, cv)
      } else {
        ""
      }
    )
  )
}

#' Build edges using kNN method
#'
#' Exact k nearest neighbors via `BiocNeighbors::findKNN` (self excluded).
#' Parallelism uses `num.threads` inside the search (fork-safe). For large
#' point sets, query rows are processed in serial chunks via `subset=` so peak
#' memory stays bounded; chunk results concatenate to the same neighbors as a
#' single full call (no downsampling).
#'
#' @param df data.table with spatial data
#' @param k number of neighbors
#' @param max_edge_length optional maximum edge length filter
#' @param n_cores parallel workers / search threads (NULL = auto)
#' @return data.table with edges
#' @keywords internal
.edges_knn <- function(df, k, max_edge_length = NULL, n_cores = NULL) {
  n <- nrow(df)
  empty <- data.table::data.table(from = character(), to = character(), dist = numeric())
  if (n < 2L) return(empty)
  k <- as.integer(k)
  if (!is.finite(k) || k < 1L) {
    stop("k must be a positive integer for knn.")
  }
  k <- min(k, n - 1L)

  if (!requireNamespace("BiocNeighbors", quietly = TRUE)) {
    stop("Package 'BiocNeighbors' is required for knn. Please install it.")
  }

  coords <- as.matrix(df[, .(X, Y)])
  storage.mode(coords) <- "double"
  ids <- as.character(df$Cell_ID)
  n_cores <- .sphinx_default_cores(n_cores)

  # Serial chunks + internal threads: exact and avoids unsafe OpenMP+fork.
  # Chunk when n is large to bound intermediate index/distance matrices.
  chunk_size <- 50000L
  if (n <= chunk_size) {
    chunks <- list(seq_len(n))
  } else {
    n_chunks <- as.integer(ceiling(n / chunk_size))
    chunks <- parallel::splitIndices(n, n_chunks)
  }

  parts <- vector("list", length(chunks))
  for (ci in seq_along(chunks)) {
    idx <- chunks[[ci]]
    nn <- BiocNeighbors::findKNN(
      coords,
      k = k,
      get.index = TRUE,
      get.distance = TRUE,
      num.threads = n_cores,
      subset = idx
    )
    from_idx <- rep(idx, each = k)
    to_idx <- as.vector(t(nn$index))
    d <- as.vector(t(nn$distance))
    parts[[ci]] <- data.table::data.table(
      from = ids[from_idx],
      to = ids[to_idx],
      dist = d
    )
  }
  ed <- data.table::rbindlist(parts, use.names = TRUE)
  if (!nrow(ed)) return(empty)

  if (!is.null(max_edge_length) && is.finite(max_edge_length)) {
    ed <- ed[dist <= max_edge_length]
  }
  if (nrow(ed) == 0L) return(ed)
  ed[, key := ifelse(from < to, paste(from, to), paste(to, from))]
  ed <- ed[order(key, dist)][!duplicated(key)][, key := NULL][]
  ed
}

#' Build edges using radius method (exact ball graph)
#'
#' Connects every pair of cells whose Euclidean distance is <= `radius` using
#' `dbscan::frNN` as the primary search. Query points are split into chunks and
#' searched against the full coordinate set (`frNN(..., query=)`). Chunks run
#' in parallel via fork on Unix to lower peak memory and wall time; neighbors
#' match a single full `frNN` call after dropping self-hits from the query
#' interface. Undirected edges are emitted once (`j > i`).
#'
#' @param df data.table with spatial data
#' @param radius search radius
#' @param n_cores parallel workers (NULL = auto)
#' @return data.table with edges
#' @keywords internal
.edges_radius <- function(df, radius, n_cores = NULL) {
  n <- nrow(df)
  empty <- data.table::data.table(from = character(), to = character(), dist = numeric())
  if (n < 2L || is.null(radius) || !is.finite(radius) || radius <= 0) {
    return(empty)
  }
  if (!requireNamespace("dbscan", quietly = TRUE)) {
    stop("Package 'dbscan' is required for radius graphs. Please install it.")
  }

  old_dt_threads <- data.table::getDTthreads()
  on.exit(data.table::setDTthreads(old_dt_threads), add = TRUE)
  data.table::setDTthreads(1L)

  coords <- as.matrix(df[, .(X, Y)])
  storage.mode(coords) <- "double"
  ids <- as.character(df$Cell_ID)
  r <- as.numeric(radius)
  n_cores <- .sphinx_default_cores(n_cores)

  # Prefer more / smaller chunks under parallel to bound per-worker frNN lists.
  n_chunks <- if (n_cores <= 1L) {
    1L
  } else {
    max(as.integer(n_cores), min(n, as.integer(n_cores) * 4L))
  }
  n_chunks <- max(1L, min(n_chunks, n))
  if (n_chunks <= 1L) {
    chunks <- list(seq_len(n))
  } else {
    chunks <- parallel::splitIndices(n, n_chunks)
  }

  worker <- function(idx) {
    # query= returns neighbors in the full `x` index; may include self.
    fr <- dbscan::frNN(
      x = coords,
      eps = r,
      query = coords[idx, , drop = FALSE],
      sort = TRUE
    )
    from_chunks <- list()
    to_chunks <- list()
    dist_chunks <- list()
    ci <- 0L
    for (j in seq_along(idx)) {
      oi <- idx[[j]]
      nbr <- fr$id[[j]]
      d <- fr$dist[[j]]
      if (!length(nbr)) next
      keep <- nbr != oi
      nbr <- nbr[keep]
      d <- d[keep]
      keep2 <- nbr > oi
      if (!any(keep2)) next
      ci <- ci + 1L
      from_chunks[[ci]] <- rep.int(oi, sum(keep2))
      to_chunks[[ci]] <- nbr[keep2]
      dist_chunks[[ci]] <- d[keep2]
    }
    if (ci == 0L) return(empty)
    data.table::data.table(
      from = ids[unlist(from_chunks, use.names = FALSE)],
      to = ids[unlist(to_chunks, use.names = FALSE)],
      dist = unlist(dist_chunks, use.names = FALSE)
    )
  }

  parts <- .sphinx_lapply(chunks, worker, n_cores = n_cores)
  ed <- data.table::rbindlist(parts, use.names = TRUE)
  if (!nrow(ed)) return(empty)
  ed
}

#' Resolve parallel worker count for network builders
#'
#' Order: explicit `n_cores` > env `SPHINX_N_CORES` > min(detected cores, 8).
#'
#' @param n_cores optional integer
#' @return integer >= 1
#' @keywords internal
.sphinx_default_cores <- function(n_cores = NULL) {
  if (!is.null(n_cores) && length(n_cores) >= 1L) {
    nc <- suppressWarnings(as.integer(n_cores[[1L]]))
    if (is.finite(nc) && nc >= 1L) return(nc)
  }
  env <- suppressWarnings(as.integer(Sys.getenv("SPHINX_N_CORES", unset = "")))
  if (length(env) == 1L && is.finite(env) && env >= 1L) return(env)
  detected <- suppressWarnings(parallel::detectCores(logical = TRUE))
  if (!is.finite(detected) || detected < 1L) return(1L)
  max(1L, min(as.integer(detected), 8L))
}

#' lapply with optional fork parallelism (Unix); serial on Windows / n_cores=1
#' @keywords internal
.sphinx_lapply <- function(X, FUN, n_cores = 1L, ...) {
  n_cores <- max(1L, as.integer(n_cores))
  if (length(X) == 0L) return(list())
  if (length(X) == 1L || n_cores <= 1L || .Platform$OS.type == "windows") {
    return(lapply(X, FUN, ...))
  }
  parallel::mclapply(X, FUN, ..., mc.cores = min(n_cores, length(X)))
}

#' Build edges using Delaunay triangulation
#'
#' For large point sets, `deldir` (Fortran) fails with
#' "long vectors are not supported in .Fortran". In that regime we tile the
#' FOV into overlapping spatial chunks, triangulate each chunk, and merge
#' undirected edges (deduplicated).
#'
#' @param df data.table with spatial data
#' @param max_edge_length optional maximum edge length filter
#' @param k neighbor count used only if `deldir` is missing (kNN fallback)
#' @param chunk_max_points max points per chunk before tiling (default 2.5e5)
#' @return data.table with edges
#' @keywords internal
.edges_delaunay <- function(df, max_edge_length = NULL, k = 6L,
                            chunk_max_points = 250000L) {
  empty <- data.table::data.table(
    from = character(), to = character(), dist = numeric()
  )
  n <- nrow(df)
  if (n < 2L) return(empty)

  chunk_max_points <- as.integer(chunk_max_points)
  if (!is.finite(chunk_max_points) || chunk_max_points < 1000L) {
    chunk_max_points <- 250000L
  }

  # Prefer Qhull (geometry::delaunayn) for large FOVs - deldir's Fortran
  # backend hits "long vectors are not supported" around ~1e6 points.
  if (n > chunk_max_points && requireNamespace("geometry", quietly = TRUE)) {
    ed <- tryCatch(
      .edges_delaunay_qhull(df, max_edge_length = max_edge_length),
      error = function(e) {
        message(
          "geometry::delaunayn failed (", conditionMessage(e),
          "); falling back to chunked deldir."
        )
        NULL
      }
    )
    if (!is.null(ed)) return(ed)
  }

  if (!requireNamespace("deldir", quietly = TRUE)) {
    kk <- if (!is.null(k) && is.finite(k) && as.integer(k) >= 1L) as.integer(k) else 6L
    warning("Package 'deldir' not installed; falling back to kNN (k=", kk, ")")
    return(.edges_knn(df, k = kk, max_edge_length = max_edge_length))
  }

  if (n <= chunk_max_points) {
    ed <- tryCatch(
      .edges_delaunay_deldir_block(df, max_edge_length = max_edge_length),
      error = function(e) {
        if (!grepl("long vectors|Fortran", conditionMessage(e), ignore.case = TRUE)) {
          stop(e)
        }
        message(
          "deldir failed on n=", n, " (", conditionMessage(e),
          "); retrying with spatial chunks."
        )
        .edges_delaunay_chunked(
          df,
          max_edge_length = max_edge_length,
          chunk_max_points = min(chunk_max_points, max(5000L, n %/% 4L))
        )
      }
    )
    return(ed)
  }

  message(
    "Delaunay: n=", n, " - using overlapping spatial deldir chunks",
    " (chunk_max_points=", chunk_max_points, ")"
  )
  .edges_delaunay_chunked(
    df,
    max_edge_length = max_edge_length,
    chunk_max_points = chunk_max_points
  )
}

#' Qhull Delaunay via geometry::delaunayn (handles ~1e6 points)
#' @keywords internal
.edges_delaunay_qhull <- function(df, max_edge_length = NULL) {
  empty <- data.table::data.table(
    from = character(), to = character(), dist = numeric()
  )
  n <- nrow(df)
  if (n < 2L) return(empty)
  message("Delaunay: n=", n, " via geometry::delaunayn (Qhull)")
  xy <- cbind(as.numeric(df$X), as.numeric(df$Y))
  # QJ: joggle inputs to avoid precision issues on near-cocircular points
  tri <- geometry::delaunayn(xy, options = "QJ")
  if (is.null(tri) || !nrow(tri)) return(empty)
  # each row = triangle vertex indices (1-based)
  e1 <- c(tri[, 1], tri[, 2], tri[, 3])
  e2 <- c(tri[, 2], tri[, 3], tri[, 1])
  lo <- pmin(e1, e2)
  hi <- pmax(e1, e2)
  keep <- lo != hi
  lo <- lo[keep]
  hi <- hi[keep]
  key <- paste(lo, hi, sep = "\r")
  uniq <- !duplicated(key)
  lo <- lo[uniq]
  hi <- hi[uniq]
  d <- sqrt((xy[lo, 1] - xy[hi, 1])^2 + (xy[lo, 2] - xy[hi, 2])^2)
  ed <- data.table::data.table(
    from = as.character(df$Cell_ID[lo]),
    to = as.character(df$Cell_ID[hi]),
    dist = as.numeric(d)
  )
  if (!is.null(max_edge_length) && is.finite(max_edge_length)) {
    ed <- ed[dist <= max_edge_length]
  }
  message(sprintf("  Qhull Delaunay edges: %d", nrow(ed)))
  ed
}

#' Single-block deldir triangulation -> undirected edges
#' @keywords internal
.edges_delaunay_deldir_block <- function(df, max_edge_length = NULL) {
  empty <- data.table::data.table(
    from = character(), to = character(), dist = numeric()
  )
  if (nrow(df) < 2L) return(empty)
  dd <- deldir::deldir(df$X, df$Y, suppressMsge = TRUE)
  # delsgs$ind1/ind2 are 1-based point indices (faster than triang.list)
  i1 <- as.integer(dd$delsgs$ind1)
  i2 <- as.integer(dd$delsgs$ind2)
  if (!length(i1)) return(empty)
  lo <- pmin(i1, i2)
  hi <- pmax(i1, i2)
  keep <- lo != hi
  lo <- lo[keep]
  hi <- hi[keep]
  if (!length(lo)) return(empty)
  # unique undirected
  key <- paste(lo, hi, sep = "\r")
  uniq <- !duplicated(key)
  lo <- lo[uniq]
  hi <- hi[uniq]
  d <- sqrt((df$X[lo] - df$X[hi])^2 + (df$Y[lo] - df$Y[hi])^2)
  ed <- data.table::data.table(
    from = as.character(df$Cell_ID[lo]),
    to = as.character(df$Cell_ID[hi]),
    dist = as.numeric(d)
  )
  if (!is.null(max_edge_length) && is.finite(max_edge_length)) {
    ed <- ed[dist <= max_edge_length]
  }
  ed
}

#' Overlapping spatial-chunk Delaunay for large FOVs
#'
#' Each grid core is triangulated with a spatial buffer. An edge is kept if its
#' midpoint falls in that core (half-open), so every edge has a unique owner
#' tile while the buffer supplies the neighbours needed for a correct local DT.
#' @keywords internal
.edges_delaunay_chunked <- function(df, max_edge_length = NULL,
                                    chunk_max_points = 250000L) {
  empty <- data.table::data.table(
    from = character(), to = character(), dist = numeric()
  )
  DT <- data.table::as.data.table(df)[, .(Cell_ID = as.character(Cell_ID), X, Y)]
  n <- nrow(DT)
  if (n < 2L) return(empty)

  nn1 <- tryCatch(
    .nn_dist_at_rank(as.matrix(DT[, .(X, Y)]), k_nn = 1L),
    error = function(e) NA_real_
  )
  med_nn <- if (length(nn1) && all(is.finite(nn1))) stats::median(nn1) else NA_real_
  xr <- range(DT$X, na.rm = TRUE)
  yr <- range(DT$Y, na.rm = TRUE)
  dx <- max(xr[2] - xr[1], .Machine$double.eps)
  dy <- max(yr[2] - yr[1], .Machine$double.eps)
  diag_len <- sqrt(dx^2 + dy^2)
  buf_candidates <- c(
    if (is.finite(med_nn) && med_nn > 0) 100 * med_nn else NA_real_,
    if (!is.null(max_edge_length) && is.finite(max_edge_length)) {
      as.numeric(max_edge_length)
    } else {
      NA_real_
    },
    0.05 * diag_len
  )
  buffer <- max(buf_candidates[is.finite(buf_candidates)], na.rm = TRUE)
  if (!is.finite(buffer) || buffer <= 0) buffer <- 0.05 * max(diag_len, 1)

  n_tile <- max(2L, as.integer(ceiling(sqrt(n / as.numeric(chunk_max_points)))))
  x_breaks <- seq(xr[1], xr[2], length.out = n_tile + 1L)
  y_breaks <- seq(yr[1], yr[2], length.out = n_tile + 1L)
  if (length(unique(x_breaks)) < 2L) x_breaks <- c(xr[1] - 1, xr[2] + 1)
  if (length(unique(y_breaks)) < 2L) y_breaks <- c(yr[1] - 1, yr[2] + 1)

  message(sprintf(
    "  Delaunay chunks: %dx%d tiles | buffer=%.3f | chunk_max=%d",
    length(x_breaks) - 1L, length(y_breaks) - 1L, buffer, chunk_max_points
  ))

  edge_parts <- list()
  part_i <- 0L
  nx <- length(x_breaks) - 1L
  ny <- length(y_breaks) - 1L
  for (ix in seq_len(nx)) {
    for (iy in seq_len(ny)) {
      x0 <- x_breaks[ix]
      x1 <- x_breaks[ix + 1L]
      y0 <- y_breaks[iy]
      y1 <- y_breaks[iy + 1L]
      x_right <- ix == nx
      y_top <- iy == ny

      block <- DT[
        X >= (x0 - buffer) & X <= (x1 + buffer) &
          Y >= (y0 - buffer) & Y <= (y1 + buffer)
      ]
      if (nrow(block) < 2L) next

      # Ownership region for edge midpoints (half-open except final row/col)
      in_core_mid <- function(mx, my) {
        okx <- if (x_right) mx >= x0 & mx <= x1 else mx >= x0 & mx < x1
        oky <- if (y_top) my >= y0 & my <= y1 else my >= y0 & my < y1
        okx & oky
      }

      block_ed <- .edges_delaunay_block_maybe_split(
        block,
        in_core_mid = in_core_mid,
        max_edge_length = max_edge_length,
        chunk_max_points = chunk_max_points,
        depth = 0L
      )
      if (nrow(block_ed)) {
        part_i <- part_i + 1L
        edge_parts[[part_i]] <- block_ed
      }
    }
  }

  if (!length(edge_parts)) return(empty)
  ed <- data.table::rbindlist(edge_parts, use.names = TRUE, fill = TRUE)
  ed[, `:=`(a = pmin(from, to), b = pmax(from, to))]
  data.table::setorder(ed, a, b, dist)
  ed <- ed[, .SD[1L], by = .(a, b)]
  ed[, `:=`(from = a, to = b, a = NULL, b = NULL)]
  if (!is.null(max_edge_length) && is.finite(max_edge_length)) {
    ed <- ed[dist <= max_edge_length]
  }
  message(sprintf("  Delaunay chunked merge: %d unique edges", nrow(ed)))
  ed[]
}

#' Triangulate one expanded block; keep edges whose midpoint is in core
#' @keywords internal
.edges_delaunay_block_maybe_split <- function(block, in_core_mid,
                                              max_edge_length = NULL,
                                              chunk_max_points = 250000L,
                                              depth = 0L) {
  empty <- data.table::data.table(
    from = character(), to = character(), dist = numeric()
  )
  if (nrow(block) < 2L || !is.function(in_core_mid)) return(empty)

  run_one <- function(blk) {
    ed <- .edges_delaunay_deldir_block(blk, max_edge_length = max_edge_length)
    if (!nrow(ed)) return(empty)
    # map endpoints to coordinates in blk
    i_from <- match(ed$from, as.character(blk$Cell_ID))
    i_to <- match(ed$to, as.character(blk$Cell_ID))
    ok <- is.finite(i_from) & is.finite(i_to)
    ed <- ed[ok]
    i_from <- i_from[ok]
    i_to <- i_to[ok]
    mx <- (blk$X[i_from] + blk$X[i_to]) / 2
    my <- (blk$Y[i_from] + blk$Y[i_to]) / 2
    ed[in_core_mid(mx, my)]
  }

  split_quad <- function() {
    if (depth >= 12L) {
      stop("Delaunay chunk split exceeded max depth at n=", nrow(block))
    }
    xm <- stats::median(block$X)
    ym <- stats::median(block$Y)
    # Sub-cores inherit parent ownership predicate (midpoint still tested
    # against the original tile core via in_core_mid).
    quads <- list(
      block[X <= xm & Y <= ym],
      block[X <= xm & Y > ym],
      block[X > xm & Y <= ym],
      block[X > xm & Y > ym]
    )
    parts <- lapply(quads, function(q) {
      if (nrow(q) < 2L) return(empty)
      .edges_delaunay_block_maybe_split(
        q, in_core_mid, max_edge_length, chunk_max_points, depth + 1L
      )
    })
    data.table::rbindlist(parts, use.names = TRUE, fill = TRUE)
  }

  if (nrow(block) > chunk_max_points) {
    return(split_quad())
  }

  tryCatch(
    run_one(block),
    error = function(e) {
      if (!grepl("long vectors|Fortran", conditionMessage(e), ignore.case = TRUE)) {
        stop(e)
      }
      message(
        "  deldir Fortran limit at n=", nrow(block),
        " depth=", depth, " - splitting"
      )
      split_quad()
    }
  )
}

#' Build edges with a Stereopy-style sliding window + local kNN
#'
#' Windows of side `tile` (`window_size`) slide over the field with stride
#' `sliding_step` (default half the window). Inside each window a kNN graph
#' is built; overlapping windows are merged (undirected, shortest edge kept).
#'
#' @param df data.table with spatial data
#' @param k number of neighbors inside each window
#' @param tile window side length (Stereopy `d`)
#' @param sliding_step window stride (Stereopy `s`; default `tile / 2`)
#' @param max_edge_length optional maximum edge length filter
#' @return data.table with edges
#' @keywords internal
.edges_window <- function(df, k, tile, sliding_step = NULL,
                          max_edge_length = NULL) {
  empty <- data.table::data.table(
    from = character(), to = character(), dist = numeric()
  )
  DT <- df[, .(Cell_ID, X, Y)]
  n <- nrow(DT)
  if (n < 2L) return(empty)

  k <- as.integer(k)
  if (!is.finite(k) || k < 1L) {
    stop("n_neighbors must be a positive integer for window.")
  }
  d <- as.numeric(tile)[1L]
  if (!is.finite(d) || d <= 0) {
    stop("window_size must be a positive number.")
  }
  s <- if (is.null(sliding_step) || !is.finite(sliding_step) || sliding_step <= 0) {
    d / 2
  } else {
    as.numeric(sliding_step)[1L]
  }
  if (s > d) {
    warning(
      "window_slide_step > window_size; clamping step to window_size ",
      "so sliding windows still cover the field."
    )
    s <- d
  }

  xmin <- min(DT$X, na.rm = TRUE)
  xmax <- max(DT$X, na.rm = TRUE)
  ymin <- min(DT$Y, na.rm = TRUE)
  ymax <- max(DT$Y, na.rm = TRUE)

  window_starts <- function(lo, hi, d, s) {
    span <- hi - lo
    if (!is.finite(span) || span <= 0 || span <= d) return(lo)
    starts <- seq(lo, hi - d, by = s)
    last <- hi - d
    if (abs(starts[length(starts)] - last) > 1e-8 * max(1, d)) {
      starts <- c(starts, last)
    }
    unique(starts)
  }
  x_starts <- window_starts(xmin, xmax, d, s)
  y_starts <- window_starts(ymin, ymax, d, s)

  res <- vector("list", length(x_starts) * length(y_starts))
  idx <- 0L
  for (wx in x_starts) {
    xhi <- wx + d
    for (wy in y_starts) {
      yhi <- wy + d
      sub <- DT[X >= wx & X <= xhi & Y >= wy & Y <= yhi]
      if (nrow(sub) < 2L) next
      coords <- as.matrix(sub[, .(X, Y)])
      kk <- min(k + 1L, nrow(sub))
      if (kk < 2L) next
      nn <- RANN::nn2(coords, coords, k = kk)
      from_idx <- rep.int(seq_len(nrow(sub)), kk - 1L)
      to_idx <- as.vector(t(nn$nn.idx[, -1L, drop = FALSE]))
      distv <- as.vector(t(nn$nn.dists[, -1L, drop = FALSE]))
      ed <- data.table::data.table(
        from = sub$Cell_ID[from_idx],
        to = sub$Cell_ID[to_idx],
        dist = distv
      )
      if (!is.null(max_edge_length) && is.finite(max_edge_length)) {
        ed <- ed[dist <= max_edge_length]
      }
      if (nrow(ed) > 0L) {
        idx <- idx + 1L
        res[[idx]] <- ed
      }
    }
  }
  if (idx < 1L) return(empty)
  ed <- data.table::rbindlist(res[seq_len(idx)], use.names = TRUE, fill = TRUE)
  ed[, key := ifelse(from < to, paste(from, to), paste(to, from))]
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
#' @examples
#' df <- prepare_data(Sphinx:::.sphinx_example_df(40))
#' edges <- build_spatial_network(df, method = "knn", n_neighbors = 5, verbose = FALSE)
#' feat <- calculate_neighborhood_features(df, edges)
#' names(feat)
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
#' @param k Number of clusters (for kmeans). This is **not** spatial neighbor count;
#'   use `n_neighbors` in `build_spatial_network()`.
#' @param use_pca Whether to use PCA for dimensionality reduction
#' @param var_threshold Variance threshold for PCA components
#' @param n_components Explicit number of PCA components (overrides var_threshold)
#' @param min_cluster_size Minimum points per cluster (for hdbscan)
#' @param cluster_colname Name for the output cluster column
#' @return Data.table with cluster assignments in cluster_colname
#' @examples
#' \donttest{
#' df <- prepare_data(Sphinx:::.sphinx_example_df(40))
#' edges <- build_spatial_network(df, method = "knn", n_neighbors = 5, verbose = FALSE)
#' feat <- calculate_neighborhood_features(df, edges)
#' cl <- cluster_neighborhoods(feat, edges, method = "kmeans", k = 3)
#' "Neighborhood_Cluster" %in% names(cl)
#' }
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

#' Purity from an undirected edge list (row index in `ids`)
#' @keywords internal
.purity_from_edges <- function(ids, celltype, edges, min_cells) {
  n <- length(ids)
  purity <- rep(NA_real_, n)
  if (is.null(edges) || !nrow(edges)) return(purity)
  id_to_i <- setNames(seq_len(n), as.character(ids))
  types <- as.character(celltype)
  ed <- data.table::data.table(
    a = c(as.character(edges$from), as.character(edges$to)),
    b = c(as.character(edges$to), as.character(edges$from))
  )
  nbr <- ed[, .(nbrs = list(unique(b))), by = a]
  for (i in seq_len(nrow(nbr))) {
    ii <- id_to_i[[nbr$a[i]]]
    if (is.null(ii) || is.na(ii)) next
    nbs_i <- unname(id_to_i[nbr$nbrs[[i]]])
    nbs_i <- nbs_i[!is.na(nbs_i) & nbs_i != ii]
    if (length(nbs_i) < min_cells) next
    tab <- table(types[nbs_i])
    purity[ii] <- max(tab) / sum(tab)
  }
  purity
}

#' Calculate neighborhood purity using flexible neighbor definitions
#'
#' Neighbor definitions match `build_spatial_network()`: knn, radius ball,
#' Delaunay triangulation, or Stereopy-style sliding window + local kNN.
#'
#' @param df A data.frame or data.table containing spatial coordinates and cell type labels.
#' @param x_col Character, name of the X-coordinate column.
#' @param y_col Character, name of the Y-coordinate column.
#' @param celltype_col Character, name of the cell type column.
#' @param method Character, neighbor definition: "window", "radius", "knn", or "delaunay".
#' @param n_neighbors Integer, neighbor count for `knn` and `window`.
#' @param k Integer, **deprecated** neighbor count. Use `n_neighbors`.
#' @param radius Numeric search radius for `method = "radius"`. `NULL` uses
#'   the same 1-NN x 3 default as `build_spatial_network()`.
#' @param window_size Sliding-window side length for `method = "window"`.
#'   `NULL` uses the same data default as `build_spatial_network()`.
#' @param window_slide_step Window stride; `NULL` uses `window_size / 2`.
#' @param min_cells Integer, minimum number of neighbors required to compute purity.
#' @param verbose Logical, print progress messages.
#' @return A data.table with an added column `Neighborhood_Purity`.
#' @examples
#' df <- prepare_data(Sphinx:::.sphinx_example_df(40))
#' out <- calculate_neighborhood_purity(df, method = "knn", n_neighbors = 5, verbose = FALSE)
#' summary(out$Neighborhood_Purity)
#' @export
calculate_neighborhood_purity <- function(df,
                                          x_col = "X",
                                          y_col = "Y",
                                          celltype_col = "celltype",
                                          method = c("window", "radius", "knn", "delaunay"),
                                          n_neighbors = 10L,
                                          k = NULL,
                                          radius = NULL,
                                          window_size = NULL,
                                          window_slide_step = NULL,
                                          min_cells = 5,
                                          verbose = TRUE) {
  if (!is.null(k)) {
    warning(
      "`k` in calculate_neighborhood_purity() is deprecated for neighbor count. ",
      "Use n_neighbors.",
      call. = FALSE
    )
    if (missing(n_neighbors)) n_neighbors <- as.integer(k)
  }
  k <- as.integer(n_neighbors)
  if (!is.finite(k) || k < 1L) k <- 20L

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
  x <- as.numeric(df[[x_col]])
  y <- as.numeric(df[[y_col]])
  celltype <- df[[celltype_col]]
  purity <- rep(NA_real_, length(x))
  ids <- if ("Cell_ID" %in% names(df)) as.character(df$Cell_ID) else as.character(seq_along(x))
  work <- data.table::data.table(Cell_ID = ids, X = x, Y = y)

  ## ---- 3. Branch by method ------------------------------------------
  if (method == "window") {
    if (is.null(window_size) || !is.finite(window_size) || window_size <= 0) {
      opt_win <- .calc_adaptive_window_size(x, y)
      window_size <- opt_win$window_size
      if (is.null(window_slide_step)) window_slide_step <- opt_win$sliding_step
    } else if (is.null(window_slide_step)) {
      window_slide_step <- as.numeric(window_size) / 2
    }
    if (verbose) {
      message(
        "window_size=", round(window_size, 3),
        " | window_slide_step=", round(window_slide_step, 3),
        " | n_neighbors=", k
      )
    }
    win_edges <- .edges_window(
      work, k = k, tile = window_size, sliding_step = window_slide_step
    )
    purity <- .purity_from_edges(ids, celltype, win_edges, min_cells)
  }

  if (method == "radius") {
    if (is.null(radius) || !is.finite(radius) || radius <= 0) {
      d1 <- RANN::nn2(cbind(x, y), k = 2)$nn.dists[, 2]
      radius <- .default_search_radius(stats::median(d1, na.rm = TRUE), multiplier = 3)
    }
    if (verbose) message("radius=", round(radius, 3))
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
#' Counts cell-type contact edges and optionally normalizes by type abundance
#' so frequent types do not dominate the interaction matrix.
#'
#' Normalization:
#' \itemize{
#'   \item `abundance_normalized`: \eqn{c_{ij} / (n_i n_j)}
#'   \item `enrichment_matrix`: observed / expected under random labeling of
#'     the same graph, where
#'     \eqn{E_{ij} = E \cdot 2 p_i p_j} (i != j) and
#'     \eqn{E_{ii} = E \cdot p_i^2}, with \eqn{p_i = n_i / N}.
#' }
#'
#' @param df Spatial data with cell types
#' @param edges Spatial network edges (data.frame with "from" and "to")
#' @param celltype_col Column name for cell types in metadata.
#' @param normalize_abundance logical; if TRUE (default), also compute
#'   abundance-normalized and enrichment matrices.
#' @return List with `interaction_matrix` (raw counts), `network`, and when
#'   `normalize_abundance = TRUE` also `abundance_normalized`,
#'   `enrichment_matrix`, and `type_abundance`.
#' @examples
#' df <- prepare_data(Sphinx:::.sphinx_example_df(40))
#' edges <- build_spatial_network(df, method = "knn", n_neighbors = 5, verbose = FALSE)
#' intx <- analyze_spatial_interactions(df, edges)
#' dim(intx$interaction_matrix)
#' @export
analyze_spatial_interactions <- function(df, edges,
                                         celltype_col = "celltype",
                                         normalize_abundance = TRUE) {
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

  # Build interaction frequency table (each undirected edge counted once)
  inter_table <- table(pmin(valid_edges$type1, valid_edges$type2),
                       pmax(valid_edges$type1, valid_edges$type2))

  # Convert to symmetric matrix of raw contact counts
  all_types <- sort(unique(as.character(df[[celltype_col]])))
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

  # Build graph from raw counts
  g <- igraph::graph_from_adjacency_matrix(int_mat, mode = "undirected", weighted = TRUE)

  out <- list(interaction_matrix = int_mat, network = g)

  if (isTRUE(normalize_abundance)) {
    type_n <- as.numeric(table(factor(df[[celltype_col]], levels = all_types)))
    names(type_n) <- all_types
    N <- sum(type_n)
    p <- type_n / N
    n_edges <- nrow(valid_edges)

    # count / (n_i * n_j): down-weights abundant pairs
    abund_norm <- int_mat
    for (a in all_types) {
      for (b in all_types) {
        denom <- type_n[a] * type_n[b]
        abund_norm[a, b] <- if (denom > 0) int_mat[a, b] / denom else NA_real_
      }
    }

    # Observed/expected under random type labels on the same edge set
    expected <- outer(p, p) * n_edges
    # off-diagonal: two orderings of an undirected edge
    expected <- expected * 2
    diag(expected) <- (p^2) * n_edges
    enrichment <- int_mat / expected
    enrichment[!is.finite(enrichment)] <- NA_real_

    out$abundance_normalized <- abund_norm
    out$enrichment_matrix <- enrichment
    out$type_abundance <- type_n
  }

  out
}


#' Generate color palette for visualizations
#'
#' Soft candy colors that stay distinct (no neon / fluorescent hues).
#' Continuous heatmaps keep the pink sequential ramp unchanged.
#'
#' @param n Number of colors needed
#' @return Vector of color codes
#' @examples
#' get_color_palette(5)
#' @export
get_color_palette <- function(n) {
  # One color per hue family; spread around the wheel for richness + transparency robustness
  candy_base <- c(
    "#E05C6E",  # rose-red
    "#4EA8DE",  # sky blue
    "#E8C04A",  # gold
    "#4CB87A",  # green
    "#8B6BC9",  # purple
    "#E8884A",  # orange
    "#3DB8A0",  # teal
    "#C45BA0",  # magenta
    "#A67C52",  # brown
    "#5B7FD6",  # indigo
    "#A8C75A",  # yellow-green
    "#6B7C85",  # slate
    "#D4A017",  # amber
    "#2E8B57",  # forest
    "#E76F51",  # terracotta
    "#45A8D0",  # cyan
    "#9B59B6",  # orchid
    "#1ABC9C",  # mint
    "#3498DB",  # clear blue
    "#F39C12"   # sunflower
  )
  if (n <= 0) {
    return(character(0))
  } else if (n <= length(candy_base)) {
    return(candy_base[seq_len(n)])
  } else {
    extra_n <- n - length(candy_base)
    # Additional hues offset so they do not land on the same families as candy_base
    extra <- scales::hue_pal(h = c(8, 368), c = 70, l = 52)(length(candy_base) + extra_n)
    extra <- setdiff(toupper(extra), toupper(candy_base))
    if (length(extra) < extra_n) {
      extra <- c(extra, scales::hue_pal(h = c(20, 380), c = 55, l = 45)(extra_n))
    }
    return(c(candy_base, extra[seq_len(extra_n)]))
  }
}

#' Stable cell-type -> color mapping (same label => same color across plots)
#'
#' By default colors are assigned from a fixed soft pastel pool via a
#' deterministic string hash, so a cell type keeps its color even when other
#' types are absent. If `all_levels` is supplied, colors are assigned by
#' alphabetical order over that full set (preferred when the complete type
#' list is known).
#'
#' @param levels Character vector of cell-type / category labels to color
#' @param all_levels Optional full set of labels for ordered assignment
#' @return Named character vector of hex colors
#' @examples
#' assign_celltype_colors(c("B", "T", "Macrophage"))
#' @export
assign_celltype_colors <- function(levels, all_levels = NULL) {
  labs <- unique(as.character(levels))
  labs <- labs[!is.na(labs) & nzchar(labs)]
  if (!length(labs)) {
    return(setNames(character(0), character(0)))
  }

  # Alphabetical assignment over the reference set -> unique colors, stable
  # whenever the same full type list is used (pass all_levels when zooming).
  ref <- if (is.null(all_levels)) {
    sort(labs)
  } else {
    ref0 <- sort(unique(as.character(all_levels)))
    ref0[!is.na(ref0) & nzchar(ref0)]
  }
  if (!length(ref)) {
    return(setNames(character(0), character(0)))
  }
  cols <- get_color_palette(length(ref))
  names(cols) <- ref
  cols
}

#' Candy sequential palette for continuous heatmaps
#' @param n Number of colors
#' @keywords internal
.sphinx_soft_sequential <- function(n = 256) {
  grDevices::colorRampPalette(c(
    "#FFF5F8", "#FFE0EC", "#FFB3D1", "#FF8FAB",
    "#F06595", "#E63980", "#C9184A"
  ))(n)
}

#' Candy diverging palette (sky blue - white - candy pink)
#' @keywords internal
.sphinx_soft_diverging <- function(n = 256) {
  grDevices::colorRampPalette(c(
    "#4CC9F0", "#90E0EF", "#CAF0F8", "#FFF5F8",
    "#FFB3D1", "#FF8FAB", "#E63980"
  ))(n)
}

#' Publication theme for Sphinx ggplot visuals
#' @param base_size Base font size (default 14 for SCI figures)
#' @param grid Show panel grid (default FALSE)
#' @keywords internal
.sphinx_sci_theme <- function(base_size = 14, grid = FALSE) {
  base_size <- max(8, as.numeric(base_size))
  th <- ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = "black", size = base_size),
      plot.title = ggplot2::element_text(
        face = "plain", size = max(8, base_size + 4), hjust = 0.5,
        margin = ggplot2::margin(b = 10)
      ),
      plot.subtitle = ggplot2::element_text(
        face = "plain", size = max(8, base_size), hjust = 0.5
      ),
      axis.title = ggplot2::element_text(face = "plain", size = max(8, base_size + 1)),
      axis.text = ggplot2::element_text(
        face = "plain", size = max(8, base_size - 1), colour = "black"
      ),
      legend.title = ggplot2::element_text(face = "plain", size = max(8, base_size)),
      legend.text = ggplot2::element_text(face = "plain", size = max(8, base_size - 2)),
      strip.text = ggplot2::element_text(face = "plain", size = max(8, base_size)),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.8),
      axis.line = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
  if (isTRUE(grid)) {
    th <- th + ggplot2::theme(
      panel.grid.major = ggplot2::element_line(colour = "grey92", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank()
    )
  } else {
    th <- th + ggplot2::theme(panel.grid = ggplot2::element_blank())
  }
  th
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
#' @param base_size Base font size in points (default: 14; minimum 8)
#' @param color_palette Optional named color vector; if NULL uses assign_celltype_colors()
#' @param all_levels Optional full label set for stable colors when subsetting / zooming
#' @param save_path Output file path (optional)
#' @param width Plot width in inches (default: 10)
#' @param height Plot height in inches (default: 8)
#' @return ggplot object and saves plot to file if save_path provided
#' @examples
#' \donttest{
#' df <- prepare_data(Sphinx:::.sphinx_example_df(40))
#' p <- visualize_spatial_distribution(df, save_path = tempfile(fileext = ".pdf"))
#' class(p)
#' }
#' @export
visualize_spatial_distribution <- function(df,
                                           x_col = "X",
                                           y_col = "Y",
                                           celltype_col = "celltype",
                                           point_size = 1.5,
                                           point_alpha = 0.85,
                                           point_shape = 16,
                                           legend_point_size = 3,
                                           title = NULL,
                                           legend.position = "right",
                                           base_size = 14,
                                           color_palette = NULL,
                                           all_levels = NULL,
                                           save_path = NULL,
                                           width = 10,
                                           height = 8) {

  set.seed(123)

  # Create local copies of columns
  df$plot_X <- df[[x_col]]
  df$plot_Y <- df[[y_col]]
  df$plot_celltype <- df[[celltype_col]]

  # Stable cell-type colors across plots
  if (is.null(color_palette)) {
    palette <- assign_celltype_colors(df$plot_celltype, all_levels = all_levels)
  } else {
    palette <- color_palette
  }

  # Create spatial plot (no internal grid)
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
    ggplot2::labs(
      title = if (!is.null(title)) title else "Spatial Distribution of Cell Types",
      x = "X Coordinate",
      y = "Y Coordinate",
      color = "Cell Type"
    ) +
    ggplot2::coord_fixed() +
    .sphinx_sci_theme(base_size = base_size, grid = FALSE) +
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
#' @param base_size Base font size in points (default: 14; minimum 8)
#' @param show_values Whether to print distance values in cells (default: TRUE)
#' @return ggplot object and saves plot to file if save_path provided
#' @examples
#' \donttest{
#' df <- prepare_data(Sphinx:::.sphinx_example_df(40))
#' dist_res <- calculate_celltype_distances(df, celltype_col = "celltype")
#' p <- visualize_distance_heatmap(dist_res, save_path = tempfile(fileext = ".pdf"))
#' class(p)
#' }
#' @export
visualize_distance_heatmap <- function(dist_result,
                                       save_path = NULL,
                                       width = 12,
                                       height = 10,
                                       base_size = 14,
                                       show_values = TRUE) {

  # Prepare data
  dist_matrix <- dist_result$distance_matrix
  dist_df <- reshape2::melt(dist_matrix, varnames = c("Source", "Target"), value.name = "Distance")
  # Self-pair / missing distances: keep tiles blank, skip labels
  dist_df$Distance_plot <- dist_df$Distance
  dist_ok <- dist_df[!is.na(dist_df$Distance), , drop = FALSE]

  # Calculate optimal text color against candy sequential fill
  mid_cut <- stats::quantile(dist_ok$Distance, 0.55, na.rm = TRUE)
  dist_ok$TextColor <- ifelse(dist_ok$Distance < mid_cut, "#4A1942", "#FFF5F8")

  soft_cols <- c(
    "#FFF5F8", "#FFE0EC", "#FFB3D1", "#FF8FAB",
    "#F06595", "#E63980", "#C9184A"
  )

  # Text size scales with figure size; floor at 8 pt (ggplot size ~= mm)
  n_lab <- max(nrow(dist_matrix), ncol(dist_matrix), 1L)
  inch_per_cell <- min(as.numeric(width), as.numeric(height)) / n_lab
  pt_to_gg <- 72.27 / 25.4  # ggplot2::.pt
  text_size_min <- 8 / pt_to_gg
  text_size <- max(text_size_min, min(5.5, inch_per_cell * 3.2))

  # Create heatmap
  p <- ggplot2::ggplot(dist_df, ggplot2::aes(x = Source, y = Target, fill = Distance)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.4) +
    ggplot2::scale_fill_gradientn(
      name = "Mean Distance",
      colours = soft_cols,
      na.value = "grey95",
      values = scales::rescale(c(
        min(dist_ok$Distance, na.rm = TRUE),
        stats::quantile(dist_ok$Distance, 0.2, na.rm = TRUE),
        stats::quantile(dist_ok$Distance, 0.4, na.rm = TRUE),
        stats::quantile(dist_ok$Distance, 0.55, na.rm = TRUE),
        stats::quantile(dist_ok$Distance, 0.7, na.rm = TRUE),
        stats::quantile(dist_ok$Distance, 0.85, na.rm = TRUE),
        max(dist_ok$Distance, na.rm = TRUE)
      ))
    ) +
    ggplot2::labs(
      title = "Mean Distance Between Cell Types",
      x = "Source Cell Type",
      y = "Target Cell Type"
    ) +
    .sphinx_sci_theme(base_size = base_size, grid = FALSE) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "right"
    )

  if (isTRUE(show_values)) {
    p <- p +
      ggplot2::geom_text(
        data = dist_ok,
        ggplot2::aes(label = round(Distance, 1), color = TextColor),
        size = text_size,
        fontface = "plain",
        show.legend = FALSE
      ) +
      ggplot2::scale_color_identity()
  }

  p <- p +
    ggplot2::geom_tile(
      data = dist_df[as.character(dist_df$Source) == as.character(dist_df$Target), ],
      ggplot2::aes(x = Source, y = Target),
      color = "grey30",
      linewidth = 0.8,
      fill = NA
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
#' @param base_size Base font size in points (default: 14; minimum 8)
#' @return ggplot object and saves plot to file if save_path provided
#' @examples
#' \donttest{
#' df <- prepare_data(Sphinx:::.sphinx_example_df(40))
#' dist_res <- calculate_celltype_distances(df, celltype_col = "celltype")
#' p <- visualize_distance_parallel(dist_res, save_path = tempfile(fileext = ".pdf"))
#' class(p)
#' }
#' @export
visualize_distance_parallel <- function(dist_result,
                                        save_path = NULL,
                                        width = 14,
                                        height = 8,
                                        base_size = 14) {

  # Prepare data
  dist_matrix <- dist_result$distance_matrix
  cell_types <- rownames(dist_matrix)

  # Convert to long format
  plot_data <- data.frame(dist_matrix)
  plot_data$Source <- cell_types
  plot_data <- reshape2::melt(plot_data, id.vars = "Source",
                              variable.name = "Target",
                              value.name = "Distance")

  # Soft candy colors, stable by cell-type name
  colors <- assign_celltype_colors(cell_types)

  # Create parallel coordinates plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Target, y = Distance, group = Source, color = Source)) +
    ggplot2::geom_line(linewidth = 1.1, alpha = 0.85) +
    ggplot2::geom_point(size = 2.8, alpha = 0.9) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(
      title = "Cell Type Distance Relationships",
      x = "Target Cell Type",
      y = "Distance to Target",
      color = "Source Cell Type"
    ) +
    .sphinx_sci_theme(base_size = base_size, grid = TRUE) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid.major.x = ggplot2::element_blank(),
      legend.position = "right"
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(linewidth = 1.2, size = 3)))

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
#' @param base_size Base font size in points (default: 14; minimum 8)
#' @param save_path Output file path (optional)
#' @param width Plot width in inches (default: 10)
#' @param height Plot height in inches (default: 8)
#' @return ggplot object and saves plot to file if save_path provided
#' @examples
#' \donttest{
#' df <- prepare_data(Sphinx:::.sphinx_example_df(40))
#' df <- calculate_neighborhood_purity(df, method = "knn", n_neighbors = 5, verbose = FALSE)
#' p <- visualize_neighborhood_purity(df, save_path = tempfile(fileext = ".pdf"))
#' class(p)
#' }
#' @export
visualize_neighborhood_purity <- function(df,
                                          x_col = "X",
                                          y_col = "Y",
                                          max_points = 10000,
                                          point_size = 1.5,
                                          point_alpha = 0.9,
                                          point_shape = 16,
                                          title = NULL,
                                          legend.position = "right",
                                          base_size = 14,
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

  # Candy sequential colors for purity
  soft_cols <- c("#FFF5F8", "#FFD6E7", "#FFB3D1", "#FF8FAB", "#E63980", "#9B1D5A")

  # Create purity visualization (no internal grid)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = plot_X, y = plot_Y, color = Neighborhood_Purity)) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha, shape = point_shape) +
    ggplot2::scale_color_gradientn(colours = soft_cols, name = "Purity") +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title = if (!is.null(title)) title else "Neighborhood Purity",
      x = "X Coordinate",
      y = "Y Coordinate"
    ) +
    .sphinx_sci_theme(base_size = base_size, grid = FALSE) +
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
#' @examples
#' df <- prepare_data(Sphinx:::.sphinx_example_df(40))
#' df$Neighborhood_Cluster <- sample(1:3, nrow(df), replace = TRUE)
#' comp <- calculate_cluster_composition(df)
#' head(comp)
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
#' @param cell_fontsize Font size for in-cell numbers (default: 11; minimum 8)
#' @return ComplexHeatmap object and saves plot to file if save_path provided
#' @examples
#' \donttest{
#' if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
#'   df <- prepare_data(Sphinx:::.sphinx_example_df(40))
#'   df$Neighborhood_Cluster <- sample(1:3, nrow(df), replace = TRUE)
#'   comp <- calculate_cluster_composition(df)
#'   plot_composition_heatmap(comp, save_path = tempfile(fileext = ".pdf"))
#' }
#' }
#' @export
plot_composition_heatmap <- function(composition_df,
                                     cluster_col = "Neighborhood_Cluster",
                                     celltype_col = "celltype",
                                     value_col = "proportion",
                                     save_path = NULL,
                                     width = 12,
                                     height = 10,
                                     cell_fontsize = 11) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE) ||
      !requireNamespace("circlize", quietly = TRUE)) {
    stop("Packages 'ComplexHeatmap' and 'circlize' are required for plot_composition_heatmap().")
  }
  cell_fontsize <- max(8, as.numeric(cell_fontsize))
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

  # Define color function for Z-scores (candy blue-white-pink)
  z_range <- stats::quantile(z_matrix, c(0.05, 0.95), na.rm = TRUE)
  max_abs <- max(abs(z_range))
  col_limits <- c(-max_abs, max_abs)
  col_fun <- circlize::colorRamp2(
    c(col_limits[1], 0, col_limits[2]),
    c("#4CC9F0", "#FFF5F8", "#E63980")
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
    row_names_gp = grid::gpar(fontsize = max(8, 12)),
    column_names_gp = grid::gpar(fontsize = max(8, 12)),
    row_title_gp = grid::gpar(fontsize = max(8, 14)),
    column_title_gp = grid::gpar(fontsize = max(8, 14)),
    heatmap_legend_param = list(
      title_position = "topcenter",
      title_gp = grid::gpar(fontsize = max(8, 12)),
      labels_gp = grid::gpar(fontsize = max(8, 11))
    ),
    right_annotation = row_ha,
    row_dend_width = grid::unit(1.5, "cm"),
    column_dend_height = grid::unit(1.5, "cm"),
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid::grid.text(
        sprintf("%.2f", comp_matrix[i, j]),
        x, y,
        gp = grid::gpar(fontsize = cell_fontsize, col = "grey20")
      )
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
#' @param base_size Base font size in points (default: 14; minimum 8)
#' @return ggplot object and saves plot to file if save_path provided
#' @examples
#' \donttest{
#' df <- prepare_data(Sphinx:::.sphinx_example_df(40))
#' df$Neighborhood_Cluster <- sample(1:3, nrow(df), replace = TRUE)
#' comp <- calculate_cluster_composition(df)
#' p <- plot_composition_barplot(comp, save_path = tempfile(fileext = ".pdf"))
#' class(p)
#' }
#' @export
plot_composition_barplot <- function(composition_df,
                                     cluster_col = "Neighborhood_Cluster",
                                     celltype_col = "celltype",
                                     value_col = "proportion",
                                     save_path = NULL,
                                     width = 12,
                                     height = 8,
                                     base_size = 14) {

  # Generate stable color palette
  celltypes <- unique(as.character(composition_df[[celltype_col]]))
  palette <- assign_celltype_colors(celltypes)

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
    ggplot2::geom_bar(stat = "identity", position = "stack", width = 0.8) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(),
                                expand = c(0, 0)) +
    ggplot2::labs(
      x = "Neighborhood Cluster",
      y = "Cell Type Proportion",
      fill = "Cell Type"
    ) +
    .sphinx_sci_theme(base_size = base_size, grid = FALSE) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      legend.position = "right",
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
#'   (raw counts, enrichment, or abundance-normalized scores)
#' @param transform Apply log2(x+1) transformation (default: TRUE). Set FALSE
#'   for enrichment / abundance-normalized matrices.
#' @param color_palette Color palette function (default: soft sequential)
#' @param main Heatmap title
#' @param display_numbers Whether to show numbers in cells (default: FALSE)
#' @param number_format sprintf format for cell labels (default: auto)
#' @param cellwidth Cell width in points (default: 14)
#' @param cellheight Cell height in points (default: 14)
#' @param fontsize Base font size (default: 12)
#' @param save_path Output file path (optional)
#' @param width Plot width in inches (default: 10)
#' @param height Plot height in inches (default: 8)
#' @return pheatmap object and saves plot to file if save_path provided
#' @examples
#' \donttest{
#' df <- prepare_data(Sphinx:::.sphinx_example_df(40))
#' edges <- build_spatial_network(df, method = "knn", n_neighbors = 5, verbose = FALSE)
#' intx <- analyze_spatial_interactions(df, edges)
#' visualize_interaction_heatmap(intx$interaction_matrix,
#'   save_path = tempfile(fileext = ".pdf"))
#' visualize_interaction_heatmap(intx$enrichment_matrix, transform = FALSE,
#'   main = "Abundance-normalized contact enrichment")
#' }
#' @export
visualize_interaction_heatmap <- function(interaction_matrix,
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
                                          height = 8) {

  fontsize <- max(8, as.numeric(fontsize))
  mat <- interaction_matrix

  # Optional transformation
  if (isTRUE(transform)) {
    mat <- log2(mat + 1)
  }

  if (is.null(number_format)) {
    number_format <- if (isTRUE(transform) || all(mat == floor(mat), na.rm = TRUE)) {
      "%.0f"
    } else {
      "%.2f"
    }
  }

  # Soft sequential palette (avoid high-saturation inferno)
  if (is.null(color_palette)) {
    col_fun <- .sphinx_soft_sequential(256)
  } else if (is.function(color_palette)) {
    col_fun <- color_palette(256)
  } else {
    col_fun <- color_palette
  }

  # Create heatmap with smaller cells and no cell numbers by default
  hm <- pheatmap::pheatmap(
    mat = mat,
    main = main,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    fontsize = fontsize,
    fontsize_row = fontsize,
    fontsize_col = fontsize,
    display_numbers = display_numbers,
    number_format = number_format,
    number_color = "grey20",
    fontsize_number = max(8, fontsize - 3),
    color = col_fun,
    border_color = "white",
    cellwidth = cellwidth,
    cellheight = cellheight,
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
#' Draws **real pairwise cell-cell edges** using spatial coordinates. Supports
#' local zoom via `xlim`/`ylim` or `zoom_center` + `zoom_radius`. When zooming,
#' edges with both endpoints inside the window are kept so local connectivity
#' is shown faithfully.
#'
#' @param df Spatial data
#' @param edges Spatial network edges
#' @param celltype_col Cell type column name (default: "celltype")
#' @param x_col X coordinate column name (default: "X")
#' @param y_col Y coordinate column name (default: "Y")
#' @param edge_mode Edge selection mode: "all", "top", or "random" (default: "all")
#' @param top_n When edge_mode = "top", keep this many strongest edges (default: 1000)
#' @param max_edges When edge_mode = "random", sample this many edges (default: 1000)
#' @param xlim Optional x-axis limits for local zoom, e.g. `c(xmin, xmax)`
#' @param ylim Optional y-axis limits for local zoom, e.g. `c(ymin, ymax)`
#' @param zoom_center Optional center `c(x, y)` for circular/rectangular zoom
#' @param zoom_radius Half-width of zoom window around `zoom_center` (same units as coords)
#' @param point_size Point size (default: 1.5)
#' @param point_alpha Point transparency (default: 1; keep opaque so cell-type colors stay clear over edges)
#' @param edge_size_range Edge size range (default: c(0.3, 2.0))
#' @param edge_alpha_range Edge alpha range (default: c(0.3, 0.9))
#' @param edge_color Edge color (default: "grey35")
#' @param show_points Whether to show points (default: TRUE)
#' @param legend_point_size Legend point size (default: 3)
#' @param base_size Base font size (default: 14)
#' @param title Plot title
#' @param save_path Output file path (optional)
#' @param width Plot width in inches (default: 12)
#' @param height Plot height in inches (default: 10)
#' @return ggplot object and saves plot to file if save_path provided
#' @examples
#' \donttest{
#' df <- prepare_data(Sphinx:::.sphinx_example_df(40))
#' edges <- build_spatial_network(df, method = "knn", n_neighbors = 5, verbose = FALSE)
#' p <- visualize_spatial_network(df, edges, save_path = tempfile(fileext = ".pdf"))
#' # Local zoom around tissue center:
#' # visualize_spatial_network(df, edges, zoom_center = c(500, 500), zoom_radius = 80)
#' class(p)
#' }
#' @export
visualize_spatial_network <- function(df, edges,
                                      celltype_col = "celltype",
                                      x_col = "X",
                                      y_col = "Y",
                                      edge_mode = c("all", "top", "random"),
                                      top_n = 1000,
                                      max_edges = 1000,
                                      xlim = NULL,
                                      ylim = NULL,
                                      zoom_center = NULL,
                                      zoom_radius = NULL,
                                      point_size = 1.5,
                                      point_alpha = 1,
                                      edge_size_range = c(0.3, 2.0),
                                      edge_alpha_range = c(0.15, 0.45),
                                      edge_color = "grey35",
                                      show_points = TRUE,
                                      legend_point_size = 3,
                                      base_size = 14,
                                      title = NULL,
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

  # Resolve zoom window (local magnification)
  if (!is.null(zoom_center) && !is.null(zoom_radius)) {
    if (length(zoom_center) != 2L || !is.finite(zoom_radius) || zoom_radius <= 0) {
      stop("zoom_center must be c(x, y) and zoom_radius a positive number")
    }
    xlim <- c(zoom_center[1] - zoom_radius, zoom_center[1] + zoom_radius)
    ylim <- c(zoom_center[2] - zoom_radius, zoom_center[2] + zoom_radius)
  }
  use_zoom <- !is.null(xlim) && !is.null(ylim)
  if (xor(is.null(xlim), is.null(ylim))) {
    stop("Provide both xlim and ylim (or zoom_center + zoom_radius) for local zoom")
  }

  # Prepare coordinates
  edges$from <- as.character(edges$from)
  edges$to <- as.character(edges$to)
  df$Cell_ID <- as.character(df$Cell_ID)
  # Keep full cell-type color map even after local zoom
  all_celltype_levels <- unique(as.character(df[[celltype_col]]))

  df_coords <- df[, .(Cell_ID, get(x_col), get(y_col))]
  data.table::setnames(df_coords, c("Cell_ID", "X", "Y"))

  # Subset cells to zoom window first so edges reflect true local links
  if (use_zoom) {
    keep_ids <- df_coords[X >= xlim[1] & X <= xlim[2] & Y >= ylim[1] & Y <= ylim[2], Cell_ID]
    df <- df[Cell_ID %in% keep_ids]
    df_coords <- df_coords[Cell_ID %in% keep_ids]
    message(
      "Local zoom: xlim=[", round(xlim[1], 1), ", ", round(xlim[2], 1),
      "], ylim=[", round(ylim[1], 1), ", ", round(ylim[2], 1),
      "] | cells=", nrow(df)
    )
  }

  edges_with_coords <- merge(edges, df_coords, by.x = "from", by.y = "Cell_ID")
  data.table::setnames(edges_with_coords, c("X", "Y"), c("from_x", "from_y"))

  edges_with_coords <- merge(edges_with_coords, df_coords, by.x = "to", by.y = "Cell_ID")
  data.table::setnames(edges_with_coords, c("X", "Y"), c("to_x", "to_y"))

  # Calculate edge strength from real distances
  if ("dist" %in% names(edges_with_coords)) {
    edges_with_coords[, strength := 1 / pmax(dist, .Machine$double.eps)]
  } else {
    edges_with_coords[, dist := sqrt((to_x - from_x)^2 + (to_y - from_y)^2)]
    edges_with_coords[, strength := 1 / pmax(dist, .Machine$double.eps)]
  }

  # Filter edges based on mode (after zoom so local real links are preserved)
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
    message("Show all ", nrow(edges_with_coords), " real cell-cell edges")
  }

  plot_title <- if (!is.null(title)) {
    title
  } else if (use_zoom) {
    "Spatial Network (local zoom)"
  } else {
    "Spatial Interaction Network"
  }

  # Create base plot (no internal grid)
  p <- ggplot2::ggplot() +
    ggplot2::labs(
      title = plot_title,
      x = x_col,
      y = y_col
    ) +
    .sphinx_sci_theme(base_size = base_size, grid = FALSE)

  # Draw pairwise cell-cell edges first (under points).
  # Do NOT map alpha via a scale -- that can wash out opaque cell-type points.
  if (nrow(edges_with_coords) > 0) {
    edge_alpha_fixed <- mean(edge_alpha_range)
    p <- p + ggplot2::geom_segment(
      data = edges_with_coords,
      ggplot2::aes(
        x = .data$from_x, y = .data$from_y,
        xend = .data$to_x, yend = .data$to_y,
        linewidth = .data$strength
      ),
      color = edge_color,
      alpha = edge_alpha_fixed,
      show.legend = FALSE
    ) +
      ggplot2::scale_linewidth_continuous(range = edge_size_range, guide = "none")
  }

  # Overlay cell points on top of edges (fully opaque + white halo so colors stay vivid)
  if (show_points && nrow(df) > 0) {
    palette <- assign_celltype_colors(
      df[[celltype_col]],
      all_levels = all_celltype_levels
    )
    pts <- as.data.frame(df)

    # White underlay keeps cell-type colors from looking washed out on dense edges
    p <- p +
      ggplot2::geom_point(
        data = pts,
        ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]]),
        size = point_size * 1.55,
        color = "white",
        alpha = 1,
        shape = 16,
        show.legend = FALSE
      ) +
      ggplot2::geom_point(
        data = pts,
        ggplot2::aes(
          x = .data[[x_col]],
          y = .data[[y_col]],
          color = .data[[celltype_col]]
        ),
        size = point_size,
        alpha = point_alpha,
        shape = 16,
        stroke = 0
      ) +
      ggplot2::scale_color_manual(
        name = "Cell Type",
        values = palette,
        guide = ggplot2::guide_legend(
          override.aes = list(
            size = legend_point_size,
            alpha = 1,
            shape = 16,
            linetype = 0,
            linewidth = 0,
            stroke = 0
          )
        )
      )
  }

  if (use_zoom) {
    p <- p + ggplot2::coord_fixed(xlim = xlim, ylim = ylim, expand = FALSE)
  } else {
    p <- p + ggplot2::coord_fixed()
  }

  # Save plot if path provided
  if (!is.null(save_path)) {
    ggplot2::ggsave(save_path, p, width = width, height = height, dpi = 150)
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
#' @param base_size Base font size in points (default: 14; minimum 8)
#' @return ggraph object and saves plot to file if save_path provided
#' @examples
#' \donttest{
#' df <- prepare_data(Sphinx:::.sphinx_example_df(40))
#' edges <- build_spatial_network(df, method = "knn", n_neighbors = 5, verbose = FALSE)
#' intx <- analyze_spatial_interactions(df, edges)
#' p <- visualize_interaction_network(intx$network, save_path = tempfile(fileext = ".pdf"))
#' class(p)
#' }
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
                                          layout = "fr",
                                          base_size = 14) {

  base_size <- max(8, as.numeric(base_size))
  # ggraph/ggplot text size is in mm; enforce >= 8 pt
  label_size <- max(8 / (72.27 / 25.4), as.numeric(label_size))

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

  # Simplify large networks
  if (igraph::vcount(network) > max_nodes) {
    node_importance <- igraph::strength(network)
    top_nodes <- order(node_importance, decreasing = TRUE)[1:max_nodes]
    network <- igraph::induced_subgraph(network, top_nodes)
    message("Network too large, filtered to top ", max_nodes, " important nodes")
  }

  # Soft cell-type colors (stable by name)
  palette <- assign_celltype_colors(igraph::V(network)$name)

  # Create network plot -- legend keeps cell type dots only
  plot <- ggraph::ggraph(network, layout = layout) +
    ggraph::geom_edge_link(
      ggplot2::aes(width = weight, alpha = weight),
      color = "grey60",
      show.legend = FALSE
    ) +
    ggraph::scale_edge_width_continuous(
      range = edge_size_range,
      guide = "none"
    ) +
    ggraph::scale_edge_alpha_continuous(
      range = c(0.25, 0.75),
      guide = "none"
    ) +
    ggraph::geom_node_point(
      ggplot2::aes(size = strength, color = name),
      alpha = 0.9,
      shape = 16
    ) +
    ggplot2::scale_size_continuous(
      range = node_size_range,
      guide = "none"
    ) +
    ggplot2::scale_color_manual(
      name = "Cell Type",
      values = palette,
      guide = ggplot2::guide_legend(
        override.aes = list(
          size = 4,
          shape = 16,
          alpha = 1,
          linetype = 0,
          linewidth = 0,
          stroke = 0
        )
      )
    ) +
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
    ggplot2::labs(title = "Spatial Cell Interaction Network") +
    ggplot2::theme_void(base_size = base_size) +
    ggplot2::theme(
      legend.position = "right",
      legend.title = ggplot2::element_text(face = "plain", size = max(8, base_size)),
      legend.text = ggplot2::element_text(face = "plain", size = max(8, base_size - 2)),
      plot.title = ggplot2::element_text(
        hjust = 0.5, face = "plain", size = max(8, base_size + 4)
      ),
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
#' @param celltype_col Cell type column name (default: "celltype")
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
#' @examples
#' \donttest{
#' df <- prepare_data(Sphinx:::.sphinx_example_df(40))
#' df$Neighborhood_Cluster <- sample(1:3, nrow(df), replace = TRUE)
#' p <- visualize_voronoi(df, celltype_col = "celltype",
#'   save_path = tempfile(fileext = ".pdf"))
#' class(p)
#' }
#' @export
visualize_voronoi <- function(df,
                              x_col = "X",
                              y_col = "Y",
                              coloring = c("celltype", "neighborhood"),
                              highlight_cluster = NULL,
                              celltype_col = "celltype",
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
  need_nhood <- identical(coloring, "neighborhood") || !is.null(highlight_cluster)

  # Check if required columns exist
  if (!x_col %in% names(df)) stop("Column '", x_col, "' not found in data frame")
  if (!y_col %in% names(df)) stop("Column '", y_col, "' not found in data frame")
  if (!celltype_col %in% names(df)) stop("Column '", celltype_col, "' not found in data frame")
  if (need_nhood && !neighborhood_col %in% names(df)) {
    stop("Column '", neighborhood_col, "' not found in data frame")
  }

  # Create local column copies
  df <- df %>%
    dplyr::mutate(
      plot_X = .data[[x_col]],
      plot_Y = .data[[y_col]],
      plot_celltype = as.factor(.data[[celltype_col]]),
      plot_neighborhood = if (need_nhood) as.factor(.data[[neighborhood_col]]) else NA
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

    if (is.null(celltype_palette)) {
      color_palette <- assign_celltype_colors(polygons$color_group)
    } else {
      color_palette <- celltype_palette
    }
  } else {
    polygons$color_group <- polygons$neighborhood
    legend_title <- "Neighborhood Cluster"
    color_palette <- assign_celltype_colors(polygons$color_group)
  }

  # Handle highlighting logic
  if (!is.null(highlight_cluster) && coloring == "neighborhood") {
    # Add highlight identifier column
    polygons$highlight <- ifelse(polygons$neighborhood == highlight_cluster,
                                 "Highlighted", "Background")

    # Get cell type color scheme (stable by name)
    if (is.null(celltype_palette)) {
      celltype_palette <- assign_celltype_colors(polygons$celltype)
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
        title = paste("Voronoi Diagram - Highlighting", highlight_cluster)
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
      if (is.null(celltype_palette)) {
        color_palette <- assign_celltype_colors(polygons$color_group)
      } else {
        color_palette <- celltype_palette
      }
    }

    polygons$border_color <- "white"

    p <- ggplot2::ggplot(polygons, ggplot2::aes(x = x, y = y, group = cell_id, fill = color_group)) +
      ggplot2::geom_polygon(
        ggplot2::aes(color = border_color),
        alpha = 0.8, linewidth = 0.2
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
    ggplot2::ylim(y_range[1] - y_margin, y_range[2] + y_margin) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "plain", size = max(8, 16)),
      legend.title = ggplot2::element_text(face = "plain", size = max(8, 12)),
      legend.text = ggplot2::element_text(face = "plain", size = max(8, 11))
    )

  # Save plot if path provided
  if (!is.null(save_path)) {
    ggplot2::ggsave(save_path, p, width = width, height = height)
    message("Voronoi diagram saved to: ", save_path)
  }

  return(p)
}
