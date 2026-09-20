#!/usr/bin/env Rscript
## Fig3 Spatch COAD-3
## a  Spatial cell-type map
## b  Cell-type pairwise distance heatmap
## c  Epithelial / endothelial / fibroblast distance ridges
## d  Four-method cellular-neighborhood comparison
## e  CN composition Z-score heatmap
## f  Neighborhood purity spatial map
## g  Purity comparison (epi- vs stromal-dominant)
## h  CN interaction network
## i  Spatial network (full and cell-type zoom)
## j  CellChat ligand–receptor bubble


rm(list = ls())
set.seed(1234)
options(bitmapType = "cairo")

root_dir <- local({
  ca <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", ca[grep("^--file=", ca)])
  if (length(f) && nzchar(f[[1]])) dirname(normalizePath(f[[1]])) else normalizePath(getwd())
})
.sphinx_r <- local({
  cands <- c(file.path(dirname(root_dir), "R"), file.path(root_dir, "Sphinx", "R"))
  hit <- cands[dir.exists(cands)]
  if (!length(hit)) stop("Cannot find Sphinx R/ sources")
  normalizePath(hit[[1]])
})
for (f in list.files(.sphinx_r, pattern = "[.]R$", full.names = TRUE)) {
  sys.source(f, envir = globalenv(), chdir = FALSE)
}

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(igraph)
  library(ggraph)
  library(scales)
  library(pheatmap)
})

source(file.path(root_dir, "style.R"), local = TRUE)
suppressPackageStartupMessages(library(cowplot))

spatch_dir <- file.path(root_dir, "data", "spatch", "COAD3")
codex_rds <- file.path(spatch_dir, "COAD_codex.rds")
cellchat_rds <- file.path(spatch_dir, "cellchat_coad.rds")
cache_dir <- file.path(root_dir, "results", "Fig3", "cache")
fig3_dir <- file.path(root_dir, "results", "Fig3")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

results_csv <- file.path(cache_dir, "results.csv")
edges_rds <- file.path(cache_dir, "spatial_edges.rds")
interaction_rds <- file.path(cache_dir, "interaction_results.rds")
dist_rds <- file.path(cache_dir, "dist_result.rds")
purity_rds <- file.path(cache_dir, "neighborhood_purity.rds")

celltype_col <- "celltype"
k_cn <- 12
methods_cmp <- c("knn", "radius", "delaunay", "window")
sub_n_methods <- 40000L

save_gg <- save_gg_dir(fig3_dir)
save_hm <- save_hm_dir(fig3_dir)

seurat_to_df <- function(obj, celltype_col = "annotation") {
  emb <- as.data.frame(Embeddings(obj, "spatial"))
  colnames(emb)[1:2] <- c("X", "Y")
  df <- data.frame(
    Cell_ID = colnames(obj),
    X = emb[colnames(obj), "X"],
    Y = emb[colnames(obj), "Y"],
    celltype = as.character(obj@meta.data[[celltype_col]]),
    stringsAsFactors = FALSE
  )
  df[!is.na(df$celltype) & is.finite(df$X) & is.finite(df$Y), , drop = FALSE]
}

plot_distance_heatmap <- function(dist_result) {
  mat <- dist_result$distance_matrix
  dd <- reshape2::melt(mat, varnames = c("Source", "Target"), value.name = "Distance")
  dd$label <- ifelse(is.na(dd$Distance), "", sprintf("%.0f", dd$Distance))
  dist_cols <- c("#2166AC", "#F7F7F7", "#B2182B")
  d_rng <- range(dd$Distance, na.rm = TRUE)
  cols_ramp <- colorRamp(dist_cols)
  dd$txt_col <- vapply(dd$Distance, function(v) {
    if (!is.finite(v)) return("black")
    t <- (v - d_rng[1]) / max(diff(d_rng), 1e-9)
    rgb_v <- as.numeric(cols_ramp(t)) / 255
    lum <- 0.299 * rgb_v[1] + 0.587 * rgb_v[2] + 0.114 * rgb_v[3]
    if (lum >= 0.55) "black" else "white"
  }, character(1))
  p <- ggplot(dd, aes(Source, Target, fill = Distance)) +
    geom_tile(color = "white", linewidth = 0.4) +
    scale_fill_gradientn(
      colours = dist_cols,
      na.value = "grey90",
      name = "Mean\nDistance"
    ) +
    labs(title = NULL, x = "Source Cell Type", y = "Target Cell Type") +
    coord_fixed() +
    panel_theme(12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 11))
  for (i in seq_len(nrow(dd))) {
    p <- p + annotate(
      "text",
      x = dd$Source[i], y = dd$Target[i], label = dd$label[i],
      colour = dd$txt_col[i], size = 3.2
    )
  }
  p
}

plot_spatial_net <- function(df, edges, celltype_col = "celltype",
                                 edge_mode = c("top", "all"), top_n = 4000,
                                 xlim = NULL, ylim = NULL,
                                 point_size = 0.12, edge_alpha = 0.45,
                                 edge_linewidth = 0.22, title = NULL,
                                 color_map = NULL, zoom_box = NULL,
                                 max_cells = NULL) {
  edge_mode <- match.arg(edge_mode)
  df <- as.data.frame(df)
  ed <- as.data.frame(edges)
  if (!is.null(xlim) && !is.null(ylim)) {
    df <- df[df$X >= xlim[1] & df$X <= xlim[2] &
               df$Y >= ylim[1] & df$Y <= ylim[2], , drop = FALSE]
  }
  if (!is.null(max_cells) && nrow(df) > max_cells) {
    set.seed(1234)
    df <- df[sample.int(nrow(df), max_cells), , drop = FALSE]
  }
  id_set <- df$Cell_ID
  ed <- ed[ed$from %in% id_set & ed$to %in% id_set, , drop = FALSE]
  if (edge_mode == "top" && nrow(ed) > top_n) {
    if ("dist" %in% names(ed)) {
      ed <- ed[order(ed$dist), , drop = FALSE][seq_len(top_n), , drop = FALSE]
    } else {
      ed <- ed[sample.int(nrow(ed), top_n), , drop = FALSE]
    }
  }
  xy <- df[, c("X", "Y"), drop = FALSE]
  rownames(xy) <- df$Cell_ID
  seg <- data.frame(
    x = xy[ed$from, "X"], y = xy[ed$from, "Y"],
    xend = xy[ed$to, "X"], yend = xy[ed$to, "Y"]
  )
  seg <- seg[is.finite(seg$x) & is.finite(seg$xend), , drop = FALSE]
  if (is.null(color_map)) color_map <- panel_cols(df[[celltype_col]])

  p <- ggplot() +
    geom_segment(
      data = seg, aes(x = x, y = y, xend = xend, yend = yend),
      color = "#333333", alpha = edge_alpha, linewidth = edge_linewidth,
      lineend = "round"
    ) +
    geom_point(
      data = df, aes(x = X, y = Y, color = .data[[celltype_col]]),
      size = point_size, alpha = 0.75, stroke = 16
    ) +
    scale_color_manual(values = color_map, name = "Cell Type") +
    labs(title = title, x = "X", y = "Y") +
    panel_theme(12) +
    guides(color = guide_legend(override.aes = list(size = 2.4, alpha = 1), ncol = 1))
  if (!is.null(xlim) && !is.null(ylim)) {
    p <- p + coord_fixed(xlim = xlim, ylim = ylim, expand = FALSE, clip = "on", ratio = 1)
  } else {
    p <- p + coord_fixed(expand = FALSE, ratio = 1)
  }
  if (!is.null(zoom_box)) {
    p <- p + annotate(
      "rect",
      xmin = zoom_box[1], xmax = zoom_box[2],
      ymin = zoom_box[3], ymax = zoom_box[4],
      fill = NA, color = "black", linewidth = 0.8
    )
  }
  p
}

find_epi_square <- function(df, celltype = "Epithelial", half = 900) {
  epi <- df[df$celltype == celltype, , drop = FALSE]
  if (nrow(epi) < 50) {
    cx <- mean(range(df$X)); cy <- mean(range(df$Y))
    return(c(cx - half, cx + half, cy - half, cy + half))
  }
  nb <- 35L
  xr <- range(epi$X); yr <- range(epi$Y)
  xi <- pmin(nb, pmax(1L, as.integer(cut(epi$X, breaks = nb, labels = FALSE))))
  yi <- pmin(nb, pmax(1L, as.integer(cut(epi$Y, breaks = nb, labels = FALSE))))
  tab <- table(factor(xi, levels = seq_len(nb)), factor(yi, levels = seq_len(nb)))
  ij <- which(tab == max(tab), arr.ind = TRUE)[1, ]
  cx <- xr[1] + (ij[1] - 0.5) * diff(xr) / nb
  cy <- yr[1] + (ij[2] - 0.5) * diff(yr) / nb
  c(cx - half, cx + half, cy - half, cy + half)
}

plot_pair_distance_ridges <- function(dist_result, pair_list, color_map = NULL,
                                      xmax = NULL) {
  rows <- lapply(names(pair_list), function(lab) {
    src <- pair_list[[lab]][1]; tgt <- pair_list[[lab]][2]
    key <- paste(src, tgt, sep = "_")
    v <- dist_result$distance_distributions[[key]]
    v <- as.numeric(v[is.finite(v)])
    if (!length(v)) return(NULL)
    if (length(v) > 40000L) {
      set.seed(1234)
      v <- sample(v, 40000L)
    }
    data.frame(pair = lab, source = src, distance = v, stringsAsFactors = FALSE)
  })
  dd <- dplyr::bind_rows(rows)
  if (!nrow(dd)) stop("No distance distributions for requested pairs")
  ord <- dd %>%
    dplyr::group_by(pair) %>%
    dplyr::summarise(m = median(distance), .groups = "drop") %>%
    dplyr::arrange(m) %>%
    dplyr::pull(pair)
  dd$pair <- factor(dd$pair, levels = rev(ord))
  if (is.null(xmax) || !is.finite(xmax)) {
    xmax <- as.numeric(stats::quantile(dd$distance, 0.98, na.rm = TRUE))
  }
  dd <- dd %>% dplyr::filter(distance <= xmax)
  fill_vals <- if (!is.null(color_map)) {
    srcs <- vapply(names(pair_list), function(nm) pair_list[[nm]][1], "")
    setNames(unname(color_map[srcs]), names(pair_list))
  } else {
    panel_cols(names(pair_list), semantic = FALSE)
  }
  fill_vals <- fill_vals[as.character(levels(dd$pair))]
  p <- ggplot2::ggplot(dd, ggplot2::aes(x = distance, y = pair, fill = pair))
  if (requireNamespace("ggridges", quietly = TRUE)) {
    p <- p + ggridges::geom_density_ridges(
      alpha = 0.88, scale = 1.05, rel_min_height = 0.01,
      color = "grey30", linewidth = 0.3, show.legend = FALSE
    )
  } else {
    p <- p + ggplot2::geom_violin(
      scale = "width", width = 0.9, alpha = 0.88,
      color = "grey30", linewidth = 0.3, show.legend = FALSE
    )
  }
  med <- dd %>% dplyr::group_by(pair) %>% dplyr::summarise(m = median(distance), .groups = "drop")
  p +
    ggplot2::geom_point(
      data = med, ggplot2::aes(x = m, y = pair),
      inherit.aes = FALSE, shape = 124, size = 3.6, color = "grey15"
    ) +
    ggplot2::scale_fill_manual(values = fill_vals, drop = FALSE) +
    ggplot2::coord_cartesian(xlim = c(0, xmax)) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.01, 0.02))) +
    ggplot2::labs(x = "Nearest-neighbor distance", y = NULL) +
    panel_theme(12) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 10),
      panel.grid.major.x = ggplot2::element_line(colour = "grey92", linewidth = 0.3)
    )
}

## ---- load COAD3 ----
stopifnot(file.exists(codex_rds))
codex.obj <- readRDS(codex_rds)
df0 <- seurat_to_df(codex.obj, celltype_col = "annotation")
df <- prepare_data(
  as.data.table(df0),
  cell_id_col = "Cell_ID", x_col = "X", y_col = "Y", celltype_col = "celltype"
)
ct_map <- panel_cols(df$celltype)

if (file.exists(results_csv) && file.exists(edges_rds) && file.exists(interaction_rds)) {
  clustered_df <- read.csv(results_csv, check.names = FALSE)
  spatial_edges <- readRDS(edges_rds)
  interaction_results <- readRDS(interaction_rds)
  dist_result <- if (file.exists(dist_rds)) readRDS(dist_rds) else {
    d <- calculate_celltype_distances(df, celltype_col = celltype_col)
    saveRDS(d, dist_rds); d
  }
} else {
  dist_result <- calculate_celltype_distances(df, celltype_col = celltype_col)
  spatial_edges <- build_spatial_network(df, method = "auto", celltype_col = celltype_col)
  feature_df <- calculate_neighborhood_features(df, spatial_edges, celltype_col = celltype_col)
  clustered_df <- cluster_neighborhoods(
    feature_df = feature_df, spatial_edges = spatial_edges, method = "kmeans", k = k_cn
  )
  interaction_results <- analyze_spatial_interactions(
    clustered_df, spatial_edges, celltype_col = celltype_col
  )
  write.csv(clustered_df, results_csv, row.names = FALSE)
  saveRDS(spatial_edges, edges_rds)
  saveRDS(interaction_results, interaction_rds)
  saveRDS(dist_result, dist_rds)
}

## record auto method selection (same decision tree as build_spatial_network)
auto_metrics <- .calc_spatial_metrics(df, celltype_col = celltype_col)
auto_sel <- .select_method(auto_metrics, nrow(df))
auto_info <- data.frame(
  n_cells = nrow(df),
  selected_method = auto_sel$method,
  reason = auto_sel$reason,
  radius = if (!is.null(auto_sel$params$radius)) auto_sel$params$radius else NA_real_,
  n_neighbors = if (!is.null(auto_sel$params$n_neighbors)) auto_sel$params$n_neighbors else NA_real_,
  median_1nn = auto_metrics$median_1nn_dist,
  cv_nn_dist = auto_metrics$cv_nn_dist,
  clark_evans_R = auto_metrics$clark_evans_R,
  cat_morans_I = auto_metrics$cat_morans_I,
  stringsAsFactors = FALSE
)
write.csv(auto_info, file.path(fig3_dir, "Fig3_auto_method_selection.csv"), row.names = FALSE)

## ---- a: spatial cell-type distribution ----
save_gg(
  visualize_spatial_distribution(
    df, celltype_col = celltype_col, point_size = 0.25, point_alpha = 0.9,
    color_palette = ct_map, title = NULL
  ) + panel_theme(12),
  "Fig3a_spatial_distribution", width = 5.8, height = 4.8, legend = "Fig3_celltype"
)

## ---- b: mean nearest-neighbor distance heatmap ----
save_gg(
  plot_distance_heatmap(dist_result),
  "Fig3b_distance_heatmap", width = 5.2, height = 4.6, legend = "Fig3b_distance"
)

## ---- c: distance ridges for epi / endo / fib pairs ----
ridge_pairs <- list(
  "Endothelial -> Fibroblast" = c("Endothelial", "Fibroblast"),
  "Fibroblast -> Endothelial" = c("Fibroblast", "Endothelial"),
  "Epithelial -> Endothelial" = c("Epithelial", "Endothelial"),
  "Endothelial -> Epithelial" = c("Endothelial", "Epithelial"),
  "Epithelial -> Fibroblast" = c("Epithelial", "Fibroblast"),
  "Fibroblast -> Epithelial" = c("Fibroblast", "Epithelial")
)
ridge_sum <- data.frame(
  pair = names(ridge_pairs),
  source = vapply(ridge_pairs, `[`, 1, FUN.VALUE = character(1)),
  target = vapply(ridge_pairs, `[`, 2, FUN.VALUE = character(1)),
  mean_nn = vapply(
    ridge_pairs,
    function(pr) as.numeric(dist_result$distance_matrix[pr[1], pr[2]]),
    numeric(1)
  ),
  stringsAsFactors = FALSE
)
ridge_sum <- ridge_sum[order(ridge_sum$mean_nn), , drop = FALSE]
write.csv(ridge_sum, file.path(fig3_dir, "Fig3c_ridge_pairs_summary.csv"), row.names = FALSE)

save_gg(
  plot_pair_distance_ridges(dist_result, ridge_pairs, color_map = ct_map, xmax = 120),
  "Fig3c_distance_ridge_epi_endo_to_fib",
  width = 5.4, height = 4.0,
  legend = NULL
)

## ---- d: CN maps under four neighborhood definitions ----
fig3d_png <- file.path(fig3_dir, "Fig3d_4methods_CN_comparison.png")
fig3d_metric_csv <- file.path(fig3_dir, "Fig3d_method_comparison_metrics.csv")
run_fig3d_maps <- !file.exists(fig3d_png)
run_fig3d_metrics <- !file.exists(fig3d_metric_csv)
if (!run_fig3d_maps && !run_fig3d_metrics) {
} else {
  set.seed(1234)
  n_sub <- if (isTRUE(run_fig3d_maps)) {
    min(sub_n_methods, nrow(df))
  } else {
    min(15000L, nrow(df))
  }
  df_sub <- df %>% sample_n(n_sub)
  spatial_panels <- list()
  metric_rows <- list()
  for (m in methods_cmp) {
    edges_m <- tryCatch(
      build_spatial_network(df_sub, method = m, celltype_col = celltype_col),
      error = function(e) NULL
    )
    if (is.null(edges_m)) next
    feat_m <- calculate_neighborhood_features(df_sub, edges_m, celltype_col = celltype_col)
    if (isTRUE(run_fig3d_maps)) {
      cn_m <- cluster_neighborhoods(
        feature_df = feat_m, spatial_edges = edges_m, method = "kmeans", k = k_cn
      )
      spatial_panels[[m]] <- visualize_spatial_distribution(
        cn_m, celltype_col = "Neighborhood_Cluster",
        point_size = 0.3, point_alpha = 0.9,
        color_palette = panel_cols(cn_m$Neighborhood_Cluster, semantic = FALSE), title = m
      ) + panel_theme(11) + theme(legend.position = "none")
    }
    if (isTRUE(run_fig3d_metrics)) {
      nbr_cols <- grep("^Nbr_", names(feat_m), value = TRUE)
      mean_pur <- if (length(nbr_cols)) {
        mat <- as.matrix(as.data.frame(feat_m)[, nbr_cols, drop = FALSE])
        mean(apply(mat, 1, max), na.rm = TRUE)
      } else {
        NA_real_
      }
      same_frac <- if ("same_type" %in% names(edges_m)) {
        mean(as.logical(edges_m$same_type), na.rm = TRUE)
      } else {
        NA_real_
      }
      metric_rows[[m]] <- data.frame(
        Method = m,
        n_edges = nrow(edges_m),
        same_type_edge_frac = same_frac,
        mean_purity = mean_pur,
        is_auto_selected = identical(m, auto_sel$method),
        stringsAsFactors = FALSE
      )
    }
  }
  if (isTRUE(run_fig3d_maps) && length(spatial_panels) > 0) {
    save_gg(
      wrap_plots(spatial_panels, ncol = 2) +
        plot_annotation(title = NULL),
      "Fig3d_4methods_CN_comparison", width = 7.5, height = 6.5, keep_title = TRUE, legend = NULL
    )
  }
  if (isTRUE(run_fig3d_metrics) && length(metric_rows) > 0) {
    met <- dplyr::bind_rows(metric_rows)
    write.csv(met, fig3d_metric_csv, row.names = FALSE)
  }
}

## ---- e: CN × cell-type composition heatmap ----
comp_df <- calculate_cluster_composition(
  clustered_df, cluster_col = "Neighborhood_Cluster", celltype_col = celltype_col
)
comp_mat <- comp_df %>%
  select(Neighborhood_Cluster, celltype, proportion) %>%
  tidyr::pivot_wider(names_from = celltype, values_from = proportion, values_fill = 0) %>%
  as.data.frame()
rownames(comp_mat) <- comp_mat$Neighborhood_Cluster
prop_mat <- as.matrix(comp_mat[, -1, drop = FALSE])
z_mat <- t(apply(prop_mat, 1, scale))
rownames(z_mat) <- rownames(prop_mat)
colnames(z_mat) <- colnames(prop_mat)
z_mat[is.na(z_mat)] <- 0
cluster_counts <- comp_df %>%
  group_by(Neighborhood_Cluster) %>%
  summarise(total = sum(count), .groups = "drop")
log_counts <- log10(cluster_counts$total[match(rownames(z_mat), cluster_counts$Neighborhood_Cluster)] + 1)
hm_e <- ComplexHeatmap::Heatmap(
  z_mat,
  name = "Z-score",
  col = circlize::colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = grid::gpar(fontsize = 11),
  column_names_gp = grid::gpar(fontsize = 11),
  column_title = NULL,
  row_title = "CN",
  width = grid::unit(75, "mm"),
  height = grid::unit(70, "mm"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    v <- prop_mat[i, j]
    z <- z_mat[i, j]
    txt_col <- if (is.finite(z) && abs(z) > 0.85) "white" else "black"
    grid::grid.text(
      sprintf("%.2f", v), x, y,
      gp = grid::gpar(fontsize = 8, col = txt_col)
    )
  },
  right_annotation = ComplexHeatmap::rowAnnotation(
    `Cell Count` = ComplexHeatmap::anno_barplot(
      log_counts, bar_width = 0.7,
      gp = grid::gpar(fill = "grey70", col = NA),
      width = grid::unit(1.6, "cm")
    )
  ),
  heatmap_legend_param = list(
    title_gp = grid::gpar(fontsize = 11),
    labels_gp = grid::gpar(fontsize = 10),
    legend_height = grid::unit(3, "cm"),
    legend_width = grid::unit(4, "mm")
  )
)
save_hm(
  hm_e, "Fig3e_CN_composition_heatmap",
  width = 9.2, height = 5.6,
  padding = grid::unit(c(5, 5, 5, 18), "mm"),
  legend = "Fig3e_Zscore"
)

## ---- f: neighborhood purity map ----
if (file.exists(purity_rds)) {
  neighborhood_purity <- readRDS(purity_rds)
} else {
  neighborhood_purity <- calculate_neighborhood_purity(
    clustered_df, celltype_col = celltype_col,
    method = "knn", n_neighbors = 10, verbose = FALSE
  )
  saveRDS(neighborhood_purity, purity_rds)
}
save_gg(
  visualize_neighborhood_purity(
    neighborhood_purity, point_size = 0.25, point_alpha = 0.9,
    title = NULL
  ) +
    scale_color_gradientn(colours = c("#2166AC", "#F7F7F7", "#B2182B"), name = "Purity") +
    panel_theme(12),
  "Fig3f_neighborhood_purity", width = 5.5, height = 4.5, legend = "Fig3f_purity"
)

## ---- g: epi- vs stromal-dominant purity ----
strom_frac <- neighborhood_purity$Nbr_Fibroblast + neighborhood_purity$Nbr_Endothelial
grp <- rep(NA_character_, nrow(neighborhood_purity))
grp[neighborhood_purity$Nbr_Epithelial >= 0.5] <- "Epithelial-dominant"
grp[is.na(grp) & strom_frac >= 0.5] <- "Stromal-dominant"
pur_cmp <- data.frame(
  Group = factor(grp, levels = c("Epithelial-dominant", "Stromal-dominant")),
  Purity = neighborhood_purity$Neighborhood_Purity
)
pur_cmp <- pur_cmp[!is.na(pur_cmp$Group), , drop = FALSE]
sum_pur <- pur_cmp %>%
  group_by(Group) %>%
  summarise(
    mean = mean(Purity),
    se = sd(Purity) / sqrt(n()),
    .groups = "drop"
  )
wt <- wilcox.test(Purity ~ Group, data = pur_cmp)
p_lab <- if (wt$p.value < 1e-4) "p < 0.0001" else sprintf("p = %.3f", wt$p.value)
ymax <- max(sum_pur$mean + sum_pur$se) * 1.18
p_pur_bar <- ggplot(sum_pur, aes(Group, mean, fill = Group)) +
  geom_col(width = 0.42, color = NA) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.1, linewidth = 0.4) +
  scale_fill_manual(values = c(
    "Epithelial-dominant" = "#FDB462",
    "Stromal-dominant" = "#80B1D3"
  ), guide = "none") +
  scale_y_continuous(limits = c(0, max(1, ymax)), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  annotate("text", x = 1.5, y = ymax * 0.95, label = p_lab, size = 3.0) +
  labs(title = NULL, x = NULL, y = "Neighborhood Purity") +
  panel_theme(11) +
  theme(
    axis.text.x = element_text(size = 9, angle = 20, hjust = 1),
    plot.margin = margin(6, 8, 4, 6)
  )
save_gg(p_pur_bar, "Fig3g_purity_method_comparison", width = 3.2, height = 3.2, legend = NULL)

## ---- h: CN interaction network ----
save_gg(
  visualize_interaction_network(
    interaction_results$network,
    node_size_range = c(5, 12), edge_size_range = c(0.4, 2.2),
    label_size = 3, show_labels = TRUE, max_nodes = 50
  ) + panel_theme(12),
  "Fig3h_interaction_network", width = 5.5, height = 4.8, legend = NULL
)

## ---- i: spatial network overview + epithelial zoom ----
zoom_box <- find_epi_square(clustered_df, celltype = "Epithelial", half = 420)
xlim_z <- zoom_box[1:2]
ylim_z <- zoom_box[3:4]
side <- max(diff(xlim_z), diff(ylim_z))
cx <- mean(xlim_z); cy <- mean(ylim_z)
xlim_z <- c(cx - side / 2, cx + side / 2)
ylim_z <- c(cy - side / 2, cy + side / 2)
zoom_box <- c(xlim_z, ylim_z)

xr_full <- range(clustered_df$X, na.rm = TRUE)
yr_full <- range(clustered_df$Y, na.rm = TRUE)
asp_full <- diff(yr_full) / max(diff(xr_full), 1e-9)

save_gg(
  plot_spatial_net(
    clustered_df, spatial_edges, celltype_col = celltype_col,
    edge_mode = "top", top_n = 8000,
    point_size = 0.05, edge_alpha = 0.2, edge_linewidth = 0.12,
    title = NULL, color_map = ct_map,
    zoom_box = zoom_box
  ) + theme(aspect.ratio = asp_full),
  "Fig3i_spatial_network_full", width = 8.5, height = max(5.5, 8.5 * asp_full * 0.72 + 1.2), legend = NULL
)
save_gg(
  plot_spatial_net(
    clustered_df, spatial_edges, celltype_col = celltype_col,
    edge_mode = "all",
    xlim = xlim_z, ylim = ylim_z,
    point_size = 0.9, edge_alpha = 0.35, edge_linewidth = 0.28,
    max_cells = 1500,
    title = NULL, color_map = ct_map
  ),
  "Fig3i_spatial_network_celltype", width = 5.2, height = 5.2, legend = NULL
)

## ---- j: CellChat ligand–receptor bubble ----
format_lr_label <- function(x) {
  vapply(x, function(s) {
    parts <- strsplit(s, "_", fixed = TRUE)[[1]]
    if (length(parts) <= 1) return(s)
    if (length(parts) == 2) return(paste(parts[1], parts[2], sep = " - "))
    paste0(parts[1], " - (", paste(parts[-1], collapse = "+"), ")")
  }, character(1), USE.NAMES = FALSE)
}

if (file.exists(cellchat_rds)) {
  if (!isClass("CellChat")) {
    setClass("CellChat", slots = c(
      data = "ANY", data.signaling = "ANY", data.project = "ANY",
      data.scale = "ANY", data.smooth = "ANY", images = "ANY",
      net = "ANY", netP = "ANY", meta = "ANY", idents = "ANY",
      var.features = "ANY", layers = "ANY", LR = "ANY", DB = "ANY", options = "ANY"
    ))
  }
  cellchat <- readRDS(cellchat_rds)
  net <- slot(cellchat, "net")
  prob <- net$prob; pval <- net$pval
  ct_names <- dimnames(prob)[[1]]; lr_names <- dimnames(prob)[[3]]
  extract_pair_df <- function(source, target, pair_order) {
    si <- match(source, ct_names); ti <- match(target, ct_names)
    dd <- data.frame(
      LR_id = lr_names,
      prob = as.numeric(prob[si, ti, ]),
      pval = as.numeric(pval[si, ti, ]),
      stringsAsFactors = FALSE
    )
    dd <- dd[dd$LR_id %in% pair_order & is.finite(dd$prob) & dd$prob > 0 & dd$pval < 0.05, ]
    dd <- dd[match(intersect(pair_order, dd$LR_id), dd$LR_id), ]
    dd$LR_label <- format_lr_label(dd$LR_id)
    dd$LR_label <- factor(dd$LR_label, levels = dd$LR_label[order(dd$prob, decreasing = TRUE)])
    dd$xlab <- paste0(source, " -> ", target)
    dd$panel <- paste0(source, " -> ", target)
    dd
  }
  lr_ef <- c(
    "IHH_PTCH1", "FGF15_FGFR1", "BMP7_ACVR1_BMPR2", "BMP7_BMPR1A_BMPR2",
    "GAS6_AXL", "IGF2_IGF1R", "PDGFA_PDGFRB", "PDGFA_PDGFRA",
    "IGF2_ITGAV_ITGB3", "IGF2_IGF2R"
  )
  lr_fe <- c(
    "IGF2_ITGAV_ITGB3", "IGF2_ITGA6_ITGB4", "ADM_CALCRL", "LGALS9_P4HB",
    "VEGFA_VEGFR2", "PLAU_PLAUR", "POSTN_ITGAV_ITGB5", "POSTN_ITGAV_ITGB3",
    "WNT5A_FZD4", "WNT5A_MCAM"
  )
  df_ef <- extract_pair_df("Epithelial", "Fibroblast", lr_ef)
  df_fe <- extract_pair_df("Fibroblast", "Endothelial", lr_fe)
  prob_lim <- range(c(df_ef$prob, df_fe$prob))
  cc_cols <- c("#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B")
  make_lr_panel <- function(dd) {
    ggplot(dd, aes(x = 1, y = LR_label, color = prob)) +
      geom_point(size = 4.0, shape = 16) +
      scale_color_gradientn(
        colours = cc_cols, limits = prob_lim,
        breaks = prob_lim, labels = c("min", "max"), name = "Commun. Prob."
      ) +
      scale_x_continuous(breaks = 1, labels = unique(dd$panel), expand = c(0, 0)) +
      labs(x = NULL, y = NULL) +
      coord_cartesian(xlim = c(0.85, 1.15), clip = "off") +
      theme_classic(base_size = 9) +
      theme(
        axis.text.x = element_text(size = 7.5, angle = 45, hjust = 1, vjust = 1, colour = "black"),
        axis.text.y = element_text(size = 8, colour = "black"),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
        legend.position = "none",
        plot.margin = margin(2, 1, 2, 1),
        aspect.ratio = 5
      )
  }
  leg_bar <- data.frame(x = 1, y = seq(0, 1, length.out = 100), z = seq(0, 1, length.out = 100))
  p_leg <- ggplot(leg_bar, aes(x, y, fill = z)) +
    geom_raster() +
    scale_fill_gradientn(colours = cc_cols, guide = "none") +
    annotate("text", x = 1, y = 1.12, label = "max", size = 2.5) +
    annotate("text", x = 1, y = -0.12, label = "min", size = 2.5) +
    annotate("text", x = 1, y = 1.28, label = "Commun.\nProb.", size = 3.2, lineheight = 0.9) +
    annotate("text", x = 1.2, y = -0.45, label = "p-value", size = 2.5) +
    annotate("point", x = 0.45, y = -0.72, size = 3.2) +
    annotate("text", x = 2.0, y = -0.72, label = "p < 0.01", size = 2.4, hjust = 0) +
    coord_cartesian(xlim = c(0.2, 3.2), ylim = c(-1.0, 1.45), clip = "off") +
    theme_void() + theme(plot.margin = margin(10, 10, 24, 4))
  p_out <- make_lr_panel(df_ef) + make_lr_panel(df_fe) + p_leg +
    plot_layout(widths = c(0.34, 0.34, 0.48))
  save_gg(p_out, "Fig3j_CellChat_LR_bubble", width = 4.8, height = 4.6, legend = NULL)
}
