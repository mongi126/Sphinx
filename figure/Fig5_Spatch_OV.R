#!/usr/bin/env Rscript
## Fig5 Spatch OV — VisiumHD (OV-7) and Xenium (OV-1)
## Spatial cell-type / CN maps, fibroblast-CN focus, composition correlation


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
  library(stringr)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
  library(ggrepel)
  library(grid)
})

source(file.path(root_dir, "style.R"), local = TRUE)
suppressPackageStartupMessages(library(cowplot))

out_dir <- file.path(root_dir, "results", "Fig5")
cache_dir <- file.path(out_dir, "cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

visium_rds <- file.path(root_dir, "data", "spatch", "OV7_VisiumHD", "visiumhd.rds")
xenium_rds <- file.path(root_dir, "data", "spatch", "OV1_Xenium", "xenium.rds")
k_cn <- 12
sub_n <- 40000L

save_gg <- save_gg_dir(out_dir)
save_hm <- save_hm_dir(out_dir)

seurat_to_df <- function(obj, celltype_col = "annotation") {
  emb <- as.data.frame(Embeddings(obj, "spatial"))
  colnames(emb)[1:2] <- c("X", "Y")
  meta <- obj@meta.data
  df <- data.frame(
    Cell_ID = colnames(obj),
    X = emb[colnames(obj), "X"],
    Y = emb[colnames(obj), "Y"],
    celltype = as.character(meta[[celltype_col]]),
    stringsAsFactors = FALSE
  )
  if ("high_quality" %in% names(meta)) {
    hq <- meta$high_quality[match(df$Cell_ID, colnames(obj))]
    df <- df[!is.na(hq) & hq == 1, , drop = FALSE]
  }
  df <- df[!is.na(df$celltype) & nzchar(df$celltype) & is.finite(df$X) & is.finite(df$Y), , drop = FALSE]
  df
}

load_or_cluster <- function(rds_path, platform, results_path, edges_path) {
  if (file.exists(results_path) && file.exists(edges_path)) {
    clustered <- read.csv(results_path, check.names = FALSE)
    edges <- readRDS(edges_path)
    return(list(df = clustered, edges = edges, platform = platform))
  }
  stopifnot(file.exists(rds_path))
  obj <- readRDS(rds_path)
  raw <- seurat_to_df(obj, celltype_col = "annotation")
  if (nrow(raw) > sub_n) {
    set.seed(1234)
    raw <- raw[sample.int(nrow(raw), sub_n), , drop = FALSE]
  }
  df <- prepare_data(
    as.data.table(raw),
    cell_id_col = "Cell_ID", x_col = "X", y_col = "Y", celltype_col = "celltype"
  )
  edges <- build_spatial_network(df, method = "auto", celltype_col = "celltype", verbose = FALSE)
  feat <- calculate_neighborhood_features(df, edges, celltype_col = "celltype")
  clustered <- cluster_neighborhoods(
    feature_df = feat, spatial_edges = edges, method = "kmeans", k = k_cn
  )
  write.csv(clustered, results_path, row.names = FALSE)
  saveRDS(edges, edges_path)
  list(df = clustered, edges = edges, platform = platform)
}

## CN with highest Fibroblast fraction (tie-break: larger n)
pick_fibroblast_cn <- function(df, fib_labels = c("Fibroblast", "CAF")) {
  tab <- df %>%
    mutate(is_fib = celltype %in% fib_labels) %>%
    group_by(Neighborhood_Cluster) %>%
    summarise(n = n(), fib_prop = mean(is_fib), .groups = "drop") %>%
    arrange(desc(fib_prop), desc(n))
  as.character(tab$Neighborhood_Cluster[[1]])
}

get_cn_composition <- function(df, cn, nbr_prefix = "Nbr_") {
  nbr_cols <- names(df)[startsWith(names(df), nbr_prefix)]
  stopifnot(length(nbr_cols) > 0)
  sub <- df %>% filter(.data$Neighborhood_Cluster == cn)
  if (nrow(sub) == 0) return(NULL)
  sub %>%
    summarise(across(all_of(nbr_cols), ~ mean(.x, na.rm = TRUE))) %>%
    pivot_longer(everything(), names_to = "celltype", values_to = "fraction") %>%
    mutate(celltype = str_remove(.data$celltype, paste0("^", nbr_prefix))) %>%
    arrange(desc(.data$fraction))
}

## CN × celltype composition heatmap (Z-score; classic red–blue)
draw_cn_comp_heatmap <- function(df, panel_title = NULL) {
  comp_df <- calculate_cluster_composition(
    df, cluster_col = "Neighborhood_Cluster", celltype_col = "celltype"
  )
  prop_wide <- comp_df %>%
    select(Neighborhood_Cluster, celltype, proportion) %>%
    tidyr::pivot_wider(names_from = celltype, values_from = proportion, values_fill = 0) %>%
    as.data.frame()
  rownames(prop_wide) <- as.character(prop_wide$Neighborhood_Cluster)
  prop_mat <- as.matrix(prop_wide[, setdiff(colnames(prop_wide), "Neighborhood_Cluster"), drop = FALSE])
  storage.mode(prop_mat) <- "numeric"
  ## numeric CN row order before clustering for stable labels
  cn_lab <- rownames(prop_mat)
  cn_num <- suppressWarnings(as.integer(sub("^CN", "", cn_lab)))
  prop_mat <- prop_mat[order(cn_num, cn_lab), , drop = FALSE]

  z_mat <- t(apply(prop_mat, 1, scale))
  if (is.null(dim(z_mat))) z_mat <- matrix(z_mat, nrow = 1)
  rownames(z_mat) <- rownames(prop_mat)
  colnames(z_mat) <- colnames(prop_mat)
  z_mat[is.na(z_mat)] <- 0

  cluster_counts <- comp_df %>%
    group_by(Neighborhood_Cluster) %>%
    summarise(total = sum(count), .groups = "drop")
  log_counts <- log10(
    cluster_counts$total[match(rownames(z_mat), as.character(cluster_counts$Neighborhood_Cluster))] + 1
  )

  ComplexHeatmap::Heatmap(
    z_mat,
    name = "Z-score",
    col = circlize::colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 10),
    column_names_gp = grid::gpar(fontsize = 9),
    column_names_rot = 45,
    column_title = panel_title,
    column_title_gp = grid::gpar(fontsize = 10, fontface = "plain"),
    row_title = "CN",
    row_title_gp = grid::gpar(fontsize = 10),
    width = grid::unit(72, "mm"),
    height = grid::unit(64, "mm"),
    cell_fun = function(j, i, x, y, width, height, fill) {
      v <- prop_mat[rownames(z_mat)[i], colnames(z_mat)[j]]
      z <- z_mat[i, j]
      txt_col <- if (is.finite(z) && abs(z) > 0.85) "white" else "black"
      grid::grid.text(
        sprintf("%.2f", v), x, y,
        gp = grid::gpar(fontsize = 7, fontface = "bold", col = txt_col)
      )
    },
    right_annotation = ComplexHeatmap::rowAnnotation(
      `Cell Count` = ComplexHeatmap::anno_barplot(
        log_counts, bar_width = 0.7,
        gp = grid::gpar(fill = "grey70", col = NA),
        width = grid::unit(1.15, "cm")
      )
    ),
    heatmap_legend_param = list(
      title_gp = grid::gpar(fontsize = 10),
      labels_gp = grid::gpar(fontsize = 9),
      legend_height = grid::unit(2.6, "cm"),
      legend_width = grid::unit(3.5, "mm")
    )
  )
}

vis <- load_or_cluster(
  visium_rds, "VisiumHD",
  file.path(cache_dir, "visiumhd_results.csv"),
  file.path(cache_dir, "visiumhd_edges.rds")
)
xen <- load_or_cluster(
  xenium_rds, "Xenium",
  file.path(cache_dir, "xenium_results.csv"),
  file.path(cache_dir, "xenium_edges.rds")
)

## Fibroblast-dominant CN per platform
visiumhd_cn <- pick_fibroblast_cn(vis$df)
xenium_cn <- pick_fibroblast_cn(xen$df)
write.csv(
  data.frame(
    Platform = c("VisiumHD", "Xenium"),
    Fibroblast_CN = c(visiumhd_cn, xenium_cn),
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "Fig5_selected_fibroblast_CN.csv"),
  row.names = FALSE
)

## shared celltype colors across platforms
all_ct <- sort(unique(c(vis$df$celltype, xen$df$celltype)))
ct_map <- panel_cols(all_ct)

## a Visium HD
save_gg(
  visualize_spatial_distribution(
    vis$df, celltype_col = "celltype", point_size = 0.3, point_alpha = 0.9,
    color_palette = ct_map, title = NULL
  ) + panel_theme(12),
  "Fig5a_VisiumHD_spatial_celltype", width = 5.8, height = 4.8, legend = "Fig5_Visium_celltype"
)
save_gg(
  visualize_spatial_distribution(
    vis$df, celltype_col = "Neighborhood_Cluster", point_size = 0.3, point_alpha = 0.9,
    color_palette = panel_cols(vis$df$Neighborhood_Cluster, semantic = FALSE),
    title = NULL
  ) + panel_theme(12),
  "Fig5a_VisiumHD_spatial_CN", width = 5.8, height = 4.8, legend = "Fig5_Visium_CN"
)

## b Xenium
save_gg(
  visualize_spatial_distribution(
    xen$df, celltype_col = "celltype", point_size = 0.3, point_alpha = 0.9,
    color_palette = ct_map, title = NULL
  ) + panel_theme(12),
  "Fig5b_Xenium_spatial_celltype", width = 5.8, height = 4.8, legend = "Fig5_Xenium_celltype"
)
save_gg(
  visualize_spatial_distribution(
    xen$df, celltype_col = "Neighborhood_Cluster", point_size = 0.3, point_alpha = 0.9,
    color_palette = panel_cols(xen$df$Neighborhood_Cluster, semantic = FALSE),
    title = NULL
  ) + panel_theme(12),
  "Fig5b_Xenium_spatial_CN", width = 5.8, height = 4.8, legend = "Fig5_Xenium_CN"
)

## c composition heatmaps (compact; bold in-cell labels; platform subtitle)
save_hm(
  draw_cn_comp_heatmap(vis$df, panel_title = "Visium HD CN"),
  "Fig5c_CN_composition_heatmap_VisiumHD",
  width = 7.0, height = 5.0,
  padding = grid::unit(c(4, 4, 4, 12), "mm"),
  legend = "Fig5c_Zscore"
)
save_hm(
  draw_cn_comp_heatmap(xen$df, panel_title = "Xenium CN"),
  "Fig5c_CN_composition_heatmap_Xenium",
  width = 7.0, height = 5.0,
  padding = grid::unit(c(4, 4, 4, 12), "mm"),
  legend = NULL
)
unlink(list.files(out_dir, pattern = "^Fig5c_composition_(VisiumHD|Xenium)\\.(pdf|png)$", full.names = TRUE))

## d Voronoi for fibroblast-dominant CN (platform-tagged filenames)
vor_jobs <- list(
  list(tag = "VisiumHD", cn = visiumhd_cn, df = vis$df),
  list(tag = "Xenium", cn = xenium_cn, df = xen$df)
)
unlink(list.files(out_dir, pattern = "^Fig5d_voronoi_CN[0-9]+\\.(pdf|png)$", full.names = TRUE))

for (i in seq_along(vor_jobs)) {
  job <- vor_jobs[[i]]
  cn <- job$cn
  src <- job$df
  if (!cn %in% unique(as.character(src$Neighborhood_Cluster))) {
    next
  }
  src_plot <- src
  if (nrow(src_plot) > 15000) {
    set.seed(1234)
    src_plot <- src_plot %>% sample_n(15000)
  }
  p_v <- visualize_voronoi(
    src_plot, coloring = "neighborhood",
    celltype_col = "celltype", neighborhood_col = "Neighborhood_Cluster",
    highlight_cluster = cn, show_composition = TRUE
  ) + ggtitle(paste0(job$tag, "  ", cn)) + panel_theme(12) +
    theme(plot.title = element_text(hjust = 0.5, size = 13, face = "plain"))
  save_gg(
    p_v,
    paste0("Fig5d_voronoi_", job$tag, "_", cn),
    width = 6.5, height = 5.5, keep_title = TRUE,
    legend = if (i == 1L) "Fig5_voronoi_celltype" else NULL
  )
}

## e fibroblast-CN composition correlation (Xenium vs VisiumHD)
comp_x <- get_cn_composition(xen$df, xenium_cn)
comp_v <- get_cn_composition(vis$df, visiumhd_cn)
unlink(list.files(out_dir, pattern = "^Fig5e_composition_CN.*_vs_CN.*\\.csv$", full.names = TRUE))
if (!is.null(comp_x) && !is.null(comp_v)) {
  merged <- full_join(
    comp_x %>% rename(xenium_fraction = fraction),
    comp_v %>% rename(visiumhd_fraction = fraction),
    by = "celltype"
  ) %>%
    mutate(
      xenium_fraction = coalesce(xenium_fraction, 0),
      visiumhd_fraction = coalesce(visiumhd_fraction, 0),
      total_fraction = xenium_fraction + visiumhd_fraction
    ) %>%
    arrange(desc(total_fraction))
  pearson_r <- cor(merged$xenium_fraction, merged$visiumhd_fraction, method = "pearson")
  spearman_r <- cor(merged$xenium_fraction, merged$visiumhd_fraction, method = "spearman")
  write.csv(
    merged,
    file.path(out_dir, paste0("Fig5e_composition_", xenium_cn, "_vs_", visiumhd_cn, ".csv")),
    row.names = FALSE
  )
  p_sc <- ggplot(merged, aes(xenium_fraction, visiumhd_fraction, label = celltype, color = total_fraction)) +
    geom_point(size = 2.4, alpha = 0.95) +
    geom_smooth(aes(group = 1), method = "lm", se = TRUE, linewidth = 0.7, color = "grey20") +
    ggrepel::geom_text_repel(size = 3.5, max.overlaps = 30, min.segment.length = 0, show.legend = FALSE) +
    scale_color_gradient(low = "#2166AC", high = "#B2182B", name = "Total\nfraction") +
    labs(
      title = NULL,
      subtitle = paste0(
        "Xenium ", xenium_cn, " vs VisiumHD ", visiumhd_cn,
        "  |  r=", signif(pearson_r, 3), "  rho=", signif(spearman_r, 3)
      ),
      x = paste0("Xenium ", xenium_cn, " mean Nbr fraction"),
      y = paste0("VisiumHD ", visiumhd_cn, " mean Nbr fraction")
    ) +
    panel_theme(12) +
    theme(plot.subtitle = element_text(hjust = 0.5, size = 11))
  save_gg(p_sc, "Fig5e_composition_correlation", width = 5.8, height = 4.6, keep_title = TRUE, legend = "Fig5e_correlation")
}
