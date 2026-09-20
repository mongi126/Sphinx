#!/usr/bin/env Rscript
## Fig2 LUAD-3
## a  UMAP by cluster
## d  Marker bubble plot
## b  UMAP by cell type
## c  Spatial cell-type map


rm(list = ls())
set.seed(1234)
options(bitmapType = "cairo")

.root <- local({
  ca <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", ca[grep("^--file=", ca)])
  if (length(f) && nzchar(f[[1]])) dirname(normalizePath(f[[1]])) else normalizePath(getwd())
})
.sphinx_r <- local({
  cands <- c(file.path(dirname(.root), "R"), file.path(.root, "Sphinx", "R"))
  hit <- cands[dir.exists(cands)]
  if (!length(hit)) stop("Cannot find Sphinx R/ sources")
  normalizePath(hit[[1]])
})
for (f in list.files(.sphinx_r, pattern = "[.]R$", full.names = TRUE)) {
  sys.source(f, envir = globalenv(), chdir = FALSE)
}

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
  library(grid)
})

source(file.path(.root, "style.R"), local = TRUE)
suppressPackageStartupMessages(library(cowplot))

data_dir <- file.path(.root, "data")
out_dir <- file.path(.root, "results", "Fig2")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

save_gg <- save_gg_dir(out_dir)

obj_rds <- file.path(data_dir, "luad3_celltype.rds")
stopifnot(file.exists(obj_rds))
codex.obj <- readRDS(obj_rds)

if (!all(c("X", "Y") %in% colnames(codex.obj@meta.data))) {
  codex.obj <- extract_spatial_coordinates(codex.obj)
}

assay_use <- if ("Akoya" %in% Assays(codex.obj)) "Akoya" else DefaultAssay(codex.obj)
DefaultAssay(codex.obj) <- assay_use

## ---- a. UMAP by cluster ----
cl_levels <- sort(unique(as.character(codex.obj$seurat_clusters)))
cl_cols <- panel_cols(cl_levels)
p_cluster <- DimPlot(
  codex.obj,
  reduction = "umap",
  group.by = "seurat_clusters",
  cols = unname(cl_cols[cl_levels]),
  label = TRUE,
  repel = TRUE,
  pt.size = 0.3
) + theme(aspect.ratio = 1) + ggtitle(NULL)
save_gg(p_cluster, "Fig2a_umap_cluster", width = 8, height = 7, legend = "Fig2a_cluster")

## ---- d. Marker bubble plot ----
marker_groups <- list(
  TAM = c("CD163", "TIGIT", "Collagen-IV", "TP53"),
  `GC B` = c("CD20", "CXCL13", "HLA-A", "HLA-DPB1"),
  `PD-L1 DC` = c("Podoplanin", "PD-L1", "DC-LAMP", "Ki67", "E-Cadherin"),
  `FOXP3 Tumor` = c("CD8", "LIF", "CD56", "CTLA4", "CD68", "CD44"),
  `CD8 T` = c("CD45RO", "CD4", "CD3e", "b-Actin"),
  APC = c("CD11c"),
  Endo = c("CD31", "b-Catenin1"),
  `Gal3 DC` = c("Galectin3", "LAG3"),
  CAF = c("a-SMA"),
  Tumor = c("Pan-Cytokeratin", "NKX2-1")
)

group_colors <- panel_cols(names(marker_groups))

cluster_to_group <- c(
  "0" = "TAM", "5" = "GC B", "6" = "GC B", "1" = "PD-L1 DC", "11" = "PD-L1 DC",
  "2" = "FOXP3 Tumor", "3" = "CD8 T", "4" = "APC", "7" = "Endo",
  "8" = "Gal3 DC", "9" = "CAF", "10" = "Tumor"
)

features <- unlist(marker_groups, use.names = FALSE)
available <- rownames(codex.obj[[assay_use]])
features <- intersect(features, available)
marker_groups <- lapply(marker_groups, function(x) intersect(x, available))
marker_groups <- marker_groups[vapply(marker_groups, length, 1L) > 0]

feature_to_group <- unlist(lapply(names(marker_groups), function(g) {
  setNames(rep(g, length(marker_groups[[g]])), marker_groups[[g]])
}))

cluster_order <- c("0", "5", "6", "1", "11", "2", "3", "4", "7", "8", "9", "10")
Idents(codex.obj) <- factor(as.character(codex.obj$seurat_clusters), levels = cluster_order)

expr_mat <- as.matrix(GetAssayData(codex.obj, assay = assay_use, layer = "data")[features, , drop = FALSE])
clusters <- as.character(codex.obj$seurat_clusters)
feat_thr <- apply(expr_mat, 1, function(x) stats::median(x, na.rm = TRUE))

plot_list <- lapply(cluster_order, function(cl) {
  in_cl <- clusters == cl
  if (!any(in_cl)) return(NULL)
  vapply(features, function(f) {
    x <- expr_mat[f, ]
    mu_in <- mean(x[in_cl], na.rm = TRUE)
    mu_out <- mean(x[!in_cl], na.rm = TRUE)
    pct <- mean(x[in_cl] > feat_thr[[f]], na.rm = TRUE) * 100
    log2fc <- log2((mu_in + 1e-6) / (mu_out + 1e-6))
    c(avg_log2FC = log2fc, pct.exp = pct, avg_in = mu_in)
  }, numeric(3)) %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("features.plot") %>%
    mutate(id = cl)
})
plot_df <- bind_rows(plot_list)
plot_df$features.plot <- factor(plot_df$features.plot, levels = features)
plot_df$id <- factor(plot_df$id, levels = cluster_order)
plot_df$x <- as.numeric(plot_df$features.plot)
plot_df$y <- as.numeric(plot_df$id)
plot_df$marker_group <- feature_to_group[as.character(plot_df$features.plot)]
plot_df$cluster_group <- cluster_to_group[as.character(plot_df$id)]

plot_df <- plot_df %>%
  group_by(features.plot) %>%
  mutate(
    rank_fc = rank(-avg_log2FC, ties.method = "min"),
    keep_extra = rank_fc <= 2 & pct.exp >= 30 & avg_log2FC >= 0.35
  ) %>%
  ungroup() %>%
  filter(marker_group == cluster_group | keep_extra)

plot_df$pct.exp <- pmin(100, pmax(plot_df$pct.exp, 0))

anno_df <- do.call(rbind, lapply(names(marker_groups), function(g) {
  data.frame(group = g, feature = marker_groups[[g]], stringsAsFactors = FALSE)
}))
anno_df$feature <- factor(anno_df$feature, levels = features)
anno_df$x <- as.numeric(anno_df$feature)
group_pos <- anno_df %>%
  group_by(group) %>%
  summarise(xmin = min(x) - 0.5, xmax = max(x) + 0.5, xmid = mean(x), .groups = "drop")
group_pos$group <- factor(group_pos$group, levels = names(marker_groups))

y_max <- length(cluster_order)
bar_ymin <- y_max + 0.85
bar_ymax <- y_max + 1.55

p_bubble <- ggplot() +
  geom_point(
    data = plot_df,
    aes(x = x, y = y, size = pct.exp, color = avg_log2FC),
    shape = 16
  ) +
  geom_rect(
    data = group_pos,
    aes(xmin = xmin, xmax = xmax, ymin = bar_ymin, ymax = bar_ymax, fill = group),
    color = NA
  ) +
  geom_text(
    data = group_pos,
    aes(x = xmid, y = (bar_ymin + bar_ymax) / 2, label = group),
    size = 3.2, fontface = "bold", color = "black"
  ) +
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_color_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 1.0, name = "Avg log2FC",
    limits = c(0.5, 2.0), oob = scales::squish,
    guide = guide_colorbar(
      title.position = "top", title.hjust = 0.5,
      barwidth = unit(0.45, "cm"), barheight = unit(3.2, "cm"), order = 1
    )
  ) +
  scale_size_continuous(
    range = c(1.5, 9), limits = c(0, 100),
    breaks = c(25, 50, 75, 100), name = "Percent Expressed",
    guide = guide_legend(
      title.position = "top", title.hjust = 0.5,
      override.aes = list(color = "black"), order = 2
    )
  ) +
  scale_x_continuous(
    breaks = seq_along(features), labels = features,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    breaks = seq_along(cluster_order), labels = cluster_order, expand = c(0, 0)
  ) +
  coord_cartesian(ylim = c(0.5, y_max + 0.5), clip = "off") +
  theme_bw(base_size = 13) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.6),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    plot.margin = margin(36, 12, 8, 8)
  )
save_gg(p_bubble, "Fig2d_marker_bubble_plot", width = 14, height = 7.2, legend = "Fig2d_marker")

## ---- Annotate ----
cluster_ids <- as.character(0:11)
celltype_labels <- c(
  "TAM", "PD-L1_Ki67_DC", "FOXP3_proliferate_Tumor", "CD8_T",
  "APC", "GC_B", "GC_B", "Endothelial",
  "Galectin3_DC", "CAF", "Tumor", "PD-L1_Ki67_DC"
)
codex.obj <- annotate_celltypes(
  codex.obj,
  cluster_ids = cluster_ids,
  celltype_labels = celltype_labels,
  cluster_column = "Akoya_snn_res.0.5"
)

## ---- b. UMAP by cell type ----
ct_levels <- sort(unique(na.omit(as.character(codex.obj$celltype))))
ct_cols <- panel_cols(ct_levels)
codex.obj$celltype <- factor(codex.obj$celltype, levels = ct_levels)
p_ct <- DimPlot(
  codex.obj, reduction = "umap", group.by = "celltype",
  cols = unname(ct_cols[ct_levels]), pt.size = 0.5
) + theme(aspect.ratio = 1) + ggtitle(NULL) + panel_theme(12)
save_gg(p_ct, "Fig2b_umap_celltype", width = 6.5, height = 5.8, legend = "Fig2_celltype")

## ---- c. Spatial distribution ----
md <- codex.obj@meta.data
p_sp <- ggplot(md, aes(x = X, y = Y, color = celltype)) +
  geom_point(size = 0.4, alpha = 0.85) +
  scale_color_manual(values = ct_cols, name = "Cell Type") +
  coord_fixed() +
  labs(x = "X", y = "Y", title = NULL) +
  panel_theme(12)
save_gg(p_sp, "Fig2c_spatial_distribution", width = 6.8, height = 5.8, legend = NULL)  # shared Fig2_celltype
