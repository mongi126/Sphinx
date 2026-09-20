#!/usr/bin/env Rscript
## Fig4 TSU-33 (CODEX)
## a  Spatial cell-type map
## b  CN composition Z-score heatmap
## d  Within-CN heterotypic contact enrichment
## e  Differential marker volcano (CN6 / CN10 / CN27)
## f  Pathway enrichment bars (top terms by CN)
## g  Kaplan–Meier survival by CN signature (High vs Low)


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
  library(forcats)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(patchwork)
  library(ggrepel)
  library(pheatmap)
})

source(file.path(root_dir, "style.R"), local = TRUE)
suppressPackageStartupMessages(library(cowplot))

data_dir <- file.path(root_dir, "data")
out_dir <- file.path(root_dir, "results", "Fig4")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

obj_rds <- file.path(data_dir, "tsu33_celltype.rds")
results_src <- file.path(data_dir, "tsu33_results.csv")
diff_src <- file.path(data_dir, "tsu33_diff_results.csv")
enrich_src <- file.path(data_dir, "tsu33_cluster_enrich.rds")
expr_csv <- file.path(data_dir, "tsu33_expression_filtered.csv")
interaction_src <- file.path(data_dir, "tsu33_interaction_results.rds")
edges_src <- file.path(data_dir, "tsu33_spatial_edges.rds")

results_csv <- file.path(out_dir, "tsu33_results.csv")
edges_rds <- file.path(out_dir, "tsu33_spatial_edges.rds")
interaction_rds <- file.path(out_dir, "tsu33_interaction_results.rds")
cn_rdata <- file.path(data_dir, "tsu33_cn.Rdata")
surv_csv <- file.path(data_dir, "surv_data.csv")
if (!file.exists(surv_csv)) {
  alt <- file.path(out_dir, "surv_data.csv")
  if (file.exists(alt)) surv_csv <- alt
}

target_cn <- c("CN6", "CN10", "CN27")
k_cn <- 30

save_gg <- save_gg_dir(out_dir)
save_hm <- save_hm_dir(out_dir)

draw_cn_heatmap <- function(meta, highlight_cn = NULL) {
  comp_raw <- meta %>%
    filter(!is.na(Neighborhood_Cluster), !is.na(celltype)) %>%
    group_by(Neighborhood_Cluster, celltype) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(Neighborhood_Cluster) %>%
    mutate(prop = count / sum(count)) %>%
    ungroup()

  all_cn <- sort(unique(as.character(comp_raw$Neighborhood_Cluster)))
  all_type <- sort(unique(as.character(comp_raw$celltype)))
  template <- expand.grid(
    Neighborhood_Cluster = all_cn, celltype = all_type, stringsAsFactors = FALSE
  )
  wide <- template %>%
    left_join(comp_raw, by = c("Neighborhood_Cluster", "celltype")) %>%
    mutate(prop = ifelse(is.na(prop), 0, prop)) %>%
    select(Neighborhood_Cluster, celltype, prop) %>%
    pivot_wider(names_from = celltype, values_from = prop) %>%
    as.data.frame()
  rownames(wide) <- as.character(wide$Neighborhood_Cluster)
  prop_mat <- as.matrix(wide[, setdiff(colnames(wide), "Neighborhood_Cluster"), drop = FALSE])
  storage.mode(prop_mat) <- "numeric"

  z_mat <- t(apply(prop_mat, 1, scale))
  if (is.null(dim(z_mat))) z_mat <- matrix(z_mat, nrow = 1)
  rownames(z_mat) <- rownames(prop_mat)
  colnames(z_mat) <- colnames(prop_mat)
  z_mat[is.na(z_mat)] <- 0

  row_hc <- hclust(dist(z_mat))
  col_hc <- hclust(dist(t(z_mat)))
  ro <- row_hc$order
  co <- col_hc$order
  z_mat <- z_mat[ro, co, drop = FALSE]
  prop_show <- prop_mat[ro, co, drop = FALSE]

  cluster_counts <- comp_raw %>%
    group_by(Neighborhood_Cluster) %>%
    summarise(total = sum(count), .groups = "drop")
  log_counts <- log10(
    cluster_counts$total[match(rownames(prop_show), as.character(cluster_counts$Neighborhood_Cluster))] + 1
  )

  col_fun <- circlize::colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))
  row_ha <- ComplexHeatmap::rowAnnotation(
    `Cell Count\n(log10+1)` = ComplexHeatmap::anno_barplot(
      log_counts, bar_width = 0.7,
      gp = grid::gpar(fill = "grey55", col = NA),
      which = "row", width = grid::unit(1.8, "cm")
    )
  )

  ComplexHeatmap::Heatmap(
    z_mat,
    name = "Z-score",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 11, fontface = "plain"),
    column_names_gp = grid::gpar(fontsize = 11, fontface = "plain"),
    column_title = NULL,
    row_title = "CN",
    row_title_gp = grid::gpar(fontsize = 12, fontface = "plain"),
    right_annotation = row_ha,
    heatmap_legend_param = list(
      title_gp = grid::gpar(fontsize = 11, fontface = "plain"),
      labels_gp = grid::gpar(fontsize = 10)
    ),
    cell_fun = function(j, i, x, y, w, h, fill) {
      grid::grid.text(
        sprintf("%.2f", prop_show[i, j]), x, y,
        gp = grid::gpar(fontsize = 8, fontface = "plain")
      )
    }
  )
}

plot_volcano_vertical <- function(diff_sub, diff_thresh = 0.15, p_thresh = 0.05, y_cap = 50,
                                  base_size = 9, label_size = 2.2, point_size = 1.35,
                                  n_up = 3L, n_down = 3L, cluster_levels = NULL) {
  ## Volcano: top 3 up / top 3 down markers per CN.
  df <- diff_sub %>%
    mutate(
      Significance = case_when(
        adj.p.value < p_thresh & MeanDiff >= diff_thresh ~ "Up",
        adj.p.value < p_thresh & MeanDiff <= -diff_thresh ~ "Down",
        TRUE ~ "NS"
      ),
      neglog10p = -log10(pmax(adj.p.value, .Machine$double.xmin)),
      neglog10p_plot = pmin(neglog10p, y_cap)
    )
  if (is.null(cluster_levels)) {
    cluster_levels <- unique(as.character(df$Cluster))
  }
  df$Cluster <- factor(as.character(df$Cluster), levels = cluster_levels)
  label_df <- df %>%
    filter(Significance %in% c("Up", "Down")) %>%
    group_by(Cluster) %>%
    group_modify(~ {
      d <- .x %>% arrange(desc(abs(MeanDiff)), adj.p.value)
      up <- d %>% filter(Significance == "Up") %>% slice_head(n = n_up)
      down <- d %>% filter(Significance == "Down") %>% slice_head(n = n_down)
      bind_rows(up, down)
    }) %>%
    ungroup() %>%
    mutate(
      Protein_lab = dplyr::recode(
        as.character(Protein),
        `Pan-Cytokeratin` = "Pan-CK",
        `b-Catenin1` = "b-Cat1",
        `Podoplanin` = "PDPN",
        `Caveolin` = "Caveolin",
        .default = as.character(Protein)
      )
    )
  x_rng <- df %>%
    group_by(Cluster) %>%
    summarise(
      xmin = min(MeanDiff, na.rm = TRUE),
      xmax = max(MeanDiff, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(xspan = pmax(xmax - xmin, 1e-6))
  ## Label positions in the mid plot band.
  place_side <- function(labs, side) {
    if (nrow(labs) == 0) return(labs)
    fx <- if (side == "Down") c(0.10, 0.25, 0.38) else c(0.62, 0.76, 0.90)
    fy <- c(0.58, 0.44, 0.30)
    labs %>%
      left_join(x_rng, by = "Cluster") %>%
      group_by(Cluster) %>%
      arrange(MeanDiff, .by_group = TRUE) %>%
      mutate(
        idx = pmin(dplyr::row_number(), length(fx)),
        lab_x = xmin + fx[idx] * xspan,
        lab_y = fy[idx] * y_cap
      ) %>%
      ungroup()
  }
  lab_all <- bind_rows(
    place_side(label_df %>% filter(Significance == "Down"), "Down"),
    place_side(label_df %>% filter(Significance == "Up"), "Up")
  )
  ggplot(df, aes(MeanDiff, neglog10p_plot, color = Significance)) +
    geom_point(size = point_size, alpha = 0.9) +
    scale_color_manual(
      values = c(Up = "#D55E00", Down = "#0072B2", NS = "#BDBDBD"),
      breaks = c("Up", "Down", "NS")
    ) +
    geom_vline(
      xintercept = c(-diff_thresh, diff_thresh),
      linetype = "dashed", color = "grey50", linewidth = 0.35
    ) +
    geom_hline(
      yintercept = -log10(p_thresh),
      linetype = "dotted", color = "grey50", linewidth = 0.35
    ) +
    geom_segment(
      data = lab_all,
      aes(x = MeanDiff, y = neglog10p_plot, xend = lab_x, yend = lab_y),
      color = "grey60", linewidth = 0.22, alpha = 0.85,
      inherit.aes = FALSE
    ) +
    geom_label(
      data = lab_all,
      aes(x = lab_x, y = lab_y, label = Protein_lab, color = Significance),
      size = label_size, hjust = 0.5, vjust = 0.5,
      fill = alpha("white", 0.95), label.size = 0,
      label.padding = unit(0.07, "lines"),
      show.legend = FALSE, inherit.aes = FALSE
    ) +
    coord_cartesian(ylim = c(0, y_cap + 2), clip = "on") +
    scale_y_continuous(breaks = seq(0, y_cap, by = 10), expand = expansion(mult = c(0.02, 0.02))) +
    scale_x_continuous(expand = expansion(mult = c(0.08, 0.08))) +
    facet_wrap(~ Cluster, ncol = 1, scales = "free_x") +
    labs(title = NULL, x = NULL, y = NULL, color = NULL) +
    panel_theme(base_size) +
    theme(
      strip.background = element_rect(fill = "grey90", colour = "grey80", linewidth = 0.25),
      strip.text = element_text(
        size = 9.5, face = "plain", colour = "black",
        margin = margin(2, 1, 2, 1)
      ),
      plot.title = element_blank(),
      axis.title = element_blank(),
      axis.text = element_text(size = 7, colour = "black"),
      axis.line = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      legend.text = element_text(size = base_size),
      plot.margin = margin(2, 4, 2, 4)
    )
}


## Short labels for enrichment terms (Reactome / GO / CORUM)
.enrich_term_short <- function(term_clean) {
  map <- c(
    "RUNX1 And FOXP3 Control Development Of Regulatory T Lymphocytes (Tregs)" =
      "RUNX1/FOXP3 Treg development",
    "RUNX1 Regulates Transcription Of Genes Involved In WNT Signaling" =
      "RUNX1-WNT signaling",
    "CD19-Vav-PI 3-kinase (p85 subunit) complex" =
      "CD19-Vav-PI3K(p85) complex",
    "LIFR-LIF-gp130 complex" =
      "LIFR-LIF-gp130 complex",
    "CD163 Mediating An Anti-Inflammatory Response" =
      "CD163 anti-inflammatory response",
    "Nef Mediated CD4 Down-regulation" =
      "Nef-mediated CD4 downregulation",
    "Binding And Entry Of HIV Virion" =
      "HIV virion binding & entry",
    "Early Phase Of HIV Life Cycle" =
      "Early HIV life cycle",
    "Scavenging Of Heme From Plasma" =
      "Plasma heme scavenging",
    "Alpha-defensins" =
      "Alpha-defensins",
    "negative regulation of regulatory T cell differentiation" =
      "Neg. regulation of Treg differentiation",
    "cellular response to granulocyte macrophage colony-stimulating factor stimulus" =
      "Response to GM-CSF",
    "epithelial tube branching involved in lung morphogenesis" =
      "Epithelial branching (lung)",
    "positive regulation of cell-cell adhesion mediated by integrin" =
      "Integrin-mediated cell-cell adhesion",
    "positive regulation of epithelial cell differentiation involved in kidney development" =
      "Epithelial differentiation (kidney)",
    "chronic inflammatory response" =
      "Chronic inflammatory response",
    "lung morphogenesis" =
      "Lung morphogenesis",
    "Translocation Of ZAP-70 To Immunological Synapse" =
      "ZAP-70 to immunological synapse",
    "IL-6-type Cytokine Receptor Ligand Interactions" =
      "IL-6-type cytokine-receptor interactions",
    "Signal Transduction" =
      "Signal transduction",
    "regulation of transcription, DNA-templated" =
      "DNA transcription",
    "regulation of transcription by RNA polymerase II" =
      "Pol II transcription",
    "regulation of gene expression" =
      "Gene-expression regulation"
  )
  out <- unname(map[term_clean])
  ifelse(is.na(out), term_clean, out)
}

plot_cn_enrichment_bar <- function(cluster_enrich, target_CN, n_top = 5L,
                                   label_inside = TRUE) {
  plot_dat <- cluster_enrich %>%
    filter(Cluster == target_CN) %>%
    distinct(Term, .keep_all = TRUE) %>%
    mutate(
      log10FDR = -log10(pmax(FDR, 1e-300)),
      Term_clean = str_replace(Term, "\\s+R-HSA-[0-9]+$", ""),
      Term_clean = str_replace(Term_clean, "\\s*\\(GO:[0-9]+\\)$", ""),
      Term_clean = str_replace(Term_clean, "\\s*\\(human\\)$", ""),
      Term_short = .enrich_term_short(Term_clean)
    )
  selected <- plot_dat %>%
    filter(Odds.Ratio > 1, is.finite(FDR), is.finite(Odds.Ratio)) %>%
    arrange(FDR, desc(Odds.Ratio)) %>%
    head(n_top)
  if (nrow(selected) == 0) return(NULL)
  ## Enrichment score: log2(odds ratio); bars ordered long to short.
  selected <- selected %>%
    mutate(
      enrich_score = log2(pmax(Odds.Ratio, 1e-6)),
      Term_ord = factor(Term_short, levels = Term_short[order(enrich_score, log10FDR)])
    )
  fill_fun <- scales::col_numeric(
    palette = c("#FEE0D2", "#FC9272", "#EF3B2C", "#A50F15"),
    domain = range(selected$enrich_score, na.rm = TRUE)
  )
  selected <- selected %>%
    mutate(
      fill_col = fill_fun(enrich_score),
      fill_lum = vapply(fill_col, function(cl) {
        rgb <- as.numeric(grDevices::col2rgb(cl))
        0.299 * rgb[1] + 0.587 * rgb[2] + 0.114 * rgb[3]
      }, numeric(1)),
      txt_col = ifelse(fill_lum < 140, "white", "grey10"),
      y_lab = 0.04 * max(enrich_score)
    )
  y_max <- max(selected$enrich_score, na.rm = TRUE)
  p <- ggplot(selected, aes(Term_ord, enrich_score, fill = enrich_score)) +
    geom_col(width = 0.72, color = NA) +
    coord_flip(clip = "off") +
    scale_fill_gradientn(
      colours = c("#FEE0D2", "#FC9272", "#EF3B2C", "#A50F15"),
      name = expression(log[2]~OR)
    ) +
    labs(title = target_CN, x = NULL, y = expression(log[2](OR))) +
    panel_theme(11) +
    theme(
      plot.title = element_text(size = 11, hjust = 0.5),
      axis.text.x = element_text(size = 9),
      axis.title.x = element_text(size = 10),
      legend.position = "right",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      plot.margin = margin(4, 10, 4, 6),
      panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3)
    ) +
    expand_limits(y = y_max * 1.08)
  if (isTRUE(label_inside)) {
    p <- p +
      geom_text(
        aes(y = y_lab, label = Term_short, color = txt_col),
        hjust = 0, size = 3.0, fontface = "plain", show.legend = FALSE
      ) +
      scale_color_identity() +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  } else {
    p <- p + theme(axis.text.y = element_text(size = 8.5, color = "grey15"))
  }
  p
}

## ---- load ----
stopifnot(file.exists(results_src))
clustered_df <- read.csv(results_src, check.names = FALSE)
if (names(clustered_df)[1] %in% c("", "X")) clustered_df <- clustered_df[, -1, drop = FALSE]
clustered_df$Cell_ID <- as.character(clustered_df$Cell_ID)
clustered_df$Neighborhood_Cluster <- as.character(clustered_df$Neighborhood_Cluster)
clustered_df$celltype <- as.character(clustered_df$celltype)
write.csv(clustered_df, results_csv, row.names = FALSE)

## Spatial edges: load and apply biological filter
filter_edges <- function(edges, max_edge_length = NULL, max_degree = 10L,
                             require_mutual = TRUE) {
  edges <- as.data.frame(edges)
  edges$from <- as.character(edges$from)
  edges$to <- as.character(edges$to)
  if (is.null(max_edge_length) || !is.finite(max_edge_length)) {
    max_edge_length <- as.numeric(stats::quantile(edges$dist, 0.95, na.rm = TRUE))
  }
  out <- .filter_edges_biological(
    data.table::as.data.table(edges),
    max_edge_length = max_edge_length,
    max_degree = max_degree,
    require_mutual = require_mutual,
    verbose = TRUE
  )
  as.data.frame(out)
}

if (file.exists(edges_src)) {
  spatial_edges <- readRDS(edges_src)
  spatial_edges <- filter_edges(spatial_edges)
} else if (file.exists(edges_rds)) {
  spatial_edges <- filter_edges(readRDS(edges_rds))
} else {
  df_net <- prepare_data(
    data.table::as.data.table(clustered_df),
    cell_id_col = "Cell_ID", x_col = "X", y_col = "Y", celltype_col = "celltype"
  )
  spatial_edges <- build_spatial_network(
    df_net, method = "auto", celltype_col = "celltype",
    n_neighbors = 10, require_mutual = TRUE, verbose = TRUE
  )
  spatial_edges <- as.data.frame(spatial_edges)
}
saveRDS(spatial_edges, edges_rds)

## Recompute interactions on filtered edges
interaction_results <- analyze_spatial_interactions(
  clustered_df, spatial_edges, celltype_col = "celltype", normalize_abundance = TRUE
)
saveRDS(interaction_results, interaction_rds)
if (file.exists(interaction_src)) {
  file.copy(interaction_src, file.path(out_dir, "tsu33_interaction_results_raw.rds"), overwrite = TRUE)
}

## Protein matrix (optional; DE prefers precomputed diff)
protein <- NULL
if (file.exists(obj_rds)) {
  codex.obj <- readRDS(obj_rds)
  assay_use <- if ("Akoya" %in% Assays(codex.obj)) "Akoya" else DefaultAssay(codex.obj)
  protein <- tryCatch(
    as.data.frame(
      t(as.matrix(GetAssayData(codex.obj, assay = assay_use, layer = "data"))),
      check.names = FALSE
    ),
    error = function(e) NULL
  )
}
if ((is.null(protein) || ncol(protein) == 0) && file.exists(expr_csv)) {
  protein <- read.csv(expr_csv, check.names = FALSE, row.names = 1)
}

ct_map <- panel_cols(clustered_df$celltype)

## ---- a: cell-type spatial ----
save_gg(
  visualize_spatial_distribution(
    clustered_df, celltype_col = "celltype", point_size = 0.35, point_alpha = 0.9,
    color_palette = ct_map, title = NULL
  ) + panel_theme(12),
  "Fig4a_spatial_distribution_celltype", width = 5.8, height = 4.8, legend = "Fig4_celltype"
)

## ---- b: CN composition heatmap ----
save_hm(
  draw_cn_heatmap(clustered_df, highlight_cn = target_cn),
  "Fig4b_CN_composition_heatmap", width = 8.5, height = 7.5, legend = "Fig4b_Zscore"
)

## ---- d: within-CN heterotypic contact enrichment ----
## Edges with both ends in the same target CN; enrichment = obs/expected
## under random labeling within that CN subgraph.
.ct_pair_short <- function(x) {
  x <- as.character(x)
  map <- c(
    "EMT_like_Tumor" = "EMT-Tum",
    "Proliferating_Tumor" = "Prolif-Tum",
    "Proliferating_DC" = "Prolif-DC",
    "PD-1_CD4_T" = "PD1-CD4T",
    "M2_Macrophage" = "M2-Mac",
    "Endothelial" = "Endo",
    "Tumor" = "Tumor",
    "CAF" = "CAF",
    "DC" = "DC",
    "B" = "B",
    "GC_B" = "GC-B"
  )
  ifelse(x %in% names(map), unname(map[x]), x)
}

within_cn_contact_tables <- function(df, edges, cns, min_obs = 10L) {
  df <- df %>% mutate(Cell_ID = as.character(Cell_ID))
  edges <- edges %>%
    mutate(from = as.character(from), to = as.character(to))
  rows <- lapply(cns, function(cn) {
    ids <- df$Cell_ID[as.character(df$Neighborhood_Cluster) == cn]
    e_in <- edges %>% filter(from %in% ids, to %in% ids)
    if (nrow(e_in) < 20L) return(NULL)
    sub <- df %>% filter(Cell_ID %in% ids)
    ir <- analyze_spatial_interactions(
      sub, e_in, celltype_col = "celltype", normalize_abundance = TRUE
    )
    en <- ir$enrichment_matrix
    obs <- ir$interaction_matrix
    idx <- which(upper.tri(en, diag = FALSE), arr.ind = TRUE)
    data.frame(
      CN = cn,
      type_a = rownames(en)[idx[, 1]],
      type_b = colnames(en)[idx[, 2]],
      enrichment = as.numeric(en[idx]),
      n_contacts = as.numeric(obs[idx]),
      stringsAsFactors = FALSE
    ) %>%
      filter(is.finite(enrichment), n_contacts >= min_obs)
  })
  bind_rows(rows)
}

contact_all_df <- within_cn_contact_tables(
  clustered_df, spatial_edges, target_cn, min_obs = 10L
)
contact_pairs_df <- contact_all_df %>%
  group_by(CN) %>%
  group_modify(~ {
    enriched <- .x %>% filter(enrichment > 1) %>% arrange(desc(enrichment))
    pick <- if (nrow(enriched) >= 2L) enriched else arrange(.x, desc(enrichment))
    head(pick, 4L)
  }) %>%
  ungroup()
write.csv(
  contact_pairs_df,
  file.path(out_dir, "Fig4d_within_CN_contact_pairs.csv"),
  row.names = FALSE
)

pair_key <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "||")
top_keys <- unique(pair_key(contact_pairs_df$type_a, contact_pairs_df$type_b))
cn_cols <- c(CN6 = "#E64B35", CN10 = "#4DBBD5", CN27 = "#00A087")

plot_dat_d <- contact_all_df %>%
  mutate(pkey = pair_key(type_a, type_b)) %>%
  filter(pkey %in% top_keys, enrichment > 1) %>%
  mutate(
    pair = paste(.ct_pair_short(type_a), .ct_pair_short(type_b), sep = "–"),
    log2OR = log2(pmax(enrichment, 1e-6)),
    CN = factor(CN, levels = target_cn)
  )
pair_ord <- plot_dat_d %>%
  group_by(pair) %>%
  summarise(m = max(log2OR), .groups = "drop") %>%
  arrange(m) %>%
  pull(pair)
plot_dat_d <- plot_dat_d %>%
  mutate(pair = factor(pair, levels = pair_ord))

p_d_contact <- ggplot(plot_dat_d, aes(log2OR, pair, fill = CN)) +
  geom_col(position = position_dodge2(width = 0.85, preserve = "single"), width = 0.78, color = NA) +
  scale_fill_manual(values = cn_cols, name = NULL) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.06)),
    breaks = scales::pretty_breaks(n = 4)
  ) +
  labs(x = expression(log[2](OR)), y = NULL) +
  panel_theme(8) +
  theme(
    axis.text.y = element_text(size = 6.5),
    axis.text.x = element_text(size = 6.5),
    axis.title.x = element_text(size = 7),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.justification = "center",
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.32, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
    plot.margin = margin(2, 4, 2, 2)
  ) +
  guides(fill = guide_legend(nrow = 1))

save_gg(
  p_d_contact + theme(legend.position = "none"),
  "Fig4d_within_CN_contact",
  width = 8, height = 5.5, units = "cm",
  legend = NULL
)

## Vertical CN color legend for panel d
tryCatch({
  p_leg <- ggplot(
    data.frame(CN = factor(target_cn, levels = target_cn), x = 1),
    aes(x, x, fill = CN)
  ) +
    geom_col() +
    scale_fill_manual(values = cn_cols, name = NULL, drop = FALSE) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.direction = "vertical",
      legend.background = element_blank(),
      legend.key = element_rect(fill = NA, colour = NA),
      legend.text = element_text(size = 10, color = "black"),
      legend.key.size = unit(0.45, "cm"),
      plot.background = element_blank()
    ) +
    guides(fill = guide_legend(ncol = 1))
  boxes <- cowplot::get_plot_component(p_leg, "guide-box", return_all = TRUE)
  sizes <- vapply(boxes, function(g) {
    if (is.null(g) || inherits(g, "zeroGrob")) return(0)
    sum(as.numeric(grid::convertWidth(g$widths, "mm")), na.rm = TRUE) *
      sum(as.numeric(grid::convertHeight(g$heights, "mm")), na.rm = TRUE)
  }, numeric(1))
  leg <- boxes[[which.max(sizes)]]
  .save_legend_grob(
    leg,
    file.path(out_dir, "Fig4d_contact_OR_legend.pdf"),
    file.path(out_dir, "Fig4d_contact_OR_legend.png"),
    width = 2.8, height = 3.4, units = "cm"
  )
}, error = function(e) invisible(NULL))

## Merge CN into Seurat object when available
if (exists("codex.obj") && inherits(codex.obj, "Seurat")) {
  nc_map <- setNames(
    as.character(clustered_df$Neighborhood_Cluster),
    as.character(clustered_df$Cell_ID)
  )
  codex.obj$Cell_ID <- as.character(colnames(codex.obj))
  codex.obj$Neighborhood_Cluster <- nc_map[codex.obj$Cell_ID]
  save(codex.obj, clustered_df, file = file.path(out_dir, "tsu33_cn.Rdata"))
  if (!file.exists(cn_rdata)) save(codex.obj, clustered_df, file = cn_rdata)
} else {
  save(clustered_df, file = file.path(out_dir, "tsu33_cn.Rdata"))
}

## ---- e: volcano (stacked + per-CN) ----
diff_csv <- file.path(out_dir, "diff_results.csv")
if (file.exists(diff_src)) {
  diff_results <- read.csv(diff_src, check.names = FALSE)
  if (names(diff_results)[1] %in% c("", "X")) diff_results <- diff_results[, -1, drop = FALSE]
  write.csv(diff_results, diff_csv, row.names = FALSE)
} else if (file.exists(diff_csv)) {
  diff_results <- read.csv(diff_csv, check.names = FALSE)
} else if (!is.null(protein) && ncol(protein) > 0) {
  protein_df <- prepare_protein_data(clustered_df, protein)
  diff_results <- perform_differential_expression(protein_df)
  write.csv(diff_results, diff_csv, row.names = FALSE)
} else {
  stop("No differential results available")
}
diff_sub <- diff_results %>% filter(Cluster %in% target_cn)
if (!"MeanDiff" %in% names(diff_sub) && all(c("Mean_Target", "Mean_Control") %in% names(diff_sub))) {
  diff_sub$MeanDiff <- diff_sub$Mean_Target - diff_sub$Mean_Control
}
if (!"MeanDiff" %in% names(diff_sub) && "Log2FC" %in% names(diff_sub)) {
  diff_sub$MeanDiff <- diff_sub$Log2FC
}
if (nrow(diff_sub) > 0) {
  p_vol <- plot_volcano_vertical(
    diff_sub, cluster_levels = target_cn
  )
  save_gg(
    p_vol,
    "Fig4e_volcano_CN6_CN10_CN27",
    width = 2.2, height = 5.7,
    legend = "Fig4e_volcano", legend_width = 1.2, legend_height = 1.8
  )
  for (cn in target_cn) {
    diff_one <- diff_sub %>% filter(Cluster == cn)
    if (nrow(diff_one) == 0) next
    p_one <- plot_volcano_vertical(
      diff_one, cluster_levels = cn
    )
    save_gg(
      p_one,
      paste0("Fig4e_volcano_", cn),
      width = 2.2, height = 2.0,
      keep_title = TRUE,
      legend = NULL
    )
  }
}

## ---- f: enrichment bars ----
cluster_enrich <- NULL
if (file.exists(enrich_src)) {
  cluster_enrich <- readRDS(enrich_src)
  saveRDS(cluster_enrich, file.path(out_dir, "cluster_enrich.rds"))
} else {
  protein_mapping <- data.frame(
    Original_Name = c(
      "PD-1", "CD19", "CD8", "CD45RO", "CD68", "Pan-Cytokeratin", "Podoplanin", "CD4", "TP53",
      "CD31", "CD11c", "CD163", "Galectin3", "CXCL13", "PD-L1", "CTLA4", "Ki67", "NKX2-1",
      "FOXP3", "a-SMA", "CD56", "LIF", "b-Catenin1", "Caveolin"
    ),
    UniProt_ID = c(
      "Q15116", "P15391", "P01732", "P08575", "P34810", "P08727", "Q86YL7", "P01730", "P04637",
      "P16284", "P20702", "Q86VB7", "P17931", "O43927", "Q9NZQ7", "P16410", "P46013", "P43699",
      "Q9BZS1", "P62736", "P13591", "P15018", "P35222", "Q03135"
    ),
    Gene_Symbol = c(
      "PDCD1", "CD19", "CD8A", "PTPRC", "CD68", "KRT19", "PDPN", "CD4", "TP53",
      "PECAM1", "ITGAX", "CD163", "LGALS3", "CXCL13", "CD274", "CTLA4", "MKI67", "NKX2-1",
      "FOXP3", "ACTA2", "NCAM1", "LIF", "CTNNB1", "CAV1"
    )
  )
  cluster_enrich <- tryCatch({
    perform_cluster_enrichment(
      diff_results = diff_results, species = "human",
      protein_databases = c("GO_Biological_Process_2023", "Reactome_2022", "CORUM"),
      pvalueCutoff = 0.1, protein_mapping = protein_mapping
    )
  }, error = function(e) {
    NULL
  })
  if (!is.null(cluster_enrich) && is.data.frame(cluster_enrich) && nrow(cluster_enrich) > 0) {
    saveRDS(cluster_enrich, file.path(out_dir, "cluster_enrich.rds"))
  }
}

enrich_plots <- list()
if (!is.null(cluster_enrich) && is.data.frame(cluster_enrich) && nrow(cluster_enrich) > 0 &&
    all(c("Cluster", "FDR", "Term", "Odds.Ratio") %in% names(cluster_enrich))) {
  first_leg <- TRUE
  for (cn in target_cn) {
    p_en <- tryCatch(
      plot_cn_enrichment_bar(cluster_enrich, cn, n_top = 5L, label_inside = TRUE),
      error = function(e) {
        NULL
      }
    )
    if (is.null(p_en)) {
      next
    }
    enrich_plots[[cn]] <- p_en
    save_gg(
      p_en,
      paste0("Fig4f_enrichment_", cn),
      width = 3.6, height = 2.6,
      keep_title = TRUE,
      legend = if (isTRUE(first_leg)) "Fig4f_enrichment" else NULL,
      legend_width = 1.8, legend_height = 3.0
    )
    first_leg <- FALSE
  }
  ## Stack order: CN6, CN10, CN27
  if (length(enrich_plots) > 0) {
    en_ord <- enrich_plots[intersect(target_cn, names(enrich_plots))]
    p_f_col <- cowplot::plot_grid(
      plotlist = lapply(unname(en_ord), function(p) p + theme(legend.position = "none")),
      ncol = 1, align = "none", rel_heights = rep(1, length(en_ord))
    )
    save_gg(p_f_col, "Fig4f_enrichment_stack", width = 3.6, height = 7.5, legend = NULL)
  }
}

## ---- g: KM survival ----
fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-4) return(format(p, scientific = TRUE, digits = 3))
  sprintf("%.4f", p)
}
km_plots <- list()
if (file.exists(surv_csv)) {
  suppressPackageStartupMessages({ library(survival); library(survminer) })
  surv_data <- read.csv(surv_csv, stringsAsFactors = FALSE)
  available <- target_cn[target_cn %in% colnames(surv_data)]
  res_list <- list()
  for (cl in available) {
    dd <- data.frame(
      Patient_ID = surv_data$Patient_ID,
      OS.time = as.numeric(surv_data$OS.time),
      OS = as.integer(surv_data$OS),
      Signature = as.numeric(surv_data[[cl]]),
      stringsAsFactors = FALSE
    )
    dd <- dd[stats::complete.cases(dd) & dd$OS.time > 0, ]
    cut <- survminer::surv_cutpoint(
      dd, time = "OS.time", event = "OS", variables = "Signature", minprop = 0.1
    )
    opt_cut <- as.numeric(cut$cutpoint[1, "cutpoint"])
    dd$Group_optimal <- factor(
      ifelse(dd$Signature > opt_cut, "High", "Low"), levels = c("High", "Low")
    )
    sdif <- survival::survdiff(Surv(OS.time, OS) ~ Group_optimal, data = dd)
    p_raw <- 1 - stats::pchisq(sdif$chisq, length(sdif$n) - 1)
    res_list[[cl]] <- list(data = dd, cutoff = opt_cut, p_raw = p_raw)
  }
  p_bh <- p.adjust(vapply(available, function(cl) res_list[[cl]]$p_raw, numeric(1)), method = "BH")
  names(p_bh) <- available
  write.csv(
    data.frame(
      Cluster = available,
      Cutoff = vapply(available, function(cl) res_list[[cl]]$cutoff, numeric(1)),
      P_Value = vapply(available, function(cl) res_list[[cl]]$p_raw, numeric(1)),
      P_BH = as.numeric(p_bh),
      stringsAsFactors = FALSE
    ),
    file.path(out_dir, "Fig4g_survival_stats.csv"), row.names = FALSE
  )
  for (cl in available) {
    dd <- res_list[[cl]]$data
    fit <- survfit(Surv(OS.time, OS) ~ Group_optimal, data = dd)
    p_lab <- paste0("p_adj = ", fmt_p(p_bh[[cl]]))
    ## KM curves with BH-adjusted p-value annotation
    g_orig <- ggsurvplot(
      fit, data = dd, pval = FALSE, risk.table = FALSE,
      palette = c("#D55E00", "#0072B2"),
      legend.labs = c("High", "Low"),
      title = cl, xlab = "Time (days)", ylab = "Survival Probability",
      surv.median.line = "hv", ggtheme = panel_theme(12),
      font.title = c(18, "plain", "black"),
      font.x = c(13, "plain", "black"),
      font.y = c(13, "plain", "black"),
      font.tickslab = c(11, "plain", "black"),
      font.legend = c(12, "plain", "black")
    )
    g_orig$plot <- g_orig$plot +
      annotate(
        "text",
        x = max(dd$OS.time, na.rm = TRUE),
        y = 0.98,
        label = p_lab,
        hjust = 1.05, vjust = 1, size = 5.2
      ) +
      theme(plot.title = element_text(size = 18, hjust = 0.5, face = "plain"))
    save_gg(
      g_orig$plot, paste0("Fig4g_KM_", cl),
      width = 4.8, height = 4.0, keep_title = TRUE,
      legend = if (cl == available[1]) "Fig4g_KM" else NULL
    )
    ## KM panels for vertical stack
    g <- ggsurvplot(
      fit, data = dd, pval = FALSE, risk.table = FALSE,
      palette = c("#D55E00", "#0072B2"),
      legend.labs = c("High", "Low"),
      title = cl, xlab = "Time (days)", ylab = "Survival Prob.",
      surv.median.line = "hv", ggtheme = panel_theme(10),
      font.title = c(15, "plain", "black"),
      font.x = c(10, "plain", "black"),
      font.y = c(10, "plain", "black"),
      font.tickslab = c(9, "plain", "black"),
      font.legend = c(10, "plain", "black")
    )
    g$plot <- g$plot +
      annotate(
        "text",
        x = max(dd$OS.time, na.rm = TRUE),
        y = 0.98,
        label = p_lab,
        hjust = 1.05, vjust = 1, size = 4.0
      ) +
      theme(
        plot.title = element_text(size = 15, hjust = 0.5, face = "plain"),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 9),
        legend.text = element_text(size = 10),
        plot.margin = margin(2, 4, 2, 2)
      )
    km_plots[[cl]] <- g$plot
  }
  if (length(km_plots) > 0) {
    km_ord <- km_plots[intersect(target_cn, names(km_plots))]
    p_g_col <- cowplot::plot_grid(
      plotlist = unname(km_ord), ncol = 1, align = "none",
      rel_heights = rep(1, length(km_ord))
    )
    save_gg(p_g_col, "Fig4g_KM_stack", width = 3.2, height = 7.5, legend = NULL)
  }
}
