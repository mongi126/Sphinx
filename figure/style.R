## Shared figure style for Fig2–Fig5
## Soft pastel categorical palette (ColorBrewer Set2/Paired family).
## Prefer mid-light tones so large tissue domains do not read as dark blobs.
## Heatmaps keep classic red–blue diverging scale.

## light-first order: when assigned by frequency, abundant types get pale colors
QUAL_PALETTE <- c(
  "#FFED6F", # butter
  "#CCEBC5", # pale mint
  "#A6CEE3", # baby blue
  "#FCCDE5", # soft pink
  "#FDB462", # apricot
  "#B3DE69", # soft lime
  "#BEBADA", # lavender
  "#8DD3C7", # mint teal
  "#FB9A99", # light rose
  "#80B1D3", # soft sky
  "#CAB2D6", # lilac
  "#FB8072", # soft coral
  "#B2DF8A", # light green
  "#FDBF6F", # peach
  "#BC80BD", # orchid
  "#D9D9D9"  # soft gray
)

## optional semantic accents for common lineages (still soft)
CT_HINT <- c(
  Epithelial = "#FDB462",
  Tumor = "#FB8072",
  Proliferating_Tumor = "#FDBF6F",
  EMT_like_Tumor = "#BC80BD",
  Fibroblast = "#8DD3C7",
  CAF = "#66C2A5",
  Endothelial = "#80B1D3",
  Macrophage = "#BEBADA",
  M2_Macrophage = "#BEBADA",
  Monocyte = "#CAB2D6",
  DC = "#B3DE69",
  Proliferating_DC = "#CCEBC5",
  B = "#FCCDE5",
  GC_B = "#FBB4AE",
  Plasma = "#F1B6DA",
  T = "#A6CEE3",
  CD4T = "#A6CEE3",
  `PD-1_CD4_T` = "#92C5DE",
  CD8T = "#4393C3",
  Treg = "#D1E5F0",
  NK = "#FFFFBF",
  Mast = "#FEE08B",
  Neutrophil = "#FDAE61",
  SMC = "#C7EAE5"
)

## classic red–blue for heatmaps / continuous diverging
HM_COLS <- c("#2166AC", "#F7F7F7", "#B2182B")
HM_COL_FUN <- function() circlize::colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))

## Assign colors to categories.
## order = "freq": most abundant -> lightest pastel (best for spatial maps)
## semantic = TRUE: use CT_HINT when names match, fill rest from QUAL_PALETTE
panel_cols <- function(x, order = c("freq", "alpha", "as_is"), semantic = TRUE) {
  order <- match.arg(order)
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return(character(0))
  tab <- sort(table(x), decreasing = TRUE)
  if (order == "freq") {
    lv <- names(tab)
  } else if (order == "alpha") {
    lv <- sort(unique(x))
  } else {
    lv <- unique(x)
  }
  n <- length(lv)
  base <- if (n <= length(QUAL_PALETTE)) {
    QUAL_PALETTE[seq_len(n)]
  } else {
    grDevices::colorRampPalette(QUAL_PALETTE)(n)
  }
  cols <- setNames(base, lv)
  if (isTRUE(semantic)) {
    hit <- intersect(lv, names(CT_HINT))
    if (length(hit)) cols[hit] <- CT_HINT[hit]
  }
  cols
}

## Larger base fonts for assembly-friendly panels
panel_theme <- function(base = 12) {
  ggplot2::theme_classic(base_size = base) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = "black", size = base, face = "plain"),
      plot.title = ggplot2::element_text(hjust = 0.5, size = base + 1, face = "plain"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = base, face = "plain"),
      axis.title = ggplot2::element_text(size = base, face = "plain"),
      axis.text = ggplot2::element_text(size = max(10, base - 1), colour = "black", face = "plain"),
      legend.title = ggplot2::element_text(size = max(10, base - 1), face = "plain"),
      legend.text = ggplot2::element_text(size = max(10, base - 1), face = "plain"),
      legend.key.size = grid::unit(0.45, "cm"),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.55),
      strip.background = ggplot2::element_rect(fill = "grey95", colour = NA),
      strip.text = ggplot2::element_text(size = max(10, base - 1), face = "plain")
    )
}

## Pretty stacked CN composition bars
plot_cn_composition_pretty <- function(composition_df,
                                       color_map,
                                       cluster_col = "Neighborhood_Cluster",
                                       celltype_col = "celltype",
                                       value_col = "proportion",
                                       title = NULL) {
  dd <- composition_df
  cn <- as.character(dd[[cluster_col]])
  cn_u <- unique(cn)
  cn_num <- suppressWarnings(as.integer(sub("^CN", "", cn_u)))
  cn_ord <- cn_u[order(cn_num, cn_u)]
  ct_ord <- dd %>%
    dplyr::group_by(.data[[celltype_col]]) %>%
    dplyr::summarise(tot = sum(.data[[value_col]], na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(.data$tot)) %>%
    dplyr::pull(1)
  dd[[cluster_col]] <- factor(cn, levels = cn_ord)
  dd[[celltype_col]] <- factor(as.character(dd[[celltype_col]]), levels = rev(ct_ord))
  fill_vals <- color_map[as.character(ct_ord)]
  fill_vals <- fill_vals[!is.na(names(fill_vals))]
  ggplot2::ggplot(
    dd,
    ggplot2::aes(
      x = .data[[cluster_col]], y = .data[[value_col]],
      fill = .data[[celltype_col]]
    )
  ) +
    ggplot2::geom_col(width = 0.78, colour = "white", linewidth = 0.25) +
    ggplot2::scale_fill_manual(values = fill_vals, breaks = ct_ord, name = "Cell type") +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      expand = c(0, 0), limits = c(0, 1.002)
    ) +
    ggplot2::labs(x = NULL, y = "Proportion", title = title, fill = "Cell type") +
    panel_theme(12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, size = 11),
      panel.grid.major.y = ggplot2::element_line(colour = "grey92", linewidth = 0.35),
      plot.margin = ggplot2::margin(8, 10, 6, 8)
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE, ncol = 1))
}

## Save ggplot: panel without legend; optional separate legend file.
## legend: NULL/FALSE | TRUE | "id"
## units: "in" or "cm"
## PNG uses transparent device background; legends have no panel box/grid.
.save_legend_grob <- function(leg, path_pdf, path_png, width, height, units = "in") {
  p_leg <- cowplot::ggdraw(leg) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "transparent", colour = NA),
      panel.background = ggplot2::element_rect(fill = "transparent", colour = NA),
      panel.border = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(0, 0, 0, 0)
    )
  ggplot2::ggsave(
    path_pdf, p_leg, width = width, height = height, units = units,
    device = grDevices::cairo_pdf, bg = "transparent"
  )
  ggplot2::ggsave(
    path_png, p_leg, width = width, height = height, units = units,
    dpi = 300, bg = "transparent"
  )
}

save_gg_dir <- function(out_dir) {
  force(out_dir)
  function(plot_obj, name, width = 5.5, height = 4.5,
           legend_width = 2.6, legend_height = NULL,
           legend = NULL, keep_title = FALSE, units = "in") {
    if (is.null(legend_height)) legend_height <- min(height, if (identical(units, "cm")) 12 else 5.0)
    p <- plot_obj
    if (!isTRUE(keep_title)) {
      p <- p + ggplot2::labs(title = NULL) + ggplot2::theme(plot.title = ggplot2::element_blank())
    }
    export_leg <- FALSE
    leg_stem <- NULL
    if (isTRUE(legend)) {
      export_leg <- TRUE
      leg_stem <- paste0(name, "_legend")
    } else if (is.character(legend) && nzchar(legend[1])) {
      export_leg <- TRUE
      leg_stem <- if (grepl("_legend$", legend[1])) legend[1] else paste0(legend[1], "_legend")
    }
    p_main <- p + ggplot2::theme(
      legend.position = "none",
      plot.background = ggplot2::element_rect(fill = "transparent", colour = NA)
    )
    ggplot2::ggsave(
      file.path(out_dir, paste0(name, ".pdf")), p_main,
      width = width, height = height, units = units, device = grDevices::cairo_pdf, bg = "transparent"
    )
    ggplot2::ggsave(
      file.path(out_dir, paste0(name, ".png")), p_main,
      width = width, height = height, units = units, dpi = 300, bg = "transparent"
    )
    if (isTRUE(export_leg)) {
      p_for_leg <- p + ggplot2::theme(
        legend.position = "right",
        legend.background = ggplot2::element_rect(fill = "transparent", colour = NA),
        legend.box.background = ggplot2::element_rect(fill = "transparent", colour = NA),
        legend.key = ggplot2::element_rect(fill = "transparent", colour = NA),
        panel.background = ggplot2::element_rect(fill = "transparent", colour = NA),
        panel.border = ggplot2::element_blank(),
        panel.grid = ggplot2::element_blank(),
        plot.background = ggplot2::element_rect(fill = "transparent", colour = NA)
      )
      leg <- tryCatch({
        boxes <- cowplot::get_plot_component(p_for_leg, "guide-box", return_all = TRUE)
        if (is.null(boxes) || !length(boxes)) return(NULL)
        sizes <- vapply(boxes, function(g) {
          if (is.null(g) || inherits(g, "zeroGrob")) return(0)
          sum(as.numeric(grid::convertWidth(g$widths, "mm")), na.rm = TRUE) *
            sum(as.numeric(grid::convertHeight(g$heights, "mm")), na.rm = TRUE)
        }, numeric(1))
        boxes[[which.max(sizes)]]
      }, error = function(e) tryCatch(cowplot::get_legend(p_for_leg), error = function(e2) NULL))
      if (!is.null(leg)) {
        .save_legend_grob(
          leg,
          file.path(out_dir, paste0(leg_stem, ".pdf")),
          file.path(out_dir, paste0(leg_stem, ".png")),
          width = legend_width, height = legend_height, units = units
        )
      }
    }
    invisible(NULL)
  }
}

## ComplexHeatmap: body without legend; optional separate transparent legend
## legend: NULL/FALSE | TRUE | "id"
save_hm_dir <- function(out_dir) {
  force(out_dir)
  function(obj, name, width = 6, height = 5, padding = NULL,
           legend_width = 2.2, legend_height = 4.0, legend = TRUE, units = "in") {
    w_in <- if (identical(units, "cm")) width / 2.54 else width
    h_in <- if (identical(units, "cm")) height / 2.54 else height
    lw_in <- if (identical(units, "cm")) legend_width / 2.54 else legend_width
    lh_in <- if (identical(units, "cm")) legend_height / 2.54 else legend_height
    pdf_f <- file.path(out_dir, paste0(name, ".pdf"))
    png_f <- file.path(out_dir, paste0(name, ".png"))
    is_ht <- inherits(obj, "Heatmap") || inherits(obj, "HeatmapList")
    if (is_ht) {
      draw_body <- function() {
        args <- list(object = obj, show_heatmap_legend = FALSE, show_annotation_legend = FALSE)
        if (!is.null(padding)) args$padding <- padding
        do.call(ComplexHeatmap::draw, args)
      }
      grDevices::cairo_pdf(pdf_f, width = w_in, height = h_in, bg = "transparent")
      draw_body(); grDevices::dev.off()
      grDevices::png(png_f, width = w_in, height = h_in, units = "in", res = 300, type = "cairo", bg = "transparent")
      draw_body(); grDevices::dev.off()
      export_leg <- FALSE
      leg_stem <- NULL
      if (isTRUE(legend)) {
        export_leg <- TRUE
        leg_stem <- paste0(name, "_legend")
      } else if (is.character(legend) && nzchar(legend[1])) {
        export_leg <- TRUE
        leg_stem <- if (grepl("_legend$", legend[1])) legend[1] else paste0(legend[1], "_legend")
      }
      if (isTRUE(export_leg)) {
        tryCatch({
          ht <- if (inherits(obj, "HeatmapList")) obj@ht_list[[1]] else obj
          lgd <- ComplexHeatmap::color_mapping_legend(ht@matrix_color_mapping, plot = FALSE)
          grDevices::cairo_pdf(
            file.path(out_dir, paste0(leg_stem, ".pdf")),
            width = lw_in, height = lh_in, bg = "transparent"
          )
          grid::grid.newpage()
          ComplexHeatmap::draw(lgd)
          grDevices::dev.off()
          grDevices::png(
            file.path(out_dir, paste0(leg_stem, ".png")),
            width = lw_in, height = lh_in, units = "in", res = 300, type = "cairo", bg = "transparent"
          )
          grid::grid.newpage()
          ComplexHeatmap::draw(lgd)
          grDevices::dev.off()
        }, error = function(e) invisible(NULL))
      }
    } else {
      draw_one <- function() {
        if (!is.null(obj$gtable)) {
          grid::grid.newpage(); grid::grid.draw(obj$gtable)
        } else {
          print(obj)
        }
      }
      grDevices::cairo_pdf(pdf_f, width = w_in, height = h_in, bg = "transparent")
      draw_one(); grDevices::dev.off()
      grDevices::png(png_f, width = w_in, height = h_in, units = "in", res = 300, type = "cairo", bg = "transparent")
      draw_one(); grDevices::dev.off()
    }
    invisible(NULL)
  }
}
