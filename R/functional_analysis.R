# R/functional_analysis.R

#' Prepare protein expression data for functional analysis
#'
#' @param cluster_df Data frame with neighborhood cluster assignments and spatial coordinates
#' @param expr_df Data frame with protein expression data (rows = cells, columns = proteins)
#' @return Merged data frame containing spatial, cluster, and protein expression data
#' @examples
#' cluster_df <- Sphinx:::.sphinx_example_df(10)
#' cluster_df$Neighborhood_Cluster <- sample(1:2, 10, replace = TRUE)
#' expr_df <- data.frame(PD.1 = rnorm(10), CD3 = rnorm(10),
#'   row.names = cluster_df$Cell_ID)
#' merged <- prepare_protein_data(cluster_df, expr_df)
#' nrow(merged)
#' @export
prepare_protein_data <- function(cluster_df, expr_df) {
  expr_df$Cell_ID <- rownames(expr_df)
  cluster_df$Cell_ID <- gsub("Cell_", "", cluster_df$Cell_ID)

  message("Cluster data ID sample: ", paste(head(cluster_df$Cell_ID), collapse = ", "))
  message("Expression data ID sample: ", paste(head(expr_df$Cell_ID), collapse = ", "))

  merged_df <- merge(
    cluster_df[, c("Cell_ID", "X", "Y", "Neighborhood_Cluster")],
    expr_df,
    by = "Cell_ID",
    all.x = TRUE
  )

  if ("X.y" %in% colnames(merged_df)) {
    merged_df$X.y <- NULL
  }
  if ("X.x" %in% colnames(merged_df)) {
    names(merged_df)[names(merged_df) == "X.x"] <- "X"
  }

  merged_df$Spatial_Zone <- as.character(merged_df$Neighborhood_Cluster)

  message("Merged dimensions: ", paste(dim(merged_df), collapse = " x "))
  message("NA values in PD.1: ", sum(is.na(merged_df$PD.1)), "/", nrow(merged_df))

  return(merged_df)
}

#' Default spatial block size for pseudo-replication-aware DE
#' @param x,y numeric coordinate vectors
#' @param target_cells_per_block approximate cells per block (default: 40)
#' @return numeric block side length
#' @keywords internal
.default_spatial_block_size <- function(x, y, target_cells_per_block = 40) {
  n <- length(x)
  xr <- diff(range(x, na.rm = TRUE))
  yr <- diff(range(y, na.rm = TRUE))
  if (!is.finite(xr) || !is.finite(yr) || xr <= 0 || yr <= 0 || n < 2) {
    return(1)
  }
  area <- xr * yr
  n_blocks <- max(4, ceiling(n / target_cells_per_block))
  # square tiling approximating n_blocks over the tissue extent
  side <- sqrt(area / n_blocks)
  max(side, min(xr, yr) / 50)
}

#' Perform differential expression analysis
#'
#' Compares each neighborhood cluster to all other cells. For CLR-normalized
#' protein values, reports mean difference (Target - Control) rather than
#' log2 fold change to avoid redundant log transforms.
#'
#' By default, tests are run on **spatial-block means** rather than individual
#' cells. Treating every cell as an independent observation overstates degrees
#' of freedom because nearby cells are spatially correlated (pseudo-replication).
#' Aggregating to tissue blocks reduces that inflation.
#'
#' @param protein_df Data frame with protein expression and cluster information
#' @param test_level `"spatial_block"` (default) or `"cell"` (legacy; not
#'   recommended due to pseudo-replication)
#' @param block_size numeric grid side length for spatial blocks. NULL = auto.
#' @param min_blocks minimum number of Target and Control blocks required
#'   (default: 3)
#' @return Data frame with columns Mean_Target, Mean_Control, MeanDiff, p.value,
#'   adj.p.value, Significance, and `n_target` / `n_control` (cells or blocks)
#' @examples
#' df <- Sphinx:::.sphinx_example_protein_df(40)
#' res <- perform_differential_expression(df)
#' head(res[, c("Protein", "Cluster", "MeanDiff")])
#' @export
perform_differential_expression <- function(protein_df,
                                            test_level = c("spatial_block", "cell"),
                                            block_size = NULL,
                                            min_blocks = 3L) {
  test_level <- match.arg(test_level)
  protein_df <- as.data.frame(protein_df)
  meta_cols <- c("Cell_ID", "X", "Y", "Neighborhood_Cluster", "Spatial_Zone",
                 "In_Cluster", "Spatial_Block")
  protein_cols <- setdiff(colnames(protein_df), meta_cols)

  if (!all(c("X", "Y") %in% names(protein_df)) && test_level == "spatial_block") {
    warning("X/Y missing; falling back to cell-level t-tests.")
    test_level <- "cell"
  }

  if (test_level == "spatial_block") {
    if (is.null(block_size)) {
      block_size <- .default_spatial_block_size(protein_df$X, protein_df$Y)
    }
    protein_df$Spatial_Block <- paste(
      floor(as.numeric(protein_df$X) / block_size),
      floor(as.numeric(protein_df$Y) / block_size),
      sep = "_"
    )
    message(
      "DE test_level=spatial_block | block_size=", round(block_size, 3),
      " | n_blocks=", length(unique(protein_df$Spatial_Block))
    )
  } else {
    warning(
      "test_level='cell' treats each cell as independent and can inflate ",
      "significance due to spatial pseudo-replication. Prefer 'spatial_block'."
    )
  }

  clusters <- unique(protein_df$Neighborhood_Cluster)
  results <- list()

  for (cluster in clusters) {
    protein_df$In_Cluster <- ifelse(
      protein_df$Neighborhood_Cluster == cluster,
      "Target",
      "Control"
    )

    cluster_results <- lapply(protein_cols, function(prot) {
      if (test_level == "spatial_block") {
        # Mean expression per spatial block within Target / Control
        blk <- stats::aggregate(
          protein_df[[prot]],
          by = list(
            block = protein_df$Spatial_Block,
            group = protein_df$In_Cluster
          ),
          FUN = function(z) mean(z, na.rm = TRUE)
        )
        names(blk)[3] <- "expr"
        target <- blk$expr[blk$group == "Target"]
        control <- blk$expr[blk$group == "Control"]
        target <- target[is.finite(target)]
        control <- control[is.finite(control)]
        min_n <- as.integer(min_blocks)
      } else {
        target <- protein_df[[prot]][protein_df$In_Cluster == "Target"]
        control <- protein_df[[prot]][protein_df$In_Cluster == "Control"]
        target <- target[is.finite(target)]
        control <- control[is.finite(control)]
        min_n <- 2L
      }

      if (length(target) < min_n || length(control) < min_n) {
        return(NULL)
      }

      test_res <- stats::t.test(target, control)
      # Effect size still reported on all cells in the cluster (not only blocks)
      mean_target <- mean(
        protein_df[[prot]][protein_df$In_Cluster == "Target"],
        na.rm = TRUE
      )
      mean_control <- mean(
        protein_df[[prot]][protein_df$In_Cluster == "Control"],
        na.rm = TRUE
      )
      data.frame(
        Protein = prot,
        Cluster = cluster,
        Mean_Target = mean_target,
        Mean_Control = mean_control,
        MeanDiff = mean_target - mean_control,
        p.value = test_res$p.value,
        n_target = length(target),
        n_control = length(control),
        test_level = test_level,
        stringsAsFactors = FALSE
      )
    })

    cluster_df <- do.call(rbind, cluster_results)
    if (!is.null(cluster_df)) {
      cluster_df$adj.p.value <- p.adjust(cluster_df$p.value, method = "BH")
      cluster_df$Significance <- ifelse(
        cluster_df$adj.p.value < 0.05,
        ifelse(cluster_df$MeanDiff > 0, "Up", "Down"),
        "NS"
      )
      results[[as.character(cluster)]] <- cluster_df
    }
  }

  final_result <- do.call(rbind, results)
  rownames(final_result) <- NULL
  return(final_result)
}

#' Plot volcano plots for differential proteins across all clusters
#'
#' Y-axis `-log10(adj.p)` is capped at `y_cap` (default 50). Points exceeding
#' the cap are drawn as upward triangles on the truncation line.
#'
#' @param diff_results Differential expression results
#' @param diff_thresh Mean difference threshold, Target - Control (default: 0.25)
#' @param p_thresh P-value threshold (default: 0.05)
#' @param y_cap Truncation for -log10(adj.p.value) (default: 50)
#' @param base_size Base font size (default: 14)
#' @param ncol Number of facet columns (default: 3)
#' @param save_plot Whether to save the plot (default: FALSE)
#' @param output_dir Output directory for saving (default: "plots")
#' @param filename Output filename (default: "volcano_plots")
#' @param width Plot width in inches (default: 12)
#' @param height Plot height in inches (default: 8)
#' @return ggplot object
#' @examples
#' \donttest{
#' df <- Sphinx:::.sphinx_example_protein_df(40)
#' res <- perform_differential_expression(df)
#' p <- plot_volcano_all_clusters(res)
#' class(p)
#' }
#' @export
plot_volcano_all_clusters <- function(diff_results, diff_thresh = 0.25, p_thresh = 0.05,
                                      y_cap = 50, base_size = 14, ncol = 3,
                                      save_plot = FALSE, output_dir = "plots",
                                      filename = "volcano_plots",
                                      width = 12, height = 8) {
  base_size <- max(8, as.numeric(base_size))
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Please install ggplot2 package")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Please install dplyr package")
  }
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Please install ggrepel package")
  }

  df <- diff_results %>%
    dplyr::mutate(
      Significance = dplyr::case_when(
        adj.p.value < p_thresh & MeanDiff >= diff_thresh ~ "Up",
        adj.p.value < p_thresh & MeanDiff <= -diff_thresh ~ "Down",
        TRUE ~ "NS"
      ),
      neglog10p = -log10(pmax(adj.p.value, .Machine$double.xmin)),
      is_capped = .data$neglog10p > y_cap,
      neglog10p_plot = pmin(.data$neglog10p, y_cap),
      point_class = ifelse(.data$is_capped, "capped", "normal")
    )

  label_df <- df %>%
    dplyr::filter(Significance %in% c("Up", "Down")) %>%
    dplyr::group_by(Cluster) %>%
    dplyr::arrange(dplyr::desc(abs(MeanDiff)), adj.p.value) %>%
    dplyr::slice_head(n = 12) %>%
    dplyr::ungroup()

  soft_cols <- c("Up" = "#E63980", "Down" = "#4CC9F0", "NS" = "#C8C8C8")

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = MeanDiff, y = neglog10p_plot)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(color = Significance, fill = Significance, shape = point_class),
      alpha = 0.85,
      size = 2.4,
      stroke = 0.35
    ) +
    ggplot2::scale_shape_manual(
      values = c(normal = 16, capped = 24),
      guide = "none"
    ) +
    ggplot2::scale_color_manual(values = soft_cols) +
    ggplot2::scale_fill_manual(values = soft_cols, guide = "none") +
    ggplot2::geom_vline(
      xintercept = c(-diff_thresh, diff_thresh),
      linetype = "dashed", color = "grey40", linewidth = 0.4
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(p_thresh),
      linetype = "dotted", color = "grey40", linewidth = 0.4
    ) +
    ggplot2::geom_hline(
      yintercept = y_cap,
      linetype = "solid", color = "grey55", linewidth = 0.35
    ) +
    ggrepel::geom_text_repel(
      data = label_df,
      ggplot2::aes(label = Protein),
      size = max(8 / (72.27 / 25.4), base_size / 4),
      max.overlaps = 40,
      box.padding = 0.35,
      point.padding = 0.25,
      segment.color = "grey60",
      min.segment.length = 0
    ) +
    ggplot2::coord_cartesian(ylim = c(0, y_cap * 1.02)) +
    ggplot2::facet_wrap(~ Cluster, scales = "free_x", ncol = ncol) +
    ggplot2::labs(
      title = "Volcano Plots of Differential Proteins by Neighborhood Cluster",
      x = "Mean Difference (Target - Control)",
      y = expression(-log[10](adjusted ~ italic(P))),
      color = "Significance",
      caption = paste0(
        "Points with -log10(adj.P) > ", y_cap,
        " are truncated and shown as triangles on the cap line"
      )
    ) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = "black", size = base_size),
      plot.title = ggplot2::element_text(face = "plain", size = max(8, base_size + 2), hjust = 0.5),
      plot.caption = ggplot2::element_text(size = max(8, base_size - 3), colour = "grey40"),
      axis.title = ggplot2::element_text(face = "plain", size = max(8, base_size)),
      axis.text = ggplot2::element_text(face = "plain", size = max(8, base_size - 2), colour = "black"),
      strip.text = ggplot2::element_text(face = "plain", size = max(8, base_size)),
      strip.background = ggplot2::element_rect(fill = "grey95", colour = "grey80"),
      legend.position = "right",
      legend.title = ggplot2::element_text(face = "plain", size = max(8, base_size)),
      legend.text = ggplot2::element_text(size = max(8, base_size - 2)),
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.7)
    )

  if (save_plot) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    ggplot2::ggsave(
      file.path(output_dir, paste0(filename, ".pdf")), p,
      width = width, height = height
    )
    ggplot2::ggsave(
      file.path(output_dir, paste0(filename, ".png")), p,
      width = width, height = height, dpi = 300
    )
    message("Plots saved to: ", output_dir)
  }

  return(p)
}

#' Perform cluster-specific protein enrichment analysis
#'
#' @param diff_results Differential expression results
#' @param protein_mapping Optional protein to gene symbol mapping
#' @param species Species for database selection (default: "human")
#' @param protein_databases Character vector of EnrichR databases (optional)
#' @param custom_databases Custom databases to include (optional)
#' @param pvalueCutoff Significance cutoff (default: 0.05)
#' @param mean_diff_cutoff Minimum mean difference, Target - Control (default: 0)
#' @param use_adj_pvalue Whether to use adjusted p-values (default: TRUE)
#' @param background Character vector of background gene symbols. If NULL
#'   (default), uses all unique proteins in `diff_results` (recommended for
#'   targeted protein panels). Genome-wide EnrichR defaults are inappropriate
#'   when only a small panel was assayed.
#' @param include_overlap logical; pass to enrichR when using background
#'   (default: TRUE)
#' @return Data frame with enrichment results
#' @examples
#' \dontrun{
#' # Requires clusterProfiler and org.Hs.eg.db:
#' # res <- perform_differential_expression(Sphinx:::.sphinx_example_protein_df(40))
#' # enrich <- perform_cluster_enrichment(res)
#' }
#' @export
perform_cluster_enrichment <- function(
    diff_results,
    protein_mapping = NULL,
    species = "human",
    protein_databases = NULL,
    custom_databases = NULL,
    pvalueCutoff = 0.05,
    mean_diff_cutoff = 0,
    use_adj_pvalue = TRUE,
    background = NULL,
    include_overlap = TRUE
) {
  required_pkgs <- c("dplyr", "enrichR")
  miss <- setdiff(required_pkgs, rownames(installed.packages()))
  if (length(miss)) stop("Please install: ", paste(miss, collapse = ", "))

  if (is.null(protein_databases)) {
    if (species == "human") {
      protein_databases <- c("CORUM", "Reactome_2022", "KEGG_2021_Human",
                             "GO_Biological_Process_2021", "GO_Molecular_Function_2021")
    } else if (species == "mouse") {
      protein_databases <- c("Mouse_Gene_Atlas", "Reactome_2022", "KEGG_2021_Mouse",
                             "GO_Biological_Process_2021", "GO_Molecular_Function_2021")
    } else if (species == "rat") {
      protein_databases <- c("Reactome_2022", "KEGG_2021_Rat",
                             "GO_Biological_Process_2021", "GO_Molecular_Function_2021")
    } else {
      protein_databases <- c("CORUM", "Reactome_2022", "KEGG_2021_Human",
                             "GO_Biological_Process_2021", "GO_Molecular_Function_2021")
      warning("Unknown species '", species, "', using human databases")
    }
  }

  if (!is.null(custom_databases)) {
    protein_databases <- unique(c(protein_databases, custom_databases))
  }

  need <- c("Protein", "Cluster", "MeanDiff", "p.value", "adj.p.value")
  if (!all(need %in% colnames(diff_results))) {
    stop("Missing columns: ", paste(setdiff(need, colnames(diff_results)), collapse = ", "))
  }

  if (!is.null(protein_mapping)) {
    diff_results <- diff_results %>%
      dplyr::mutate(Original_Protein = Protein) %>%
      dplyr::left_join(protein_mapping, by = c("Protein" = "Original_Name")) %>%
      dplyr::mutate(Protein = dplyr::coalesce(Gene_Symbol, Protein))
  }

  # Background: all tested proteins (panel), not EnrichR genome-wide default
  if (is.null(background)) {
    background <- unique(as.character(diff_results$Protein))
    background <- background[!is.na(background) & nzchar(background)]
  } else {
    background <- unique(as.character(background))
    background <- background[!is.na(background) & nzchar(background)]
  }
  if (length(background) < 2L) {
    stop("EnrichR background must contain at least 2 gene symbols")
  }
  message(
    "EnrichR background: ", length(background),
    " genes (tested protein panel / user-supplied)"
  )

  clusters <- unique(diff_results$Cluster)
  long_tbl <- dplyr::tibble()

  dbs_avail <- tryCatch(
    enrichR::listEnrichrDbs()$libraryName,
    error = function(e) {
      warning("Cannot get available databases: ", e$message)
      return(character(0))
    }
  )

  dbs_use <- base::intersect(protein_databases, dbs_avail)

  if (length(dbs_use) == 0) {
    warning("No available databases, please check database names or network connection")
    return(dplyr::tibble())
  }

  message("Using databases: ", paste(dbs_use, collapse = ", "))

  for (cl in clusters) {
    sig_genes <- diff_results %>%
      dplyr::filter(
        Cluster == cl,
        abs(MeanDiff) > mean_diff_cutoff,
        if (use_adj_pvalue) adj.p.value < pvalueCutoff else p.value < pvalueCutoff
      ) %>%
      dplyr::pull(Protein) %>% unique()

    # Only genes present in the declared background
    sig_genes <- intersect(sig_genes, background)

    if (length(sig_genes) < 2) {
      message("Cluster ", cl, " has insufficient genes, skipping")
      next
    }

    res_list <- tryCatch(
      enrichR::enrichr(
        genes = sig_genes,
        databases = dbs_use,
        background = background,
        include_overlap = include_overlap
      ),
      error = function(e) {
        warning("Cluster ", cl, " enrichment failed: ", e$message)
        NULL
      }
    )

    if (is.null(res_list)) next

    for (db in names(res_list)) {
      db_res <- res_list[[db]] %>% as_tibble()

      if (nrow(db_res) == 0) next

      if (!"Term" %in% names(db_res)) {
        db_res$Term <- character(nrow(db_res))
      } else if (all(is.na(db_res$Term)) || is.logical(db_res$Term)) {
        db_res$Term <- as.character(db_res$Term)
      }

      required_cols <- c("Adjusted.P.value", "Overlap")
      for (col in required_cols) {
        if (!col %in% names(db_res)) {
          db_res[[col]] <- NA
        }
      }

      db_res <- db_res %>%
        mutate(
          Cluster = cl,
          Database = db,
          FDR = ifelse("Adjusted.P.value" %in% names(.), Adjusted.P.value, NA),
          GenesN = ifelse("Overlap" %in% names(.), as.numeric(gsub(".*/", "", Overlap)), NA),
          Term = as.character(Term)
        ) %>%
        filter(FDR < pvalueCutoff)

      if (nrow(long_tbl) == 0) {
        long_tbl <- db_res[0, ]
      }

      long_tbl <- bind_rows(long_tbl, db_res)
    }
  }

  class(long_tbl) <- c("tbl_df", "tbl", "data.frame")
  return(long_tbl)
}

#' Generate publication-quality enrichment visualization plots
#'
#' @param cluster_enrich Enrichment results data frame
#' @param top_n Number of top terms to show per cluster (default: 5)
#' @param fdr_cutoff FDR cutoff for filtering (default: 0.05)
#' @param term_trunc_length Length to truncate term names (default: 40)
#' @param base_font_size Base font size for plots (default: 7)
#' @param cluster_colors Optional vector of colors for clusters
#' @param plot_types Types of plots to generate (default: all)
#' @param save_plot Whether to save plots (default: FALSE)
#' @param output_dir Output directory for saving (default: "enrichment_plots")
#' @param width Plot width in inches when saving (default: 9)
#' @param height Plot height in inches when saving (default: 10)
#' @return List of ggplot objects
#' @examples
#' \dontrun{
#' # Requires enrichment results from perform_cluster_enrichment():
#' # plot_enrichment_results(enrich)
#' }
#' @export
plot_enrichment_results <- function(
    cluster_enrich,
    top_n = 5,
    fdr_cutoff = 0.05,
    term_trunc_length = 40,
    base_font_size = 12,
    cluster_colors = NULL,
    plot_types = c("bar", "heatmap", "dot"),
    save_plot = FALSE,
    output_dir = "enrichment_plots",
    width = 9,
    height = 10
) {

  required_pkgs <- c("dplyr", "ggplot2", "scales", "stringr")
  miss <- setdiff(required_pkgs, rownames(installed.packages()))
  if (length(miss)) stop("Please install required packages: ", paste(miss, collapse = ", "))

  required_cols <- c("Cluster", "FDR", "Term", "GenesN")
  if (!all(required_cols %in% colnames(cluster_enrich))) {
    stop("Missing required columns: ", paste(setdiff(required_cols, colnames(cluster_enrich)), collapse = ", "))
  }

  enrich_top <- cluster_enrich %>%
    dplyr::filter(FDR < fdr_cutoff) %>%
    dplyr::group_by(Cluster) %>%
    dplyr::arrange(FDR) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      log10FDR = -log10(FDR),
      Term_short = stringr::str_replace(Term, " \\(.*", ""),
      Term_short = stringr::str_trunc(Term_short, term_trunc_length)
    )

  if (nrow(enrich_top) == 0) {
    warning("No significant terms found after FDR filtering.")
    return(list())
  }

  base_font_size <- max(8, as.numeric(base_font_size))
  sci_theme <- ggplot2::theme_classic(base_size = base_font_size) +
    ggplot2::theme(
      text            = ggplot2::element_text(colour = "black", size = base_font_size),
      plot.title      = ggplot2::element_text(
        face = "plain", size = max(8, base_font_size + 2),
        hjust = 0.5, margin = ggplot2::margin(b = 6)
      ),
      axis.title      = ggplot2::element_text(face = "plain", size = max(8, base_font_size + 1)),
      axis.text.x     = ggplot2::element_text(
        size = max(8, base_font_size), angle = 45, hjust = 1, vjust = 1, colour = "black"
      ),
      axis.text.y     = ggplot2::element_text(size = max(8, base_font_size), colour = "black"),
      axis.ticks      = ggplot2::element_line(linewidth = 0.35),
      panel.grid      = ggplot2::element_blank(),
      panel.border    = ggplot2::element_rect(fill = NA, colour = "black", linewidth = 0.7),
      axis.line       = ggplot2::element_blank(),
      legend.text     = ggplot2::element_text(size = max(8, base_font_size - 1)),
      legend.title    = ggplot2::element_text(size = max(8, base_font_size), face = "plain"),
      legend.position = "bottom",
      strip.background= ggplot2::element_rect(fill = "grey95", colour = "grey80"),
      strip.text      = ggplot2::element_text(face = "plain", size = max(8, base_font_size)),
      plot.margin     = ggplot2::margin(6, 6, 6, 6, "pt")
    )

  n_clusters <- length(unique(enrich_top$Cluster))
  if (is.null(cluster_colors)) {
    soft_base <- c(
      "#E05C6E", "#4EA8DE", "#E8C04A", "#4CB87A", "#8B6BC9",
      "#E8884A", "#3DB8A0", "#C45BA0", "#A67C52", "#5B7FD6"
    )
    if (n_clusters <= length(soft_base)) {
      cluster_colors <- soft_base[seq_len(n_clusters)]
    } else {
      cluster_colors <- get_color_palette(n_clusters)
    }
  }
  names(cluster_colors) <- unique(as.character(enrich_top$Cluster))

  soft_seq <- c("#FFF5F8", "#FFE0EC", "#FFB3D1", "#FF8FAB", "#F06595", "#E63980", "#C9184A")

  plot_list <- list()

  if ("bar" %in% plot_types) {
    p_bar <- ggplot2::ggplot(
      enrich_top,
      ggplot2::aes(x = log10FDR, y = stats::reorder(Term_short, log10FDR), fill = Cluster)
    ) +
      ggplot2::geom_col(width = 0.75, alpha = 0.95) +
      ggplot2::facet_wrap(~ Cluster, scales = "free_y", ncol = 2) +
      ggplot2::scale_fill_manual(values = cluster_colors) +
      ggplot2::labs(
        x = expression(-log[10](FDR)),
        y = "Enriched Biological Terms",
        title = "Top Enriched Terms by Cluster"
      ) +
      sci_theme +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = base_font_size, lineheight = 0.9),
        legend.position = "none"
      )
    plot_list$bar_plot <- p_bar
  }

  if ("heatmap" %in% plot_types) {
    p_heatmap <- ggplot2::ggplot(
      enrich_top,
      ggplot2::aes(x = Cluster, y = stats::reorder(Term_short, log10FDR), fill = log10FDR)
    ) +
      ggplot2::geom_tile(color = "white", linewidth = 0.6, height = 0.92, width = 0.92) +
      ggplot2::scale_fill_gradientn(
        colours = soft_seq,
        name = expression(-log[10](FDR))
      ) +
      ggplot2::labs(
        x = "Cluster",
        y = "Enriched Biological Terms",
        title = "Enrichment Significance Heatmap"
      ) +
      sci_theme +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = base_font_size)
      )
    plot_list$heatmap <- p_heatmap
  }

  if ("dot" %in% plot_types) {
    p_dot <- ggplot2::ggplot(
      enrich_top,
      ggplot2::aes(
        x = Cluster, y = stats::reorder(Term_short, log10FDR),
        size = GenesN, color = log10FDR
      )
    ) +
      ggplot2::geom_point(alpha = 0.9) +
      ggplot2::scale_color_gradientn(
        colours = soft_seq,
        name = expression(-log[10](FDR))
      ) +
      ggplot2::scale_size_continuous(
        range = c(3.5, 11),
        name = "Gene Count",
        breaks = pretty(range(enrich_top$GenesN), n = 4)
      ) +
      ggplot2::labs(
        x = "Cluster",
        y = "Enriched Biological Terms",
        title = "Enrichment Overview"
      ) +
      sci_theme +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(size = base_font_size),
        legend.box = "horizontal"
      )
    plot_list$dot_plot <- p_dot
  }

  if (save_plot && length(plot_list) > 0) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    for (plot_name in names(plot_list)) {
      ggplot2::ggsave(
        file.path(output_dir, paste0(plot_name, ".pdf")),
        plot_list[[plot_name]], width = width, height = height
      )
      ggplot2::ggsave(
        file.path(output_dir, paste0(plot_name, ".png")),
        plot_list[[plot_name]], width = width, height = height, dpi = 300
      )
    }
    message("Plots saved to: ", output_dir)
  }

  return(plot_list)
}
