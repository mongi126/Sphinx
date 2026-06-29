# R/functional_analysis.R

#' Prepare protein expression data for functional analysis
#'
#' @param cluster_df Data frame with neighborhood cluster assignments and spatial coordinates
#' @param expr_df Data frame with protein expression data (rows = cells, columns = proteins)
#' @return Merged data frame containing spatial, cluster, and protein expression data
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

#' Perform differential expression analysis
#'
#' Compares each neighborhood cluster to all other cells. For CLR-normalized
#' protein values, reports mean difference (Target - Control) rather than
#' log2 fold change to avoid redundant log transforms.
#'
#' @param protein_df Data frame with protein expression and cluster information
#' @return Data frame with columns Mean_Target, Mean_Control, MeanDiff, p.value,
#'   adj.p.value, and Significance
#' @export
perform_differential_expression <- function(protein_df) {
  meta_cols <- c("Cell_ID", "X", "Y", "Neighborhood_Cluster", "Spatial_Zone")
  protein_cols <- setdiff(colnames(protein_df), meta_cols)

  clusters <- unique(protein_df$Neighborhood_Cluster)
  results <- list()

  for (cluster in clusters) {
    protein_df$In_Cluster <- ifelse(
      protein_df$Neighborhood_Cluster == cluster,
      "Target",
      "Control"
    )

    cluster_results <- lapply(protein_cols, function(prot) {
      target <- protein_df[protein_df$In_Cluster == "Target", prot]
      control <- protein_df[protein_df$In_Cluster == "Control", prot]

      if (length(target) < 2 || length(control) < 2) {
        return(NULL)
      }

      test_res <- t.test(target, control)
      mean_target <- mean(target, na.rm = TRUE)
      mean_control <- mean(control, na.rm = TRUE)
      data.frame(
        Protein = prot,
        Cluster = cluster,
        Mean_Target = mean_target,
        Mean_Control = mean_control,
        MeanDiff = mean_target - mean_control,
        p.value = test_res$p.value
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
#' @param diff_results Differential expression results
#' @param diff_thresh Mean difference threshold, Target - Control (default: 0.25)
#' @param p_thresh P-value threshold (default: 0.05)
#' @param save_plot Whether to save the plot (default: FALSE)
#' @param output_dir Output directory for saving (default: "plots")
#' @param filename Output filename (default: "volcano_plots")
#' @return ggplot object
#' @export
plot_volcano_all_clusters <- function(diff_results, diff_thresh = 0.25, p_thresh = 0.05,
                                      save_plot = FALSE, output_dir = "plots",
                                      filename = "volcano_plots") {
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
    dplyr::mutate(Significance = dplyr::case_when(
      adj.p.value < p_thresh & MeanDiff >= diff_thresh ~ "Up",
      adj.p.value < p_thresh & MeanDiff <= -diff_thresh ~ "Down",
      TRUE ~ "NS"
    ))

  label_df <- df %>%
    dplyr::filter(Significance %in% c("Up", "Down"))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = MeanDiff, y = -log10(adj.p.value))) +
    ggplot2::geom_point(ggplot2::aes(color = Significance), alpha = 0.8, size = 2) +
    ggplot2::scale_color_manual(values = c("Up" = "firebrick", "Down" = "steelblue", "NS" = "grey70")) +
    ggplot2::geom_vline(xintercept = c(-diff_thresh, diff_thresh), linetype = "dashed", color = "black") +
    ggplot2::geom_hline(yintercept = -log10(p_thresh), linetype = "dotted", color = "black") +
    ggrepel::geom_text_repel(
      data = label_df,
      ggplot2::aes(label = Protein),
      size = 3,
      max.overlaps = 100,
      box.padding = 0.5,
      point.padding = 0.3,
      segment.color = "grey50"
    ) +
    ggplot2::facet_wrap(~ Cluster, scales = "free_y", ncol = 3) +
    ggplot2::labs(
      title = "Volcano Plots of Differential Proteins by Neighborhood Cluster",
      x = "Mean Difference (Target - Control)", y = "-log10(Adjusted p-value)",
      color = "Significance"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )

  if (save_plot) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    ggplot2::ggsave(file.path(output_dir, paste0(filename, ".pdf")), p, width = 12, height = 8)
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
#' @return Data frame with enrichment results
#' @export
perform_cluster_enrichment <- function(
    diff_results,
    protein_mapping = NULL,
    species = "human",
    protein_databases = NULL,
    custom_databases = NULL,
    pvalueCutoff = 0.05,
    mean_diff_cutoff = 0,
    use_adj_pvalue = TRUE
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

  clusters <- unique(diff_results$Cluster)
  long_tbl <- tibble()

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
    return(tibble())
  }

  message("Using databases: ", paste(dbs_use, collapse = ", "))

  for (cl in clusters) {
    sig_genes <- diff_results %>%
      filter(
        Cluster == cl,
        abs(MeanDiff) > mean_diff_cutoff,
        if (use_adj_pvalue) adj.p.value < pvalueCutoff else p.value < pvalueCutoff
      ) %>%
      pull(Protein) %>% unique()

    if (length(sig_genes) < 2) {
      message("Cluster ", cl, " has insufficient genes, skipping")
      next
    }

    res_list <- tryCatch(
      enrichR::enrichr(sig_genes, dbs_use),
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
#' @return List of ggplot objects
#' @export
plot_enrichment_results <- function(
    cluster_enrich,
    top_n = 5,
    fdr_cutoff = 0.05,
    term_trunc_length = 40,
    base_font_size = 7,
    cluster_colors = NULL,
    plot_types = c("bar", "heatmap", "dot"),
    save_plot = FALSE,
    output_dir = "enrichment_plots"
) {

  required_pkgs <- c("dplyr", "ggplot2", "scales", "stringr", "viridis")
  miss <- setdiff(required_pkgs, rownames(installed.packages()))
  if (length(miss)) stop("Please install required packages: ", paste(miss, collapse = ", "))

  required_cols <- c("Cluster", "FDR", "Term", "GenesN")
  if (!all(required_cols %in% colnames(cluster_enrich))) {
    stop("Missing required columns: ", paste(setdiff(required_cols, colnames(cluster_enrich)), collapse = ", "))
  }

  enrich_top <- cluster_enrich %>%
    filter(FDR < fdr_cutoff) %>%
    group_by(Cluster) %>%
    arrange(FDR) %>%
    slice_head(n = top_n) %>%
    ungroup() %>%
    mutate(
      log10FDR = -log10(FDR),
      Term_short = str_replace(Term, " \\(.*", ""),
      Term_short = str_trunc(Term_short, term_trunc_length)
    )

  if (nrow(enrich_top) == 0) {
    warning("No significant terms found after FDR filtering.")
    return(list())
  }

  sci_theme <- theme_bw(base_size = base_font_size) +
    theme(
      text            = element_text(family = "Arial", colour = "black"),
      plot.title      = element_text(face = "bold", size = base_font_size + 1,
                                     hjust = 0.5, margin = margin(b = 3)),
      axis.title      = element_text(face = "bold", size = base_font_size),
      axis.text.x     = element_text(size = base_font_size, angle = 45, hjust = 1, vjust = 1),
      axis.text.y     = element_text(size = base_font_size),
      axis.ticks      = element_line(linewidth = 0.3),
      panel.grid      = element_blank(),
      panel.border    = element_rect(fill = NA, linewidth = 0.3),
      legend.text     = element_text(size = base_font_size - 0.5),
      legend.title    = element_text(size = base_font_size, face = "bold"),
      legend.box      = "horizontal",
      legend.position = "bottom",
      strip.background= element_rect(fill = NA, colour = NA),
      strip.text      = element_text(face = "bold", size = base_font_size),
      plot.margin     = margin(2, 2, 2, 2, "pt")
    )

  n_clusters <- length(unique(enrich_top$Cluster))
  if (is.null(cluster_colors)) {
    cluster_colors <- hue_pal()(max(n_clusters, 30))
  }

  plot_list <- list()

  if ("bar" %in% plot_types) {
    p_bar <- ggplot(enrich_top, aes(x = log10FDR, y = reorder(Term_short, log10FDR), fill = Cluster)) +
      geom_col(width = 0.8, alpha = 0.9) +
      facet_wrap(~ Cluster, scales = "free_y", ncol = 2) +
      scale_fill_manual(values = setNames(cluster_colors, unique(enrich_top$Cluster))) +
      labs(
        x = expression(-log[10](FDR)),
        y = "Enriched Biological Terms",
        title = "Top Enriched Terms by Cluster"
      ) +
      sci_theme +
      theme(
        axis.text.y = element_text(size = base_font_size + 2, lineheight = 0.8),
        legend.position = "none"
      )
    plot_list$bar_plot <- p_bar
  }

  if ("heatmap" %in% plot_types) {
    p_heatmap <- enrich_top %>%
      ggplot(aes(x = Cluster, y = reorder(Term_short, log10FDR), fill = log10FDR)) +
      geom_tile(color = "white", linewidth = 0.8, height = 0.9, width = 0.9) +
      scale_fill_viridis_c(option = "inferno", name = expression(-log[10](FDR))) +
      labs(
        x = "Cluster",
        y = "Enriched Biological Terms",
        title = "Enrichment Significance Heatmap"
      ) +
      sci_theme +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = base_font_size + 2)
      )
    plot_list$heatmap <- p_heatmap
  }

  if ("dot" %in% plot_types) {
    p_dot <- enrich_top %>%
      ggplot(aes(x = Cluster, y = reorder(Term_short, log10FDR), size = GenesN, color = log10FDR)) +
      geom_point(alpha = 0.85) +
      scale_color_viridis_c(option = "turbo", name = expression(-log[10](FDR))) +
      scale_size_continuous(
        range = c(3, 10),
        name = "Gene Count",
        breaks = pretty(range(enrich_top$GenesN), n = 4)
      ) +
      labs(
        x = "Cluster",
        y = "Enriched Biological Terms",
        title = "Enrichment Overview"
      ) +
      sci_theme +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = base_font_size + 2),
        legend.box = "horizontal"
      )
    plot_list$dot_plot <- p_dot
  }

  if (save_plot && length(plot_list) > 0) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

    for (plot_name in names(plot_list)) {
      ggsave(file.path(output_dir, paste0(plot_name, ".pdf")),
             plot_list[[plot_name]], width = 8, height = 10)
    }
    message("Plots saved to: ", output_dir)
  }

  return(plot_list)
}
