# R/imports.R

#' @importFrom Seurat CreateSeuratObject DefaultAssay NormalizeData ScaleData
#' @importFrom Seurat RunPCA FindNeighbors FindClusters RunUMAP
#' @importFrom Seurat FindAllMarkers VlnPlot AverageExpression DimPlot
#' @importFrom Seurat ElbowPlot Embeddings CreateDimReducObject
#' @importFrom Seurat GetAssayData SetAssayData VariableFeatures
#' @importFrom dplyr %>% group_by top_n filter mutate case_when arrange
#' @importFrom dplyr slice_head ungroup pull left_join coalesce bind_rows
#' @importFrom dplyr summarise n
#' @importFrom ggplot2 ggplot aes geom_point theme_classic scale_color_manual
#' @importFrom ggplot2 coord_fixed labs theme element_text guide_legend
#' @importFrom ggplot2 guides ggsave geom_col facet_wrap geom_tile
#' @importFrom ggplot2 scale_fill_manual scale_size_continuous
#' @importFrom ggplot2 scale_y_continuous geom_bar geom_line geom_segment
#' @importFrom ggplot2 element_blank element_rect geom_text
#' @importFrom pheatmap pheatmap
#' @importFrom RColorBrewer brewer.pal
#' @importFrom data.table as.data.table setnames copy .I .N .SD :=
#' @importFrom igraph graph_from_data_frame V E degree strength
#' @importFrom igraph cluster_louvain induced_subgraph vcount
#' @importFrom utils write.csv head tail str
#' @importFrom grDevices colorRampPalette dev.off pdf
#' @importFrom stats dist hclust kmeans prcomp complete.cases quantile
#' @importFrom stats sd var p.adjust t.test median mad setNames na.omit reorder
#' @importFrom utils installed.packages read.csv
#' @importFrom ggplot2 theme_bw element_line margin scale_fill_viridis_c scale_color_viridis_c
#' @importFrom viridis scale_fill_viridis scale_color_viridis viridis
#' @importFrom methods as is new
#' @importFrom scattermore geom_scattermore
#' @importFrom progress progress_bar
#' @importFrom RANN nn2
#' @importFrom Matrix Matrix
#' @importFrom gridExtra grid.arrange
#' @importFrom rlang .data sym
#' @importFrom tibble as_tibble tibble
#' @importFrom tidyr pivot_wider
#' @importFrom stringr str_replace str_trunc
#' @importFrom scales hue_pal percent_format rescale
#' @importFrom deldir deldir tile.list triang.list

NULL
