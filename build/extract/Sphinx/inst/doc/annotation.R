## ----load-packages, eval=FALSE------------------------------------------------
# library(Sphinx)
# library(Seurat)
# library(ggplot2)
# library(pheatmap)
# library(dplyr)

## ----load-data, eval=FALSE----------------------------------------------------
# # Load preprocessed spatial data
# codex.obj <- readRDS("cycif1_processed.rds")
# 
# # Check object structure
# print(codex.obj)

## ----find-markers, eval=FALSE-------------------------------------------------
# # Identify top marker proteins for each cluster
# top5 <- find_top_markers(
#   codex.obj,
#   assay = "RNA",
#   save_path = "top5proteins.csv"
# )
# 
# # View top markers
# head(top5)

## ----violin-plots, eval=FALSE-------------------------------------------------
# # Create violin plots for all detected proteins
# plot_marker_violin(
#   codex.obj,
#   markers = rownames(codex.obj[["RNA"]]),
#   assay = "RNA",
#   save_path = "allmarker.png"
# )

## ----marker-heatmap, eval=FALSE-----------------------------------------------
# # Create heatmap of top marker expression
# plot_marker_heatmap(
#   codex.obj,
#   top_markers = top5,
#   save_path = "marker_heatmap.png"
# )

## ----umap-markers, eval=FALSE-------------------------------------------------
# plot_umap_markers(codex.obj,
#                   markers = c("CD31", "CD19", "Pan-Cytokeratin",
#                               "a-SMA","CD8","CD4"))

## ----spatial-markers, eval=FALSE----------------------------------------------
# plot_spatial_markers(codex.obj,
#                   markers = c("CD31", "CD19", "Pan-Cytokeratin",
#                               "a-SMA","CD8","CD4"),
#                   point_size = 0.5)

## ----annotate_celltypes, eval=FALSE-------------------------------------------
# # Define cluster to cell type mapping
# cluster_ids <- c(0, 1, 2, 3, 4, 5, 6, 7，8，
#                  9，10，11，12，13，14，15)
# celltype_labels <- c(
#   "Endothelial", "Tumor", "DC", "M2_Macrophage","Tumor","CAF","PD-1_CD4_T","B","Endothelial",
#   "GC_B", "DC", "EMT_like_Tumor", "Proliferating_DC","Endothelial","Endothelial","Proliferating_Tumor"
# )
# 
# # Apply cell type annotations
# codex.obj <- annotate_celltypes(
#   codex.obj,
#   cluster_ids = cluster_ids,
#   celltype_labels = celltype_labels,
#   cluster_column = "RNA_snn_res.0.5"
# )
# 
# # Check annotation results
# table(codex.obj$celltype)

## ----annotated-umap, eval=FALSE-----------------------------------------------
# # Create UMAP colored by cell types
# plot_annotated_umap(
#   codex.obj,
#   save_path = "umap_by_celltype.png"
# )

## ----spatial-distribution, eval=FALSE-----------------------------------------
# # Visualize spatial distribution of cell types
# plot_spatial_distribution(
#   codex.obj,
#   save_path = "celltype_spatial_plot.png"
# )

## ----save-results, eval=FALSE-------------------------------------------------
# # Save annotated object
# saveRDS(codex.obj, file = "tsu33_celltype.rds")
# 
# # Export metadata for downstream analysis
# write.csv(codex.obj@meta.data, file = "tsu33_metadata.csv")

