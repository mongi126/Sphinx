## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, comment = "#>")

## ----install, eval=FALSE------------------------------------------------------
# devtools::install_github("mongi126/Sphinx")
# library(Sphinx)
# packageVersion("Sphinx")

## ----quick-start, eval=FALSE--------------------------------------------------
# library(Sphinx)
# library(Seurat)
# 
# # 1. Preprocessing
# obj <- load_spatial_data("your_data.csv")
# obj <- filter_data(obj)
# obj <- process_data(obj)
# 
# # 2. Cell annotation
# markers <- find_top_markers(obj)
# obj     <- annotate_celltypes(obj, markers)
# 
# # 3. Spatial network
# obj <- prepare_data(obj)
# net <- build_spatial_network(obj, method = "knn")
# obj <- calculate_neighborhood_features(obj, net)
# 
# # 4. Functional analysis
# de_results <- perform_differential_expression(obj)
# enrichment <- perform_cluster_enrichment(de_results)

## ----workflow-figure, echo=FALSE, out.width="100%"----------------------------
knitr::include_graphics("../man/figures/workflow.jpg")

## ----help, eval=FALSE---------------------------------------------------------
# # Package documentation
# help(package = "Sphinx")
# 
# # Function reference
# ?load_spatial_data
# ?build_spatial_network

