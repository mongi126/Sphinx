# Internal helpers for package examples and tests
# @keywords internal

.sphinx_example_df <- function(n = 40L, seed = 1L) {
  set.seed(seed)
  data.table::data.table(
    Cell_ID = paste0("cell", seq_len(n)),
    X = stats::runif(n, 0, 100),
    Y = stats::runif(n, 0, 100),
    celltype = sample(c("T", "B", "Myeloid"), n, replace = TRUE)
  )
}

.sphinx_example_seurat <- function(n = 50L, nfeatures = 10L, seed = 1L) {
  set.seed(seed)
  cell_lambdas <- sample(c(3, 8, 15), n, replace = TRUE)
  counts <- Matrix::Matrix(
    stats::rpois(n * nfeatures, lambda = rep(cell_lambdas, each = nfeatures)),
    nrow = nfeatures,
    sparse = TRUE,
    dimnames = list(
      paste0("prot", seq_len(nfeatures)),
      paste0("cell", seq_len(n))
    )
  )
  obj <- Seurat::CreateSeuratObject(counts = counts, min.features = 0, min.cells = 0)
  obj$X <- stats::runif(n, 0, 100)
  obj$Y <- stats::runif(n, 0, 100)
  emb <- as.matrix(obj@meta.data[, c("X", "Y")])
  colnames(emb) <- c("SPATIAL_1", "SPATIAL_2")
  obj[["spatial"]] <- Seurat::CreateDimReducObject(
    embeddings = emb,
    assay = Seurat::DefaultAssay(obj),
    key = "SPATIAL_"
  )
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  Seurat::VariableFeatures(obj) <- rownames(obj)
  obj <- Seurat::ScaleData(obj, verbose = FALSE)
  obj <- Seurat::RunPCA(obj, verbose = FALSE, npcs = min(5L, nfeatures - 1L))
  obj$seurat_clusters <- factor(sample(0:2, n, replace = TRUE))
  Seurat::Idents(obj) <- obj$seurat_clusters
  umap_mat <- matrix(
    stats::rnorm(n * 2),
    ncol = 2,
    dimnames = list(colnames(obj), c("UMAP_1", "UMAP_2"))
  )
  obj[["umap"]] <- Seurat::CreateDimReducObject(
    embeddings = umap_mat,
    assay = Seurat::DefaultAssay(obj),
    key = "UMAP_"
  )
  obj
}

.sphinx_example_protein_df <- function(n = 40L, seed = 1L) {
  set.seed(seed)
  data.frame(
    Cell_ID = paste0("cell", seq_len(n)),
    X = stats::runif(n, 0, 100),
    Y = stats::runif(n, 0, 100),
    Neighborhood_Cluster = sample(1:3, n, replace = TRUE),
    PD.1 = stats::rnorm(n, mean = 2, sd = 0.5),
    CD3 = stats::rnorm(n, mean = 1.5, sd = 0.5),
    Spatial_Zone = as.character(sample(1:3, n, replace = TRUE)),
    stringsAsFactors = FALSE
  )
}
