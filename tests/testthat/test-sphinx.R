test_that("example spatial data helper works", {
  df <- Sphinx:::.sphinx_example_df(20)
  expect_equal(nrow(df), 20L)
  expect_true(all(c("Cell_ID", "X", "Y", "celltype") %in% names(df)))
})

test_that("prepare_data standardizes spatial input", {
  df <- Sphinx:::.sphinx_example_df(20)
  out <- prepare_data(df)
  expect_true(data.table::is.data.table(out))
  expect_true("Cell_ID" %in% names(out))
})

test_that("build_spatial_network returns edges", {
  df <- prepare_data(Sphinx:::.sphinx_example_df(30))
  edges <- build_spatial_network(df, method = "knn", k = 5, verbose = FALSE)
  expect_true(nrow(edges) > 0)
  expect_true(all(c("from", "to", "dist") %in% names(edges)))
})

test_that("annotate_celltypes adds celltype metadata", {
  obj <- Sphinx:::.sphinx_example_seurat(20)
  clusters <- levels(obj@meta.data$seurat_clusters)
  obj <- annotate_celltypes(
    obj,
    cluster_ids = clusters,
    celltype_labels = paste0("Type", seq_along(clusters))
  )
  expect_true("celltype" %in% colnames(obj@meta.data))
})

test_that("perform_differential_expression returns results", {
  df <- Sphinx:::.sphinx_example_protein_df(30)
  res <- perform_differential_expression(df)
  expect_true(nrow(res) > 0)
  expect_true(all(c("Protein", "Cluster", "MeanDiff") %in% names(res)))
})
