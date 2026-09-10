# Analyze spatial interactions between cell types

Counts cell-type contact edges and optionally normalizes by type
abundance so frequent types do not dominate the interaction matrix.

## Usage

``` r
analyze_spatial_interactions(
  df,
  edges,
  celltype_col = "celltype",
  normalize_abundance = TRUE
)
```

## Arguments

  - df:
    
    Spatial data with cell types

  - edges:
    
    Spatial network edges (data.frame with "from" and "to")

  - celltype\_col:
    
    Column name for cell types in metadata.

  - normalize\_abundance:
    
    logical; if TRUE (default), also compute abundance-normalized and
    enrichment matrices.

## Value

List with `interaction_matrix` (raw counts), `network`, and when
`normalize_abundance = TRUE` also `abundance_normalized`,
`enrichment_matrix`, and `type_abundance`.

## Details

Normalization:

  - `abundance_normalized`: \\(c\_{ij} / (n\_i n\_j)\\)

  - `enrichment_matrix`: observed / expected under random labeling of
    the same graph, where \\(E\_{ij} = E \\cdot 2 p\_i p\_j\\) (i \!= j)
    and \\(E\_{ii} = E \\cdot p\_i^2\\), with \\(p\_i = n\_i / N\\).

## Examples

``` r
df <- prepare_data(Sphinx:::.sphinx_example_df(40))
edges <- build_spatial_network(df, method = "knn", n_neighbors = 5, verbose = FALSE)
intx <- analyze_spatial_interactions(df, edges)
dim(intx$interaction_matrix)
```
