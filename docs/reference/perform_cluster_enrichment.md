# Perform cluster-specific protein enrichment analysis

Perform cluster-specific protein enrichment analysis

## Usage

``` r
perform_cluster_enrichment(
  diff_results,
  protein_mapping = NULL,
  species = "human",
  protein_databases = NULL,
  custom_databases = NULL,
  pvalueCutoff = 0.05,
  log2FC_cutoff = 0,
  use_adj_pvalue = TRUE
)
```

## Arguments

- diff_results:

  Differential expression results

- protein_mapping:

  Optional protein to gene symbol mapping

- species:

  Species for database selection (default: "human")

- protein_databases:

  Character vector of EnrichR databases (optional)

- custom_databases:

  Custom databases to include (optional)

- pvalueCutoff:

  Significance cutoff (default: 0.05)

- log2FC_cutoff:

  Minimum log2 fold change (default: 0)

- use_adj_pvalue:

  Whether to use adjusted p-values (default: TRUE)

## Value

Data frame with enrichment results
