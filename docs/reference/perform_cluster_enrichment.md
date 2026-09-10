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
  mean_diff_cutoff = 0,
  use_adj_pvalue = TRUE,
  background = NULL,
  include_overlap = TRUE
)
```

## Arguments

  - diff\_results:
    
    Differential expression results

  - protein\_mapping:
    
    Optional protein to gene symbol mapping

  - species:
    
    Species for database selection (default: "human")

  - protein\_databases:
    
    Character vector of EnrichR databases (optional)

  - custom\_databases:
    
    Custom databases to include (optional)

  - pvalueCutoff:
    
    Significance cutoff (default: 0.05)

  - mean\_diff\_cutoff:
    
    Minimum mean difference, Target - Control (default: 0)

  - use\_adj\_pvalue:
    
    Whether to use adjusted p-values (default: TRUE)

  - background:
    
    Character vector of background gene symbols. If NULL (default), uses
    all unique proteins in `diff_results` (recommended for targeted
    protein panels). Genome-wide EnrichR defaults are inappropriate when
    only a small panel was assayed.

  - include\_overlap:
    
    logical; pass to enrichR when using background (default: TRUE)

## Value

Data frame with enrichment results

## Examples

``` r
if (FALSE) { # \dontrun{
# Requires clusterProfiler and org.Hs.eg.db:
# res <- perform_differential_expression(Sphinx:::.sphinx_example_protein_df(40))
# enrich <- perform_cluster_enrichment(res)
} # }
```
