# Stable cell-type -\> color mapping (same label =\> same color across plots)

By default colors are assigned from a fixed soft pastel pool via a
deterministic string hash, so a cell type keeps its color even when
other types are absent. If `all_levels` is supplied, colors are assigned
by alphabetical order over that full set (preferred when the complete
type list is known).

## Usage

``` r
assign_celltype_colors(levels, all_levels = NULL)
```

## Arguments

  - levels:
    
    Character vector of cell-type / category labels to color

  - all\_levels:
    
    Optional full set of labels for ordered assignment

## Value

Named character vector of hex colors

## Examples

``` r
assign_celltype_colors(c("B", "T", "Macrophage"))
```
