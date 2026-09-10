# Calculate an adaptive window size from cell density

Iteratively adjusts window side length so that the mean number of cells
per window falls between `min_cells` and `max_cells`.

## Usage

``` r
calculate_optimal_window_size(
  df,
  x_col = "X",
  y_col = "Y",
  min_cells = 30L,
  max_cells = 50L,
  range_divisor = 100L,
  max_iter = 100L
)
```

## Arguments

  - df:
    
    Spatial data with coordinates

  - x\_col:
    
    Column name for X coordinates (default: "X")

  - y\_col:
    
    Column name for Y coordinates (default: "Y")

  - min\_cells:
    
    Target minimum mean cells per window (default: 30)

  - max\_cells:
    
    Target maximum mean cells per window (default: 50)

  - range\_divisor:
    
    Divisor for initial window size from spatial range (default: 100)

  - max\_iter:
    
    Maximum adjustment iterations (default: 100)

## Value

List with window\_size, sliding\_step, mean\_cells, x\_range, y\_range

## Examples

``` r
df <- prepare_data(Sphinx:::.sphinx_example_df(50))
ws <- calculate_optimal_window_size(df)
ws$window_size
```
