# Radar plot

This is a copy of gs_radar from GeneTonic where the labels of gene sets
are converted to parameters

## Usage

``` r
gs_radar(
  res_enrich,
  res_enrich2 = NULL,
  label1 = "scenario 1",
  label2 = "scenario 2",
  n_gs = 20,
  p_value_column = "gs_pvalue"
)
```

## Arguments

- res_enrich:

  DE results from comparison 1

- res_enrich2:

  DE results from comparison 2

- label1:

  label for comparison 1

- label2:

  label for comparison 2

- n_gs:

  number of gene sets (default = 20)

- p_value_column:

  column to use as p-value (default = 'gs_pvalue')

## Value

ggplot handle
