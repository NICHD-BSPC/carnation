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

  GeneTonic object for comparison 1

- res_enrich2:

  GeneTonic object for comparison 2 (default = NULL)

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

## Examples

``` r
library(GeneTonic)

# get DESeqResults object
data(res_dex, package='carnation')

# get enrichResult object
data(eres_dex, package='carnation')

# convert to GeneTonic object
gt <- shake_enrichResult(eres_dex)
#> Found 2186 gene sets in `enrichResult` object, of which 2186 are significant.
#> Converting for usage in GeneTonic...

# get annotation df
idx <- match(c('gene','symbol'), tolower(colnames(res_dex)))
anno_df <- res_dex[,idx]
colnames(anno_df) <- c('gene_id', 'gene_name')

# add aggregate score columns
gt <- get_aggrscores(gt, res_dex, anno_df)

# make radar plot
p <- gs_radar(gt)
```
