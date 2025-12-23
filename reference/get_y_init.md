# Get initial y-axis limits

Get initial y-axis limits

## Usage

``` r
get_y_init(df, y_delta, pseudocount)
```

## Arguments

- df:

  data.frame with counts. Must have column 'count'

- y_delta:

  y-axis padding for visualization, must be between 0 and 1

- pseudocount:

  pseudo-count to add to the data.frame

## Value

min and max limits for count column, padded for visualization

## Examples

``` r
# make example DESeq dataset
dds <- DESeq2::makeExampleDESeqDataSet()

# get gene counts
df <- get_gene_counts(dds, gene = c('gene1', 'gene2'))

# get y axis limits
get_y_init(df, y_delta = 0.01, pseudocount = 1)
#> [1]   2 120
```
