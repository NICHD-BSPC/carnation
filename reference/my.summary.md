# Summarize DESeq2 results into a dataframe

summary(res) prints out info; this function captures it into a dataframe

## Usage

``` r
my.summary(res, dds, alpha, lfc.thresh = 0)
```

## Arguments

- res:

  DESeq2 results object

- dds:

  DEseq2 object

- alpha:

  Alpha level at which to call significantly changing genes

- lfc.thresh:

  log2FoldChange threshold

## Value

Dataframe of summarized results

## Examples

``` r
library(DESeq2)

# make example DESeq data set
dds <- makeExampleDESeqDataSet()

# run DESeq2
dds <- DESeq(dds)
#> estimating size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing

# make comparisons
res <- results(dds, contrast=c('condition', 'A', 'B'))

# get summary
df <- my.summary(res, dds, alpha=0.1)
```
