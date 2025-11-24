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
