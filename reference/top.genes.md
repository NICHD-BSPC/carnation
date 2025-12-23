# Get top DE genes by log2FoldChange or adjusted p-value

Get top DE genes by log2FoldChange or adjusted p-value

## Usage

``` r
top.genes(res, fdr.thres = 0.01, fc.thres = 0, n = 10, by = "log2FoldChange")
```

## Arguments

- res:

  data.frame with DE analysis results

- fdr.thres:

  FDR threshold

- fc.thres:

  log2FoldChange threshold

- n:

  number of genes to return

- by:

  metric to determine top genes ('log2FoldChange' or 'padj')

## Value

vector of gene symbols

## Examples

``` r
library(DESeq2)

# make example DESeq dataset
dds <- makeExampleDESeqDataSet()

# run DE analysis
dds <- DESeq(dds)
#> estimating size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing

# extract comparison of interest
res <- results(dds, contrast = c("condition", "A", "B"))

g <- top.genes(res)
```
