# Create a labeled MA plot

This function creates an MA plot from a data.frame containing DE
analysis results.

## Usage

``` r
plotMA.label(
  res,
  fdr.thres = 0.01,
  fc.thres = 0,
  fc.lim = NULL,
  lab.genes = NULL,
  tolower.cols = c("SYMBOL", "ALIAS")
)
```

## Arguments

- res:

  data.frame with DE analysis results. Must contain "padj" &
  "log2FoldChange" columns

- fdr.thres:

  False discovery rate (FDR) threshold

- fc.thres:

  log2FoldChange threshold

- fc.lim:

  y-axis limits

- lab.genes:

  genes to label on MA plot

- tolower.cols:

  column names that will be converted to lower case

## Value

ggplot handle

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

# add gene and symbol column
res$gene <- rownames(res)
res$symbol <- rownames(res)

plotMA.label(res, lab.genes = c("gene1", "gene2"))
#> Warning: log-10 transformation introduced infinite values.

```
