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
