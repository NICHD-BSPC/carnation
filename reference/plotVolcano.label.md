# Create a static labeled Volcano plot

This function creates a static ggplot Volcano plot from a data.frame
containing DE analysis results. Mirrors plotVolcano.label_ly but renders
with ggplot2 for use in downloads.

## Usage

``` r
plotVolcano.label(
  res,
  fdr.thres = 0.1,
  fc.thres = 0,
  y.axis.lim = NULL,
  fc.lim = c(-4, 4),
  lab.genes = NULL,
  tolower.cols = c("SYMBOL", "ALIAS"),
  color_by = c("baseMean", "significance"),
  colorscale = "viridis",
  alpha = 0.6,
  pvalue_padj_switch = c("pvalue", "padj")
)
```

## Arguments

- res:

  data.frame with DE analysis results. Must contain "pvalue", "padj",
  "baseMean" & "log2FoldChange" columns

- fdr.thres:

  False discovery rate (FDR) threshold

- fc.thres:

  log2FoldChange threshold

- y.axis.lim:

  y-axis limits for negative log10 p-values

- fc.lim:

  x-axis limits

- lab.genes:

  genes to label on plot

- tolower.cols:

  column names that will be converted to lower case

- color_by:

  variable to color points by ('baseMean' or 'significance')

- colorscale:

  ggplot-compatible continuous colorscale name (viridis option), only
  used when color_by='baseMean'

- alpha:

  point opacity between 0 and 1

- pvalue_padj_switch:

  Switches between having the y-axis display raw or adjusted p-values.
  Significance is always determined by adjusted values.

## Value

ggplot handle

## Examples

``` r
n_genes <- 100
res <- data.frame(
         baseMean = runif(n_genes, 10, 1000),
         log2FoldChange = rnorm(n_genes, 0, 2),
         lfcSE = runif(n_genes, 0.1, 0.5),
         stat = rnorm(n_genes, 0, 3),
         pvalue = runif(n_genes, 0, 1),
         padj = runif(n_genes, 0, 1),
         symbol = paste0("GENE", 1:n_genes),
         row.names = paste0("gene", 1:n_genes)
       )

plotVolcano.label(res, lab.genes = c("gene1", "gene2"))
#> Error in plotVolcano.label(res, lab.genes = c("gene1", "gene2")): could not find function "plotVolcano.label"
```
