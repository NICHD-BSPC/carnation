# Plot an interactive PCA plot

Plot an interactive PCA plot

## Usage

``` r
plotPCA.ly(rld, intgroup)
```

## Arguments

- rld:

  DESeqTransform object output by varianceStabilizingTransformation() or
  rlog()

- intgroup:

  character vector of names in colData(x) to use for grouping

## Value

Handle to ggplot with added label field in aes_string() for plotting
with ggplotly()
