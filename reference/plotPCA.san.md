# Adjustable PCA plot

Create a PCA plot with specified PCs on x- and y-axis

## Usage

``` r
plotPCA.san(
  object,
  intgroup = "group",
  pcx,
  pcy,
  pcz = NULL,
  ntop = 500,
  samples = NULL,
  loadings = FALSE,
  loadings_ngenes = 10
)
```

## Arguments

- object:

  normalized DESeqDataSet object

- intgroup:

  metadata variable to use for grouping samples

- pcx:

  principal component to plot on x-axis

- pcy:

  principal component to plot on y-axis

- pcz:

  principal component to plot on z-axis. If not NULL, function returns a
  3-D PCA plot.

- ntop:

  number of most-variable genes to use

- samples:

  vector of sample names to show on plot

- loadings:

  boolean, show gene loadings? Default is FALSE.

- loadings_ngenes:

  integer, \# genes to show loadings for (default=10)

## Value

ggplot handle

## Examples

``` r
# make example dds object
dds <- DESeq2::makeExampleDESeqDataSet()

# normalize
rld <- DESeq2::varianceStabilizingTransformation(dds, blind=TRUE)

# make pca plot
p <- plotPCA.san(rld, intgroup='condition', pcx='PC1', pcy='PC2')
#> Warning: NAs introduced by coercion
#> Warning: NAs introduced by coercion

```
