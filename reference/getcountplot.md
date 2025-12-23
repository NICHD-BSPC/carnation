# Create gene plot

This function creates the gene plot.

## Usage

``` r
getcountplot(
  df,
  intgroup = "group",
  factor.levels,
  title = NULL,
  ylab = "Normalized counts",
  color = "gene",
  nrow = 2,
  ymin = NULL,
  ymax = NULL,
  log = TRUE,
  freey = FALSE,
  trendline = "smooth",
  facet = NULL,
  legend = TRUE,
  boxes = TRUE,
  rotate_x_labels = 30
)
```

## Arguments

- df:

  data.frame with gene counts

- intgroup:

  metadata variable to plot on x-axis

- factor.levels:

  levels of intgroup to show on x-axis

- title:

  title of plot

- ylab:

  y-axis label

- color:

  metadata variable to color by

- nrow:

  number of rows to plot if faceting

- ymin:

  y-axis lower limit

- ymax:

  y-axis upper limit

- log:

  should y-axis be log10-transformed?

- freey:

  should y-axes of faceted plots have independent scales?

- trendline:

  type of trendline to draw

- facet:

  metadata variable to facet by

- legend:

  show legend?

- boxes:

  show boxes?

- rotate_x_labels:

  angle to rotate x-axis labels (default=30)

## Value

ggplot handle

## Examples

``` r
# make example DESeq dataset
dds <- DESeq2::makeExampleDESeqDataSet()

# get gene counts
df <- get_gene_counts(dds, gene = c('gene1', 'gene2'))

# standard gene plot
p <- getcountplot(df, intgroup = "condition", factor.levels = c("A", "B"))

# with genes faceted
p1 <- getcountplot(df, intgroup = "condition", factor.levels = c("A", "B"), facet = "gene")

```
