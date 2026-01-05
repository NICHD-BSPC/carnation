# Plot a scatterplot to compare two contrasts

Plot a scatterplot to compare two contrasts

## Usage

``` r
plotScatter.label(
  compare,
  df,
  label_x,
  label_y,
  lim.x,
  lim.y,
  color.palette,
  lab.genes = NULL,
  plot_all = "no",
  name.col = "geneid",
  lines = c("yes", "yes", "yes"),
  alpha = 1,
  size = 4,
  show.grid = "yes"
)
```

## Arguments

- compare:

  string, what values to plot? can be 'log2FoldChange' or 'P-adj'

- df:

  data frame with log2FoldChange & padj values to plot from 2 contrasts

- label_x:

  string, label for x-axis

- label_y:

  string, label for y-axis

- lim.x:

  x-axis limits

- lim.y:

  y-axis limits

- color.palette:

  character vector of colors to use for significance categories 'Both -
  same LFC sign', 'Both - opposite LFC sign', 'None', label_x, label_y

- lab.genes:

  genes to label (default=NULL)

- plot_all:

  string, can be 'yes' or 'no'. if 'yes', points outside axis limits are
  plotted along x/y axis lines (default='no').

- name.col:

  gene name column to merge the 2 results, also used for labeling points

- lines:

  3-element character vector to plot gridlines in the order (x=0, y=0,
  x=y), with 'yes' or 'no' values. E.g. ('yes', 'yes', 'no') will plot
  dotted lines for x = 0 & y = 0, but not the x = y diagonal.

- alpha:

  float, marker opacity (default=1).

- size:

  float, marker size (default=4).

- show.grid:

  string, can be 'yes' (default) or 'no'.

## Value

ggplot handle

## Examples

``` r
library(DESeq2)

# make example dataset
dds <- makeExampleDESeqDataSet()

# run DE analysis
dds <- DESeq(dds)
#> estimating size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing

# extract comparisons of interest
res1 <- results(dds, contrast = c("condition", "A", "B"))
res2 <- results(dds, contrast = c("condition", "B", "A"))

# add geneid column
res1 <- cbind(geneid=rownames(res1), res1)
res2 <- cbind(geneid=rownames(res2), res2)

# make merged df from the two comparisons
cols.sub <- c('log2FoldChange', 'padj', 'geneid')
df_full <- dplyr::inner_join(
  dplyr::select(as.data.frame(res1), all_of(cols.sub)),
  dplyr::select(as.data.frame(res2), all_of(cols.sub)),
  by = 'geneid',
  suffix = c('.x', '.y')
)

# calculate x & y limits for log2FoldChange
xlim <- range(df_full[[ 'log2FoldChange.x' ]])
ylim <- range(df_full[[ 'log2FoldChange.y' ]])

# get color palette
color.palette <- RColorBrewer::brewer.pal(n=5, name='Set2')

# add significance column
sig.x <- df_full$padj.x < 0.1 & !is.na(df_full$padj.x)
sig.y <- df_full$padj.y < 0.1 & !is.na(df_full$padj.y)
up.x <- df_full$log2FoldChange.x >= 0
up.y <- df_full$log2FoldChange.y >= 0
significance <- rep('None', nrow(df_full))
significance[ sig.x & sig.y & ((up.x & up.y) | (!up.x & !up.y)) ] <- 'Both - same LFC sign'
significance[ sig.x & sig.y & ((up.x & !up.y) | (!up.x & up.y)) ] <- 'Both - opposite LFC sign'
significance[ sig.x & !sig.y ] <- 'A vs B'
significance[ !sig.x & sig.y ] <- 'B vs A'
df_full$significance <- significance

# generate scatter plot
p <- plotScatter.label(compare = 'log2FoldChange',
                       df = df_full,
                       label_x = 'A vs B',
                       label_y = 'B vs A',
                       lim.x = xlim,
                       lim.y = ylim,
                       color.palette = color.palette)
```
