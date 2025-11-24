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

  string, what values to plot? can be 'LFC' or 'P-adj'

- df:

  data frame with LFC & padj values to plot from 2 contrasts

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
