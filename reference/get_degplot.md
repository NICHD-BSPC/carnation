# Plot a degPatterns object

This function plots a degPatterns object.

## Usage

``` r
get_degplot(
  obj,
  time,
  color = NULL,
  cluster_column = "cluster",
  cluster_to_show,
  x_order,
  points = TRUE,
  boxes = TRUE,
  smooth = "smooth",
  lines = TRUE,
  facet = TRUE,
  prefix_title = "Cluster ",
  genes_to_label = NULL
)
```

## Arguments

- obj:

  degPatterns object

- time:

  metadata variable to plot on x-axis

- color:

  variable to color plot

- cluster_column:

  column to use for grouping genes

- cluster_to_show:

  which clusters to show in plot

- x_order:

  order of x-axis values

- points:

  boolean, show samples on plot? Default: TRUE

- boxes:

  boolean, show boxes on plot? Default: TRUE

- smooth:

  what type of trendline to use? can be 'smooth' (default) or 'line'.

- lines:

  show lines joining samples? Default: TRUE

- facet:

  boolean, should plot be faceted? Default: TRUE

- prefix_title:

  string, prefix for facet titles

- genes_to_label:

  genes to label on plot

## Value

ggplot handle

## Examples

``` r
# get degpatterns object
data(degpatterns_dex, package = 'carnation')

# get pattern plot
all_clusters <- unique(degpatterns_dex$normalized$cluster)

dp <- get_degplot(degpatterns_dex, time='dex',
                  cluster_to_show=all_clusters,
                  x_order=c('untrt','trt'))
```
