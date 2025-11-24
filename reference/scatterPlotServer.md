# Scatterplot module server function

Scatterplot module server function

## Usage

``` r
scatterPlotServer(id, obj, plot_args, config)
```

## Arguments

- id:

  ID string used to match the ID used to call the module UI function

- obj:

  reactiveValues object containing carnation object

- plot_args:

  reactive containing 'fdr.thres' (padj threshold), 'fc.thres' (log2FC)
  & 'gene.to.plot' (genes to be labeled)

- config:

  reactive list with config settings
