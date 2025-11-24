# MA plot module server function

MA plot module server function

## Usage

``` r
maPlotServer(id, obj, plot_args, config)
```

## Arguments

- id:

  ID string used to match the ID used to call the module UI function

- obj:

  reactiveValues object containing carnation object

- plot_args:

  reactive containing 'fdr.thres' (padj threshold), 'fc.thres' (log2FC
  threshold) & 'gene.to.plot' (genes selected in scratchpad)

- config:

  reactive list with config settings
