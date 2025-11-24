# Upset plot module server function

Upset plot module server function

## Usage

``` r
upsetPlotServer(id, obj, plot_args, gene_scratchpad, reset_genes, config)
```

## Arguments

- id:

  ID string used to match the ID used to call the module UI function

- obj:

  reactiveValues object containing carnation object

- plot_args:

  reactive containing 'fdr.thres' (padj threshold) & 'fc.thres' (log2FC)

- gene_scratchpad:

  reactiveValues object containing genes selected in scratchpad

- reset_genes:

  reactive to reset gene scratchpad selection

- config:

  reactive list with config settings
