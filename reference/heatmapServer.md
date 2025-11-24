# Heatmap module server function

Heatmap module server function

## Usage

``` r
heatmapServer(id, obj, coldata, plot_args, gene_scratchpad, config)
```

## Arguments

- id:

  ID string used to match the ID used to call the module UI function

- obj:

  reactiveValues object containing carnation object

- coldata:

  reactiveValues object containing object metadata

- plot_args:

  reactive containing 'fdr.thres' (padj threshold), 'fc.thres' (log2FC)
  & 'upset_data' (list containing data from upset plot module)

- gene_scratchpad:

  reactiveValues object containing genes selected in scratchpad which
  will be labeled

- config:

  reactive list with config settings
