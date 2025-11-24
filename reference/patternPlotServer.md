# Pattern plot module server function

Pattern plot module server function

## Usage

``` r
patternPlotServer(id, obj, coldata, plot_args, config)
```

## Arguments

- id:

  ID string used to match the ID used to call the module UI function

- obj:

  reactiveValues object containing carnation object

- coldata:

  reactiveValues object containing object metadata

- plot_args:

  reactive containing 'gene_scratchpad' (genes selected in scratchpad) &
  'upset_data' (list containing data from upset plot module)

- config:

  reactive list with config settings
