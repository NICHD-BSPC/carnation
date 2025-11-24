# Gene plot module server function

Gene plot module server function

## Usage

``` r
genePlotServer(id, obj, coldata, plot_args, config)
```

## Arguments

- id:

  ID string used to match the ID used to call the module UI function

- obj:

  reactiveValues object containing carnation object

- coldata:

  reactiveValues object containing object metadata

- plot_args:

  reactive list with 3 elements: 'gene.id' (all gene IDs) &
  'gene_scratchpad' (genes selected in scratchpad) & 'comp_all'
  (selected comparison)

- config:

  reactive list with config settings
