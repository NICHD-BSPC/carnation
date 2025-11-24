# Functional enrichment module server function

Functional enrichment module server function

## Usage

``` r
enrichServer(id, obj, upset_table, gene_scratchpad, reset_genes, config)
```

## Arguments

- id:

  ID string used to match the ID used to call the module UI function

- obj:

  reactiveValues object containing carnation object

- upset_table:

  reactive, data from upset plot module

- gene_scratchpad:

  reactive, genes selected in gene scratchpad

- reset_genes:

  reactive to reset genes in scratchpad

- config:

  reactive list with config settings
