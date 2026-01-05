# Upset plot module

Module UI & server to generate upset plots.

## Usage

``` r
upsetPlotUI(id, panel)

upsetPlotServer(id, obj, plot_args, gene_scratchpad, reset_genes, config)
```

## Arguments

- id:

  Module id

- panel:

  string, can be 'sidebar' or 'main'

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

## Value

UI returns tagList with upset plot UI. Server returns reactive with list
containing upset table, intersections & selected genes.

## Examples

``` r
if (FALSE) { # interactive()
library(shiny)

oobj <- make_example_carnation_object()

obj <- reactiveValues(
   dds = oobj$dds,
   rld = oobj$rld,
   res = oobj$res,
   all_dds = oobj$all_dds,
   all_rld = oobj$all_rld,
   dds_mapping = oobj$dds_mapping
)

plot_args <- reactive({
  list(
    fdr.thres=0.1,
    fc.thres=0
  )
})

gene_scratchpad <- reactive({ c('gene1', 'gene2') })

reset_genes <- reactiveVal()

config <- reactiveVal(get_config())

shinyApp(
  ui = fluidPage(
         sidebarPanel(upsetPlotUI('p', 'sidebar')),
         mainPanel(upsetPlotUI('p', 'sidebar'))
       ),
  server = function(input, output, session){
             upset_data <- upsetPlotServer('p', obj, plot_args,
                                           gene_scratchpad,
                                           reset_genes, config)
           }
)
}
```
