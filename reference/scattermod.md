# Scatterplot module

Module UI + server for generating scatter plots.

## Usage

``` r
scatterPlotUI(id, panel)

scatterPlotServer(id, obj, plot_args, config)
```

## Arguments

- id:

  Module id

- panel:

  string, can be 'sidebar' or 'main' passed to UI

- obj:

  reactiveValues object containing carnation object passed to server

- plot_args:

  reactive containing 'fdr.thres' (padj threshold), 'fc.thres' (log2FC)
  & 'gene.to.plot' (genes to be labeled) passed to server

- config:

  reactive list with config settings passed to server

## Value

UI returns tagList with scatter plot UI. Server invisibly returns NULL
(used for side effects).

## Examples

``` r
if (FALSE) { # interactive()
library(shiny)

# Create reactive values to simulate app state
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
    fc.thres=0,
    gene.to.plot=c('gene1', 'gene2')
  )
})

config <- reactiveVal(get_config())

shinyApp(
  ui = fluidPage(
         sidebarPanel(scatterPlotUI('p', 'sidebar')),
         mainPanel(scatterPlotUI('p', 'sidebar'))
       ),
  server = function(input, output, session){
             scatterPlotServer('p', obj, plot_args, config)
           }
)
}
```
