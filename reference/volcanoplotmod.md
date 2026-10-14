# Volcano plot module

UI & server for module to create volcano plot

## Usage

``` r
volcanoPlotUI(id, panel)

volcanoPlotServer(id, obj, plot_args, config)
```

## Arguments

- id:

  Module id

- panel:

  string, can be 'sidebar' or 'main'

- obj:

  reactiveValues object containing carnation object

- plot_args:

  reactive containing 'fdr.thres' (FDR threshold), 'fc.thres' (log2FC
  threshold) & 'gene.to.plot' (genes selected in scratchpad)

- config:

  reactive list with config settings

## Value

UI returns tagList with volcano plot UI. Server invisibly returns NULL
(used for side effects).

## Examples

``` r
if (FALSE) { # interactive()
library(shiny)
library(DESeq2)

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

# Set up coldata structure that the module expects
coldata <- reactiveValues(
  curr = list(
    all_samples = colData(oobj$dds$main),
    main = colData(oobj$dds$main)
  )
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
         sidebarPanel(volcanoPlotUI('p', 'sidebar')),
         mainPanel(volcanoPlotUI('p', 'main'))
       ),
  server = function(input, output, session){
             volcanoPlotServer('p', obj, plot_args, config)
           }
)
}
```
