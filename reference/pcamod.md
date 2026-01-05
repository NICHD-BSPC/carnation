# PCA plot module

Module UI + server to generate a pca plot.

## Usage

``` r
pcaPlotUI(id, panel)

pcaPlotServer(id, obj, coldata, config)
```

## Arguments

- id:

  Module id

- panel:

  string, can be 'sidebar' or 'main'

- obj:

  reactiveValues object containing carnation object

- coldata:

  reactiveValues object containing object metadata

- config:

  reactive list with config settings

## Value

UI returns tagList with PCA plot UI. Server invisibly returns NULL (used
for side effects).

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

config <- reactiveVal(get_config())

shinyApp(
  ui = fluidPage(
         sidebarPanel(pcaPlotUI('p', 'sidebar')),
         mainPanel(pcaPlotUI('p', 'main'))
       ),
  server = function(input, output, session){
             pcaPlotServer('p', obj, coldata, config)
           }
)
}
```
