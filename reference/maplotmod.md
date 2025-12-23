# MA plot module

UI & server for module to create MA plot

## Usage

``` r
maPlotUI(id, panel)

maPlotServer(id, obj, plot_args, config)
```

## Arguments

- id:

  Module id

- panel:

  string, can be 'sidebar' or 'main'

- obj:

  reactiveValues object containing carnation object

- plot_args:

  reactive containing 'fdr.thres' (padj threshold), 'fc.thres' (log2FC
  threshold) & 'gene.to.plot' (genes selected in scratchpad)

- config:

  reactive list with config settings

## Value

UI returns tagList with MA plot UI. Server invisibly returns NULL (used
for side effects).

## Examples

``` r
library(shiny)
library(DESeq2)

# Create reactive values to simulate app state
oobj <- make_example_carnation_object()
#> estimating size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing

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

if(interactive()){
  shinyApp(
    ui = fluidPage(
           sidebarPanel(maPlotUI('p', 'sidebar')),
           mainPanel(maPlotUI('p', 'main'))
         ),
    server = function(input, output, session){
               maPlotServer('p', obj, plot_args, config)
             }
  )
}
```
