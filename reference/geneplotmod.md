# Gene plot module

UI & server for module to create gene plot

## Usage

``` r
genePlotUI(id, panel)

genePlotServer(id, obj, coldata, plot_args, config)
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

- plot_args:

  reactive list with 3 elements: 'gene.id' (all gene IDs) &
  'gene_scratchpad' (genes selected in scratchpad) & 'comp_all'
  (selected comparison)

- config:

  reactive list with config settings

## Value

UI returns tagList with gene plot UI. Server invisibly returns NULL
(used for side effects).

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
    gene.to.plot = c("gene1", "gene2"),
    gene.id = rownames(oobj$dds$main),
    comp_all = "comp1"
  )
})

config <- reactiveVal(get_config())

if(interactive()){
  shinyApp(
    ui = fluidPage(
           sidebarPanel(genePlotUI('p', 'sidebar')),
           mainPanel(genePlotUI('p', 'main'))
         ),
    server = function(input, output, session){
               genePlotServer('p', obj, coldata, plot_args, config)
             }
  )
}
```
