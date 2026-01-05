# Heatmap module

Module UI & server to generate heatmap.

## Usage

``` r
heatmapUI(id, panel)

heatmapServer(id, obj, coldata, plot_args, gene_scratchpad, config)
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

  reactive containing 'fdr.thres' (padj threshold), 'fc.thres' (log2FC)
  & 'upset_data' (list containing data from upset plot module)

- gene_scratchpad:

  reactiveValues object containing genes selected in scratchpad which
  will be labeled

- config:

  reactive list with config settings

## Value

UI returns tagList with heatmap UI. Server invisibly returns NULL (used
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

cdata <- lapply(oobj$rld, function(x) colData(x))

coldata <- reactiveValues( all=cdata, curr=cdata )

plot_args <- reactive({
  list(
    fdr.thres=0.1,
    fc.thres=0,
    upset_data=list(genes=NULL, labels=NULL)
  )
})

gene_scratchpad <- reactive({ c('gene1', 'gene2') })

config <- reactiveVal(get_config())

shinyApp(
  ui = fluidPage(
         sidebarPanel(heatmapUI('p', 'sidebar')),
         mainPanel(heatmapUI('p', 'sidebar'))
       ),
  server = function(input, output, session){
             heatmapServer('p', obj, coldata,
                           plot_args, gene_scratchpad, config)
           }
)
}
```
