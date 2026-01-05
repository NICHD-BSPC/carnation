# Functional enrichment module

UI & module to show functional enrichment tables & plots.

## Usage

``` r
enrichUI(id, panel, tab = "none")

enrichServer(id, obj, upset_table, gene_scratchpad, reset_genes, config)
```

## Arguments

- id:

  ID string used to match the ID used to call the module UI function

- panel:

  string, can be 'sidebar' or 'main'

- tab:

  string, if 'table' show table settings, if 'plots' show plot settings;
  if 'compare_results', show comparison settings.

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

## Value

UI returns tagList with plot UI server returns reactive with gene
selected from functional enrichment tables.

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

upset_table <- reactiveValues(tbl=NULL, intersections=NULL, set_labels=NULL)

gene_scratchpad <- reactive({ c('gene1', 'gene2') })

config <- reactiveVal(get_config())

shinyApp(
  ui = fluidPage(
         sidebarPanel(
           conditionalPanel(condition = "input.func == 'Table'",
             enrichUI('p', 'sidebar', 'table')
           ),
           conditionalPanel(condition = "input.func == 'Plot'",
             enrichUI('p', 'sidebar', 'plot')
           ),
           conditionalPanel(condition = "input.func == 'Compare results'",
             enrichUI('p', 'sidebar', 'compare_results')
           )
         ),
         mainPanel(
             tabsetPanel(id='func',
               tabPanel('Table',
                 enrichUI('p', 'main', 'table')
               ), # tabPanel table

               tabPanel('Plot',
                 enrichUI('p', 'main', 'plot')
               ), # tabPanel plot

               tabPanel('Compare results',
                 enrichUI('p', 'main', 'compare_results')
               ) # tabPanel compare_results

             ) # tabsetPanel func
           ) # tabPanel
         ),
  server = function(input, output, session){
             enrich_data <- enrichServer('p', obj,
                                         upset_table,
                                         gene_scratchpad,
                                         reactive({ FALSE }),
                                         config)
           }
)
}
```
