# Metadata module

This module generates the metadata tab that allows users to view the
metadata associated with the loaded carnation object.

## Usage

``` r
metadataUI(id, panel)

metadataServer(id, obj, cols.to.drop)
```

## Arguments

- id:

  Module id

- panel:

  context for generating ui elements ('sidebar' or 'main')

- obj:

  reactiveValues object containing carnation object

- cols.to.drop:

  columns to hide from table

## Value

UI returns tagList with metadata UI. Server returns reactive object with
metadata.

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

config <- get_config()
cols.to.drop <- config$server$cols.to.drop

shinyApp(
  ui = fluidPage(
         sidebarPanel(metadataUI('p', 'sidebar')),
         mainPanel(metadataUI('p', 'main'))
       ),
  server = function(input, output, session){
             # reactiveVal to save updates
             saved_data <- reactiveVal()

             cdata <- metadataServer('p', obj, cols.to.drop)

             observeEvent(cdata(), {
               saved_data(cdata())
             })
           }
)
}
```
