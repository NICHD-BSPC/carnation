# Cnetplot module

UI & module to generate Cnetplots.

## Usage

``` r
cnetPlotUI(id, panel)

cnetPlotServer(id, obj, config)
```

## Arguments

- id:

  Module id

- panel:

  string, can be 'sidebar' or 'main'

- obj:

  reactive, dataframe containing enrichment results

- config:

  reactive list with config settings

## Value

UI returns tagList with plot UI server invisibly returns NULL (used for
side effects)

## Examples

``` r
library(shiny)

# get DESeqResults object
data(res_dex, package='carnation')

obj <- reactive({ res })

config <- reactiveVal(get_config())

# run simple shiny app with plot
if(interactive()){
  shinyApp(
    ui = fluidPage(
           sidebarPanel(cnetPlotUI('p', 'sidebar')),
           mainPanel(cnetPlotUI('p', 'main'))
         ),
    server = function(input, output, session){
               cnetPlotServer('p', obj, config)
             }
  )
}
```
