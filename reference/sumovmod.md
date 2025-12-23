# Summary overview plot module

UI & module to generate summary overview plots.

## Usage

``` r
sumovPlotUI(id, panel, type = "")

sumovPlotServer(id, obj, config, type = "")
```

## Arguments

- id:

  Module id

- panel:

  string, can be 'sidebar' or 'main'

- type:

  string, if 'comp' then show the comparison view

- obj:

  reactiveValues object containing GeneTonic object

- config:

  reactive list with config settings

## Value

UI returns tagList with plot UI server invisibly returns NULL (used for
side effects)

## Examples

``` r
library(shiny)

# get enrichResult object
data(eres_dex, package='carnation')

# convert to GeneTonic object
gt <- GeneTonic::shake_enrichResult(eres_dex)
#> Found 2186 gene sets in `enrichResult` object, of which 2186 are significant.
#> Converting for usage in GeneTonic...

obj <- reactive({
  list(l_gs = gt$l_gs,
       anno_df = gt$anno_df,
       label = 'comp1')
})

config <- reactiveVal(get_config())

# run simple shiny app with plot
if(interactive()){
  shinyApp(
    ui = fluidPage(
           sidebarPanel(sumovPlotUI('p', 'sidebar')),
           mainPanel(sumovPlotUI('p', 'main'))
         ),
    server = function(input, output, session){
               sumovPlotServer('p', obj, config)
             }
  )
}
```
