# Alluvial plot module

UI & module to generate alluvial plots.

## Usage

``` r
alluvialUI(id, panel)

alluvialServer(id, obj, res_obj, config)
```

## Arguments

- id:

  Module id

- panel:

  string, can be 'sidebar' or 'main'

- obj:

  reactiveValues object containing GeneTonic object

- res_obj:

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

# get enrichResult object
data(eres_dex, package='carnation')

# convert to GeneTonic object

gt <- GeneTonic::shake_enrichResult(eres_dex)
#> Found 2483 gene sets in `enrichResult` object, of which 2483 are significant.
#> Converting for usage in GeneTonic...

obj <- reactive({
  list(l_gs = gt$l_gs,
       anno_df = gt$anno_df,
       label = 'comp1')
})

res_obj <- reactive({ res })

config <- reactiveVal(get_config())

# run simple shiny app with plot
if(interactive()){
  shinyApp(
    ui = fluidPage(
           sidebarPanel(alluvialUI('p', 'sidebar')),
           mainPanel(alluvialUI('p', 'main'))
         ),
    server = function(input, output, session){
               alluvialServer('p', obj, res_obj, config)
             }
  )
}
```
