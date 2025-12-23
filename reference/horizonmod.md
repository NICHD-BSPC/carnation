# Horizon plot module

UI & module to generate horizon plots.

## Usage

``` r
horizonUI(id, panel)

horizonServer(id, obj, config)
```

## Arguments

- id:

  Module id

- panel:

  string, can be 'sidebar' or 'main'

- obj:

  reactiveValues object containing two GeneTonic objects

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

# get second enrichResult object
data(eres_cell, package='carnation')

# convert to GeneTonic object
gt1 <- GeneTonic::shake_enrichResult(eres_cell)
#> Found 2020 gene sets in `enrichResult` object, of which 2020 are significant.
#> Converting for usage in GeneTonic...

obj <- reactive({
  list(
    obj1 = list(l_gs = gt$l_gs,
             anno_df = gt$anno_df,
             label = 'comp1'),
    obj2 = list(l_gs = gt1$l_gs,
             anno_df = gt1$anno_df,
             label = 'comp2')
  )
})

config <- reactiveVal(get_config())

# run simple shiny app with plot
if(interactive()){
  shinyApp(
    ui = fluidPage(
           sidebarPanel(horizonUI('p', 'sidebar')),
           mainPanel(horizonUI('p', 'main'))
         ),
    server = function(input, output, session){
               horizonServer('p', obj, config)
             }
  )
}
```
