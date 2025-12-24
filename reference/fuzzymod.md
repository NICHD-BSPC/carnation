# Fuzzy enrichment map module

UI & module to generate fuzzy enrichment map plots.

## Usage

``` r
fuzzyPlotUI(id, panel)

fuzzyPlotServer(id, obj, args, config)
```

## Arguments

- id:

  Module id

- panel:

  string, can be 'sidebar' or 'main'

- obj:

  reactive containing 'distilled' enrichment results

- args:

  reactive, list with plot arguments, 'numcat' (number of categories to
  plot)

- config:

  reactive list with config settings

## Value

UI returns tagList with plot UI server returns reactive with number of
plotted terms

## Examples

``` r
library(shiny)

# get enrichResult object
data(eres_dex, package='carnation')

# preprocess & convert to GeneTonic object
gt <- GeneTonic::shake_enrichResult(eres_dex)
#> Found 2483 gene sets in `enrichResult` object, of which 2483 are significant.
#> Converting for usage in GeneTonic...

# get distilled results
df <- GeneTonic::gs_fuzzyclustering(gt[seq_len(10),],
        similarity_threshold = 0.35,
        fuzzy_seeding_initial_neighbors = 3,
        fuzzy_multilinkage_rule = 0.5)

# number of plotted terms
args <- reactive({ list(numcat=10) })

config <- reactiveVal(get_config())

# run simple shiny app with plot
if(interactive()){
  shinyApp(
    ui = fluidPage(
           sidebarPanel(fuzzyPlotUI('p', 'sidebar')),
           mainPanel(fuzzyPlotUI('p', 'main'))
         ),
    server = function(input, output, session){
               numcat <- observe({
                 fuzzyPlotServer('p',
                                 reactive({ df }),
                                 args,
                                 config)
               })
             }
  )
}
```
