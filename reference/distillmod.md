# Distilled enrichment map module

UI & module to generate distill enrichment map plots.

## Usage

``` r
distillPlotUI(id, panel)

distillPlotServer(id, obj, args, config)
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
library(GeneTonic)
#> Welcome to GeneTonic v3.0.0
#> 
#> If you use GeneTonic in your work, please cite:
#> 
#>   GeneTonic: an R/Bioconductor package for streamlining the interpretation of RNA-seq data
#>   Federico Marini, Annekathrin Ludt, Jan Linke, Konstantin Strauch
#>   BMC Bioinformatics, 2021 - https://doi.org/10.1186/s12859-021-04461-5
#> and/or (if adopting the series of protocols as a whole)
#>   Interactive and Reproducible Workflows for Exploring and Modeling RNA-seq Data with pcaExplorer, ideal, and GeneTonic
#>   Annekathrin Ludt, Arsenij Ustjanzew, Harald Binder, Konstantin Strauch, Federico Marini
#>   Current Protocols, 2022 - https://doi.org/10.1002/cpz1.411
#> 
#> Attaching package: ‘GeneTonic’
#> The following object is masked from ‘package:carnation’:
#> 
#>     gs_radar
library(shiny)

# get DESeqResults object
data(res_dex, package='carnation')

# get enrichResult object
data(eres_dex, package='carnation')

# preprocess & convert to GeneTonic object
eres2 <- GeneTonic::shake_enrichResult(eres_dex)
#> Found 2483 gene sets in `enrichResult` object, of which 2483 are significant.
#> Converting for usage in GeneTonic...
gt <- enrich_to_genetonic(eres_dex, res_dex)
#> Found 2483 gene sets in `enrichResult` object, of which 2483 are significant.
#> Converting for usage in GeneTonic...

# get distilled results
df <- distill_enrichment(
        eres2,
        res_dex,
        gt$anno_df,
        n_gs = 10,
        cluster_fun = "cluster_markov"
      )

# number of plotted terms
args <- reactive({ list(numcat=10) })

config <- reactiveVal(get_config())

# run simple shiny app with plot
if(interactive()){
  shinyApp(
    ui = fluidPage(
           sidebarPanel(distillPlotUI('p', 'sidebar')),
           mainPanel(distillPlotUI('p', 'main'))
         ),
    server = function(input, output, session){
               numcat <- observe({
                 distillPlotServer('p',
                                   reactive({ df }),
                                   args,
                                   config)
               })
             }
  )
}
```
