# Download button module

Module UI & server for download buttons.

## Usage

``` r
downloadButtonUI(id)

downloadButtonServer(id, outplot, plot_type)
```

## Arguments

- id:

  Module id

- outplot:

  reactive plot handle

- plot_type:

  reactive/static value used for output filename

## Value

UI returns tagList with download button UI. Server invisibly returns
NULL (used for side effects).

## Examples

``` r
library(shiny)
library(ggplot2)

# get example object
data(res_dex, package='carnation')
res <- as.data.frame(res_dex)

# make MA plot
p <- ggplot(res, aes(x=baseMean, y=log2foldChange)) +
       geom_point(color='black', alpha=0.5)

outplot <- reactive({ p })

# app with a single button to download a plot
if(interactive()){
  shinyApp(
    ui = fluidPage(
           downloadButtonUI('p')
         ),
    server = function(input, output, session){
               downloadButtonServer('p', outplot, 'maplot')
             }
  )
}
```
