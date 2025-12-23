# Help button module

Module UI & server for help buttons.

## Usage

``` r
helpButtonUI(id)

helpButtonServer(id, ...)
```

## Arguments

- id:

  Module id. This also doubles as prefixes for help text files.

- ...:

  other params passed to helpModal()

## Value

UI returns tagList with help button UI. Server invisibly returns NULL
(used for side effects).

## Examples

``` r
library(shiny)

# app with a single help button to show DE summary table details
if(interactive()){
  shinyApp(
    ui = fluidPage(
           helpButtonUI('de_summary_help')
         ),
    server = function(input, output, session){
               helpButtonServer('de_summary_help')
             }
  )
}
```
