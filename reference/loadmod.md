# Load data module

Module UI & server to load new data

## Usage

``` r
loadDataUI(id)

loadDataServer(id, username, config, rds = NULL)
```

## Arguments

- id:

  Module id

- username:

  user name

- config:

  reactive list with config settings

- rds:

  Object to be edited

## Value

UI returns tagList with module UI Server returns reactive with app
reload trigger

## Examples

``` r
library(shiny)

username <- 'admin'

config <- reactiveVal(get_config())

obj <- make_example_carnation_object()
#> estimating size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing

rds <- reactive({ obj=obj })

if(interactive()){
  shinyApp(
    ui = fluidPage(
           loadDataUI('p')
         ),
    server = function(input, output, session){
               loadDataServer('p', username=username, config, rds)
             }
  )
}

```
