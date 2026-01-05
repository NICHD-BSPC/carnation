# Save object module UI

Module UI & server to save carnation object.

## Usage

``` r
saveUI(id)

saveServer(id, original, current, coldata, pattern, username, config)
```

## Arguments

- id:

  Module id

- original:

  original carnation object

- current:

  current carnation object

- coldata:

  reactiveValues object containing object metadata

- pattern:

  regex pattern for finding carnation data

- username:

  user name

- config:

  reactive list with config settings

## Value

UI returns actionButton Server returns reactive with trigger to refresh
the app

## Examples

``` r
if (FALSE) { # interactive()
library(shiny)
library(DESeq2)

# default username
username <- reactive({ NULL })

# internal carnation config
config <- reactiveVal(get_config())

# regex to find carnation files
pattern <- reactive({ config()$server$pattern })

# get example object
obj <- make_example_carnation_object()

# make reactive with obj & path
original <- reactiveValues( obj = obj, path = "/path/to/carnation/obj.rds" )

# extract metadata
coldata <- reactive({ lapply(obj$dds, colData) })

# edit metadata
coldata_edit <- lapply(coldata, function(x){
                  x$type <- 'new'; x
                })

# add to object
edit_obj <- obj
for(name in names(edit_obj$dds)){
  colData(edit_obj$dds[[ name ]]) <- coldata_edit[[ name ]]
}

# run simple shiny app with plot
shinyApp(
  ui = fluidPage(
         saveUI('p')
       ),
  server = function(input, output, session){
             save_event <- saveServer('save_object',
                                      original=original,
                                      current=reactive({ edit_obj }),
                                      coldata=coldata,
                                      pattern=pattern(),
                                      username=username,
                                      config)
           }
)
}
```
