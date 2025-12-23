# Settings module

Module UI & server for user access details interface.

Server code for settings module

## Usage

``` r
settingsUI(id, panel, username)

settingsServer(id, details, depth, end_offset, assay_fun, config)
```

## Arguments

- id:

  Module id

- panel:

  context for generating ui elements ('sidebar' or 'main')

- username:

  user name

- details:

  reactive list with user name & app location details

- depth:

  project name depth

- end_offset:

  project name end offset

- assay_fun:

  function to parse assay names from file path

- config:

  reactive list with config settings

## Value

UI returns tagList with module UI Server returns reactive with list
containing user access details

## Examples

``` r
library(shiny)

# default username
username <- reactive({ NULL })

# internal carnation config
config <- reactiveVal(get_config())

# regex to find carnation files
pattern <- reactive({ config()$server$pattern })

# access permissions
assay.list <- reactiveValues(l=read_access_yaml())
#> Environment variable "CARNATION_ACCESS_YAML" not found.Using default location to save access yaml:/home/runner

if(interactive()){
  shinyApp(
    ui = fluidPage(
           sidebarPanel(uiOutput('settings_sidebar')),
           mainPanel(uiOutput('settings_main'))
         ),
    server = function(input, output, session){
               output$settings_main <- renderUI({
                 settingsUI('settings', panel='main', username=username)
               })

               output$settings_sidebar <- renderUI({
                 settingsUI('settings', panel='sidebar', username=username)
               })

               settings <- settingsServer('p',
                                          details=reactive({
                                                    list(username=username,
                                                         where=NULL)
                                                  }),
                                          depth=2,
                                          end_offset=0,
                                          assay_fun=function(x)
                                            sub(paste0(pattern(), '\\.rds$'), '',
                                                basename(x),
                                                ignore.case=TRUE),
                                          config
                                          )
             }
  )
}

```
