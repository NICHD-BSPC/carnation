# Pattern plot module

Module UI & server to generate pattern plots.

## Usage

``` r
patternPlotUI(id, panel, tab)

patternPlotServer(id, obj, coldata, plot_args, config)
```

## Arguments

- id:

  Module id

- panel:

  string, can be 'sidebar' or 'main'

- tab:

  string, if 'plot' show plot settings, if 'table' show table settings;
  if 'both', show settings for both.

- obj:

  reactiveValues object containing carnation object

- coldata:

  reactiveValues object containing object metadata

- plot_args:

  reactive containing 'gene_scratchpad' (genes selected in scratchpad) &
  'upset_data' (list containing data from upset plot module)

- config:

  reactive list with config settings

## Value

UI returns tagList with module UI server invisibly returns NULL (used
for side effects)

## Examples

``` r
if (FALSE) { # interactive()
library(shiny)
library(DESeq2)

# Create reactive values to simulate app state
oobj <- make_example_carnation_object()

obj <- reactiveValues(
   dds = oobj$dds,
   rld = oobj$rld,
   res = oobj$res,
   all_dds = oobj$all_dds,
   all_rld = oobj$all_rld,
   dds_mapping = oobj$dds_mapping
)

cdata <- lapply(oobj$rld, function(x) colData(x))

coldata <- reactiveValues( all=cdata, curr=cdata )

plot_args <- reactive({
  list(
    gene_scratchpad=c('gene1', 'gene2'),
    upset_data=list(genes=NULL, labels=NULL)
  )
})

config <- reactiveVal(get_config())

shinyApp(
  ui = fluidPage(
         sidebarPanel(
           patternPlotUI('p', 'sidebar', 'both'),
           conditionalPanel(condition = "input.pattern_mode == 'Plot'",
             patternPlotUI('p', 'sidebar', 'plot')
           ),
           conditionalPanel(condition = "input.pattern_mode == 'Table'",
             patternPlotUI('p', 'sidebar', 'table')
           )
         ),
         mainPanel(
           tabsetPanel(id='pattern_mode',
             tabPanel('Plot',
               patternPlotUI('p', 'plot')
             ), # tabPanel plot

             tabPanel('Cluster membership',
               patternPlotUI('p', 'table')
             ) # tabPanel cluster_membership

           ) # tabsetPanel pattern_mode
         ) # tabPanel pattern_analysis
       ),
  server = function(input, output, session){
             patternPlotServer('deg_plot', obj, coldata,
                               plot_args, config)
           }
)
}
```
