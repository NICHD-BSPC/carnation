#' PCA plot module
#'
#' @description
#' Module UI + server to generate a pca plot.
#'
#' @param id Module id
#' @param panel string, can be 'sidebar' or 'main'
#' @param obj reactiveValues object containing carnation object
#' @param coldata reactiveValues object containing object metadata
#' @param config reactive list with config settings
#'
#' @returns
#' UI returns tagList with PCA plot UI.
#' Server invisibly returns NULL (used for side effects).
#'
#' @examplesIf interactive()
#' library(shiny)
#' library(DESeq2)
#'
#' # Create reactive values to simulate app state
#' oobj <- make_example_carnation_object()
#'
#' obj <- reactiveValues(
#'    dds = oobj$dds,
#'    rld = oobj$rld,
#'    res = oobj$res,
#'    all_dds = oobj$all_dds,
#'    all_rld = oobj$all_rld,
#'    dds_mapping = oobj$dds_mapping
#' )
#'
#' # Set up coldata structure that the module expects
#' coldata <- reactiveValues(
#'   curr = list(
#'     all_samples = colData(oobj$dds$main),
#'     main = colData(oobj$dds$main)
#'   )
#' )
#'
#' config <- reactiveVal(get_config())
#'
#' shinyApp(
#'   ui = fluidPage(
#'          sidebarPanel(pcaPlotUI('p', 'sidebar')),
#'          mainPanel(pcaPlotUI('p', 'main'))
#'        ),
#'   server = function(input, output, session){
#'              pcaPlotServer('p', obj, coldata, config)
#'            }
#' )
#'
#' @rdname pcamod
#' @name pcamod
NULL

#' @rdname pcamod
#' @export
pcaPlotUI <- function(id, panel){
  ns <- NS(id)

  # load default config
  config <- get_config()
  dims <- config$server$de_analysis$pca_plot$dims
  pc_choices <- setNames(seq_len(dims), paste0('PC', seq_len(dims)))

  if(panel == 'sidebar'){
    tagList(
      fluidRow(
        column(6, align='left',
          tags$label(class='control-label',
                     'Sample group')
        ), # column
        column(6, align='right',
          helpButtonUI(ns('pca_samples_help'))
        ) # column
      ), # fluidRow

      selectInput(ns('pca_samples'), label=NULL,
                  choices=NULL,
                  selected=NULL
      ), # selectInput


      conditionalPanel(paste0('input["', ns('pca_samples'), '"] == "subsets"'),
        fluidRow(
          column(6, 'choose subset'),
          column(6,
            selectizeInput(ns('comp_pca'),
                           label=NULL,
                           choices=NULL,
                           selected=NULL
            ) # selectInput
          ) # column
        ) # fluidRow
      ), # conditionalPanel

      fluidRow(
        column(4, 'color by'),
        column(8,
          selectInput(ns('pca_color'), label=NULL,
                      choices=NULL, selected=NULL,
                      multiple=TRUE
          ) # selectInput
        ) # column
      ), # fluidRow

      fluidRow(
        column(4, 'x-axis'),
        column(8,
          selectInput(ns('pcx'),
                      label=NULL,
                      choices=pc_choices,
                      selected=pc_choices[1])
        ) # column
      ), # fluidRow

      fluidRow(
        column(4, 'y-axis'),
        column(8,
        selectInput(ns('pcy'),
                    label=NULL,
                    choices=pc_choices,
                    selected=pc_choices[2])
        ) # column
      ), # fluidRow

      fluidRow(
        column(4, 'z-axis'),
        column(8,
        selectInput(ns('pcz'),
                    label=NULL,
                    choices=c('none', pc_choices))
        ) # column
      ), # fluidRow

      fluidRow(
        column(4, '# of genes'),
        column(8,
        numericInput(ns('ntop'),
                     label=NULL,
                     value=500, min=0, max=1000, step=100)
        ) # column
      ), # fluidRow

      bsCollapse(
        bsCollapsePanel('subset samples',
          fluidRow(
            column(4, 'group by'),
            column(8,
              selectInput(ns('pca_cols'),
                          label=NULL,
                          choices=NULL,
                          selected=NULL)
            ) # column
          ), # fluidRow
          selectizeInput(ns('pca_col_levels'),
                         label=NULL,
                         choices=NULL,
                         selected=NULL,
                         multiple=TRUE),

          actionButton(ns('reset_pca_cols'),'Reset',
                       class='btn-primary')
        ), # bsCollapsePanel
      #), # bsCollapse

      #bsCollapse(
        bsCollapsePanel('gene loadings',
          fluidRow(
            column(4, 'show loadings?'),
            column(8,
              selectInput(ns('pca_loadings'),
                          label=NULL,
                          choices=c('yes', 'no'),
                          selected='no')
            ) # column
          ), # fluidRow

          fluidRow(
            column(4, '# of top genes'),
            column(8,
              numericInput(ns('pca_loadings_ngenes'),
                           label=NULL,
                           value=10, min=0, max=20)
            ) # column
          ) # fluidRow
        ) # bsCollapsePanel
      ), # bsCollapse

      fluidRow(align='center',
        actionButton(ns('plot_do'),
                     label='Refresh plot',
                     icon=icon('arrows-rotate'),
                     class='btn-primary',
                     style='margin-bottom: 10px;')
      ) # fluidRow
    ) # tagList
  } else if(panel == 'main'){
    tagList(
      fluidRow(
        column(6, align='left',
          helpButtonUI(ns('de_pca_help'))
        ), # column
        column(6, align='right',
          downloadButtonUI(ns('pcaplot_download'))
        ) # column
      ), # fluidRow

      withSpinner(
        plotlyOutput(ns('pcaplot'), height='500px')
      ) # withSpinner
    )
  }
} # pcaPlotUI

#' PCA plot module server function
#'
#' @rdname pcamod
#' @export
pcaPlotServer <- function(id, obj, coldata, config){

  moduleServer(
    id,

    function(input, output, session){

      ns <- NS(id)

      coldata_pca <- reactiveValues(all=NULL, current=NULL)
      plot_data <- reactiveValues(rld=NULL)

      data_loaded <- reactiveValues(flag=0)

      coldata.all <- reactive({
        list(curr=coldata$curr)
      })

      app_object <- reactive({
        list(dds=obj$dds,
             rld=obj$rld,
             all_dds=obj$all_dds,
             all_rld=obj$all_rld,
             dds_mapping=obj$dds_mapping)
      })

      reset_data <- function(){
        coldata_pca$all <- NULL
        coldata_pca$current <- NULL
        plot_data$rld <- NULL

        updateSelectInput(session,
                          'pca_samples',
                          choices=NULL,
                          selected=NULL)

        updateSelectInput(session,
                          'comp_pca',
                          choices=NULL,
                          selected=NULL)
      }

      # update min menus from reactive config
      observeEvent(config(), {
        dims <- config()$server$de_analysis$pca_plot$dims
        pc_choices <- setNames(seq_len(dims), paste0('PC', seq_len(dims)))

        updateSelectInput(session, 'pcx',
                          choices=pc_choices)
        updateSelectInput(session, 'pcy',
                          choices=pc_choices,
                          selected=pc_choices[2])
        updateSelectInput(session, 'pcz',
                          choices=c('none', pc_choices))

      })
 
      # get plot data & update menus at first load
      observeEvent(c(app_object()$dds_mapping, coldata.all()), {
        validate(
          need(!is.null(app_object()$rld) & !is.null(app_object()$dds_mapping),
               'Waiting for selection')
        )

        validate(
          need(all(unique(unlist(app_object()$dds_mapping)) %in% names(app_object()$rld)), 'Waiting for selection')
        )

        reset_data()

        # update names of comparisons for pca
        comp_choices <- unlist(unique(app_object()$dds_mapping))
        updateSelectizeInput(session,
                             'comp_pca',
                             choices=comp_choices,
                             selected=comp_choices[1])

        if(is.null(app_object()$all_dds)){
            if(length(app_object()$dds) == 1){
              choices <- 'all_samples'
            } else {
              choices <- 'subsets'
            }
        } else {
            if(length(app_object()$dds) > 1){
                choices <- c('all_samples', 'subsets')
            } else {
                choices <- 'all_samples'
            }
        }
        updateSelectInput(session,
                          'pca_samples',
                          choices=choices,
                          selected=choices[1])

        # update plot data & metadata
        if(choices[1] == 'all_samples'){
            cdata <- coldata.all()$curr[[ choices[1] ]]
            if(!is.null(app_object()$all_rld)){
                rld <- app_object()$all_rld
            } else if(length(app_object()$rld) == 1){
                rld <- app_object()$rld[1]
            }

            colData(rld) <- cdata
        } else if(choices[1] == 'subsets'){
            cdata <- coldata.all()$curr[[ comp_choices[1] ]]
            rld <- app_object()$rld[[ comp_choices[1] ]]
            colData(rld) <- cdata
        }

        plot_data$rld <- rld
        coldata_pca$all <- cdata
        coldata_pca$current <- cdata
      })

      # observer to get PCA plot data and update menus
      observeEvent(c(input$pca_samples, input$comp_pca, coldata.all()), {

        validate(
          need(!is.null(app_object()$rld) & !is.null(input$comp_pca) & input$comp_pca != '' & !is.null(app_object()$dds_mapping), 'Waiting for selection')
        )

        validate(
          need(all(unique(unlist(app_object()$dds_mapping)) %in% names(app_object()$rld)), 'Waiting for selection')
        )

        validate(
          need(input$comp_pca %in% app_object()$dds_mapping, 'Waiting for selection')
        )

        validate(
          need(input$pca_samples != '', 'Waiting for selection')
        )

        if(input$pca_samples == 'all_samples'){
            cdata <- coldata.all()$curr[[ input$pca_samples ]]
            if(!is.null(app_object()$all_rld)){
                rld <- app_object()$all_rld
            } else if(length(app_object()$rld) == 1){
                rld <- app_object()$rld[1]
            }

            colData(rld) <- cdata
        } else if(input$pca_samples == 'subsets'){
            cdata <- coldata.all()$curr[[ input$comp_pca ]]
            rld <- app_object()$rld[[ input$comp_pca ]]
            colData(rld) <- cdata
        }

        plot_data$rld <- rld
        coldata_pca$all <- cdata
        coldata_pca$current <- cdata

        column.names <- setdiff(colnames(cdata), config()$server$cols.to.drop)

        if(!is.null(input$pca_cols) & input$pca_cols %in% column.names)
          selected <- input$pca_cols
        else selected <- column.names[1]
        updateSelectInput(session, 'pca_cols',
                          choices=column.names,
                          selected=selected)

        if(is.null(input$pca_color)) selected <- column.names[1]
        else if (input$pca_color %in% column.names)
          selected <- input$pca_color
        else selected <- column.names[1]
        updateSelectInput(session, 'pca_color',
                          choices=column.names,
                          selected=selected)

        all_levels <- unique(cdata[, column.names[1]])
        validate(
            need(length(all_levels) > 0, 'Waiting for selection')
        )

        updateSelectizeInput(session,
                             'pca_col_levels',
                             choices=all_levels,
                             selected=all_levels)

        data_loaded$flag <- data_loaded$flag + 1
      }) # observeEvent

      # observer to update pca column levels
      observeEvent(c(input$pca_cols, coldata_pca$all, coldata_pca$current), {
        validate(
            need(!is.null(coldata_pca$current) & input$pca_cols != '' & input$pca_cols %in% colnames(coldata_pca$current), 'Waiting for selection')
        )

        all_levels <- unique(coldata_pca$all[, input$pca_cols])
        curr_levels <- unique(coldata_pca$current[, input$pca_cols])

        updateSelectizeInput(session,
                             'pca_col_levels',
                             choices=all_levels,
                             selected=curr_levels)
      }) # observeEvent

      # observer to update pca coldata if pca col levels change
      observeEvent(input$pca_col_levels, {
        # get full & current metadata
        df_all <- coldata_pca$all
        df_current <- coldata_pca$current

        # if current levels are present in current metadata
        if(all(input$pca_col_levels %in% df_current[,input$pca_cols])){
            # subset current metadata based on current levels
            idx <- df_current[,input$pca_cols] %in% input$pca_col_levels
            coldata_pca$current <- df_current[idx,]
        } else if(all(input$pca_col_levels %in% df_all[,input$pca_cols])){
            # otherwise, update current metadata from full metadata
            idx <- df_all[,input$pca_cols] %in% input$pca_col_levels
            coldata_pca$current <- df_all[idx,]
        } else if(input$pca_col_levels == ''){
            coldata_pca$current <- NULL
        }
      }, ignoreNULL=FALSE) # observeEvent

      # observer to reset PCA cols
      observeEvent(input$reset_pca_cols, {
        coldata_pca$current <- coldata_pca$all
      })

      # reactive to generate pca plot
      pca_plot <- eventReactive(c(data_loaded$flag,
                                  input$plot_do), {
        validate(
          need(!is.null(plot_data$rld), 'Waiting for selection')
        )

        rld <- plot_data$rld

        validate(
          need(input$pcx != input$pcy, 'Please choose different PCs to plot on x and y axes')
        )
        pcx <- input$pcx
        pcy <- input$pcy

        if(input$pca_loadings == 'yes') loadings <- TRUE
        else loadings <- FALSE

        if(is.na(input$pca_loadings_ngenes) | input$pca_loadings_ngenes < 0){
          loadings <- FALSE

          showNotification(
            'Warning: Number of top gene loadings must be positive. Turning off gene loadings',
            type='warning'
          )
        }

        ntop <- input$ntop
        if(is.na(ntop) | ntop < 0){
          showNotification(
            'Warning: Number of genes to use for PCA must be positive. Resetting to default (500)',
            type='warning'
          )

          ntop <- 500
        } else if(ntop > nrow(rld)){
          showNotification(
            paste0('Warning: Number of genes to use for PCA must be < number of genes (',
                   nrow(rld), '). Resetting to default (500)'),
            type='warning'
          )

          ntop <- 500
        }

        samples <- rownames(coldata_pca$current)
        if(length(samples) == 0){
            samples <- coldata_pca$current$samplename
        }

        validate(
          need(input$pca_color != '', 'No coloring variables selected')
        )

        validate(
          need(length(samples) >= 2, 'Must have at least two samples to plot PCA')
        )

        if(input$pcz == 'none'){
          p <- plotPCA.san(rld, input$pca_color, pcx, pcy, ntop=ntop,
                           samples=samples,
                           loadings=loadings,
                           loadings_ngenes=input$pca_loadings_ngenes)
        } else {

          validate(
            need(!input$pcz %in% c(input$pcx, input$pcy),
                 'Please choose different PC to plot on z axis')
          )
          pcz <- input$pcz

          p <- plotPCA.san(rld, input$pca_color, pcx, pcy, pcz=pcz,
                           ntop=ntop,
                           samples=samples,
                           loadings=loadings,
                           loadings_ngenes=input$pca_loadings_ngenes)
        }

        p
      }) # eventReactive pca_plot

      output$pcaplot <- renderPlotly({
        pca_plot()
      })

      ######################### Help buttons #######################

      helpButtonServer('pca_samples_help', size='l')
      helpButtonServer('de_pca_help', size='l')
      downloadButtonServer('pcaplot_download', pca_plot, 'pcaplot')
    }
  ) # moduleServer
} # pcaPlotServer
