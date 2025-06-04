########################## individual modules for plots #########################
#
# General structure:
#
# UI -> sidebar (plot options), main panel (plot)
# server -> plot code
#
########################## blank template ######################################
#
# plotUI <- function(id, panel){
#   ns <- NS(id)
#
#   config <- get_config()
#
#   # NOTE: enter plottype here
#   plottype <- **
#
#   defaults <- config$ui$functional_enrichment$plots[[ plottype ]]
#
#   if(panel == 'sidebar'){
#
#   } else if(panel == 'main'){
#
#   }
# }
#
# plotServer <- function(id, obj){
#
#   moduleServer(
#     id,
#
#     function(input, output, session){
#
#       ns <- NS(id)
#
#       config <- get_config()
#
#       # NOTE: enter plottype here
#       plottype <- **
#
#       defaults <- config$ui$functional_enrichment$plots[[plottype]]
#
#       plot_args <- reactive({
#         # NOTE: list containing plot args
#         # - each element should correspond to inputs listed above
#         list(
#
#         )
#       })
#
#       enrichplot <- eventReactive(c(obj(), plot_args()), {
#         l_gs <- obj()$l_gs
#
#         # NOTE: enter ploting code here
#
#         p
#       })
#
#       output$plot_out <- renderPlot({
#         enrichplot()
#       })
#
#       helpButtonServer(paste0(plottype, '_help'), size='l')
#       downloadButtonServer('enrichplot_download', enrichplot,
#                            plottype)
#     } # function
#   ) # moduleServer
# } # function

######################### Summary overview #######################################

#' Summary overview plot ui
#'
#' @param id ID string used to match the ID used to call the module server function
#' @param panel string, can be 'sidebar' or 'main'
#' @param type string, if 'comp' then show the comparison view
#'
#' @export
sumovPlotUI <- function(id, panel, type=''){
  ns <- NS(id)

  # get default config
  config <- get_config()

  if(panel == 'sidebar'){
    tagList(

      fluidRow(
        column(6, h5('# of terms')),
        column(6,
          numericInput(ns('numcat'), label=NULL,
            value=config$ui$functional_enrichment$plots$summary_overview$numcat
          ) # numericInput
        ) # column
      ), # fluidRow

      fluidRow(
        column(6, h5('p-value column')),
        column(6,
          selectInput(ns('pval'), label=NULL,
            choices=unlist(config$ui$functional_enrichment$plots$summary_overview$pval)
          ) # selectInput
        ) # column
      ), # fluidRow

      fluidRow(
        column(6, h5('color by')),
        column(6,
          selectInput(ns('color'), label=NULL,
            choices=config$ui$functional_enrichment$plots$summary_overview$color
          ) # selectInput
        ) # column
      ), # fluidRow

      fluidRow(
        column(6, h5('max name length')),
        column(6,
          numericInput(ns('catlen'), label=NULL,
            value=config$ui$functional_enrichment$plots$summary_overview$catlen
          ) # numericInput
        ) # column
      ), # fluidRow

      fluidRow(align='center',
        actionButton(ns('plot_do'),
                     label='Refresh plot',
                     icon=icon('arrows-rotate'),
                     class='btn-primary',
                     style='margin-bottom: 10px;')
      ) # fluidRow
    ) # tagList
  } else if(panel == 'main_btns'){
    if(type == 'comp'){
      hlp <- 'comp_summary_overview_help'
    } else {
      hlp <- 'summary_overview_help'
    }

    tagList(
      fluidRow(
        column(6, align='right',
          helpButtonUI(ns(hlp))
        ), # column
        column(6,
          downloadButtonUI(ns('enrichplot_download'))
        ) # column
      ) # fluidRow
    ) # tagList
  } else if(panel == 'main'){
    tagList(
      withSpinner(
        plotOutput(ns('plot_out'), height='600px', width='100%')
      ) # withSpinner
    ) # tagList
  }
} # function


#' Summary overview plot server function
#'
#' @param id ID string used to match the ID used to call the module UI function
#' @param obj reactiveValues object containing GeneTonic object
#' @param config reactive list with config settings
#' @param type string, if 'comp' then show the comparison view
#'
#' @export
sumovPlotServer <- function(id, obj, config, type=''){

  moduleServer(
    id,

    function(input, output, session){

      ns <- NS(id)

      plottype <- 'summary_overview'

      # update from reactive config
      observeEvent(config(), {
        updateNumericInput(session, 'numcat',
                           value=config()$ui$functional_enrichment$plots[[ plottype ]]$numcat)
        updateNumericInput(session, 'catlen',
                           value=config()$ui$functional_enrichment$plots[[ plottype ]]$catlen)
        updateSelectInput(session, 'pval',
                          choices=config()$ui$functional_enrichment$plots[[ plottype ]]$pval)
        updateSelectInput(session, 'color',
                          choices=config()$ui$functional_enrichment$plots[[ plottype ]]$color)
      })

      enrichplot <- eventReactive(c(obj(),
                                    input$plot_do), {
        # get defaults from reactive config
        defaults <- config()$ui$functional_enrichment$plots[[ plottype ]]

        # this is needed to handle reactive ui
        if(is.null(input$numcat)){
          numcat <- defaults$numcat
          catlen <- defaults$catlen
          pcol <- unlist(defaults$pval)[1]
          color_by <- defaults$color[1]
        } else {
          numcat <- input$numcat
          catlen <- input$catlen
          pcol <- input$pval
          color_by <- input$color
        }

        text_size <- config()$server$functional_enrichment$plot[[ plottype ]]$text_size
        validate(
            need(!is.null(numcat) & numcat > 0,
                 'Number of terms must be > 0')
        )

        validate(
            need(!is.null(catlen) & catlen > 0,
                 'Max term name length must be > 0')
        )

        # if a single object is passed
        if('l_gs' %in% names(obj())){
          l_gs <- obj()$l_gs

          validate(
            need(nrow(l_gs) > 0, '')
          )

          l_gs$gs_description <- substr(l_gs$gs_description, 1, catlen)

          n_gs <- min(nrow(l_gs), numcat)

          validate(
              need(sum(duplicated(l_gs$gs_description[1:n_gs])) == 0,
                   'Duplicate category names after truncation. Please increase max name length')
          )

          # find p-value column by matching
          pcol_idx <- grep(pcol, colnames(l_gs))

          # make sure pvalue column exists
          if(length(pcol_idx) == 0){
            showNotification(
              paste0('P-value column: ', pcol, ' does not exist in enrichment table!')
            )

            validate(
              need(pcol %in% colnames(l_gs), '')
            )
          }
          pcol <- colnames(l_gs)[pcol_idx]

          # TODO: expose more params
          p <- gs_summary_overview(l_gs, n_gs = n_gs,
                p_value_column = pcol,
                color_by = color_by) +

               # change y-axis label to show selected pvalue column
               # without a 'gs_' prefix
               ylab(paste('-log10', sub('gs_', '', pcol))) +
               theme(text=element_text(size=text_size))

        } else {
          # when multiple objects are passed
          l_gs <- obj()$obj1$l_gs
          label1 <- obj()$obj1$label

          l_gs$gs_description <- substr(l_gs$gs_description, 1, catlen)

          n_gs1 <- min(nrow(l_gs), numcat)

          validate(
              need(sum(duplicated(l_gs$gs_description[1:n_gs1])) == 0,
                   'Duplicate category names for comparison 1 after truncation. Please increase max name length')
          )

          l_gs2 <- obj()$obj2$l_gs
          label2 <- obj()$obj2$label

          l_gs2$gs_description <- substr(l_gs2$gs_description, 1, catlen)

          n_gs2 <- min(nrow(l_gs2), numcat)

          validate(
              need(sum(duplicated(l_gs2$gs_description[1:n_gs2])) == 0,
                   'Duplicate category names for comparison 2 after truncation. Please increase max name length')
          )

          n_gs <- min(nrow(l_gs), nrow(l_gs2), numcat)

          # find p-value column by matching
          pcol_idx <- grep(pcol, colnames(l_gs))
          # make sure pvalue column exists
          if(length(pcol_idx) == 0){
            showNotification(
              paste0('P-value column: ', pcol, ' does not exist in enrichment table!')
            )

            validate(
              need(pcol %in% colnames(l_gs), '')
            )
          }
          pcol <- colnames(l_gs)[pcol_idx]

          # make sure pvalue column also exists for comparison 2
          if(!pcol %in% colnames(l_gs2)){
            showNotification(
              paste0('P-value column: ', pcol, ' does not exist for comparison 2!')
            )

            validate(
              need(pcol %in% colnames(l_gs2), '')
            )
          }

          p <- gs_summary_overview_pair(l_gs, l_gs2,
                           n_gs = n_gs,
                           p_value_column = pcol,
                           color_by = color_by) +

                # change y-axis label to show selected pvalue column
                # without a 'gs_' prefix
                ylab(paste('-log10', sub('gs_', '', pcol))) +
                theme(text=element_text(size=text_size))

        }
        p
      })

      output$plot_out <- renderPlot({
        enrichplot()
      })

      if(type == 'comp'){
        helpButtonServer(paste0('comp_', plottype, '_help'), size='l')
      } else {
        helpButtonServer(paste0(plottype, '_help'), size='l')
      }
      downloadButtonServer('enrichplot_download', enrichplot,
                           'summary_overview')

    } # function
  ) # moduleServer
} # function

######################### Enrichment Map #######################################

#' Enrichment map plot UI
#'
#' @param id ID string used to match the ID used to call the module server function
#' @param panel string, can be 'sidebar' or 'main'
#'
#' @export
enrichmapUI <- function(id, panel){
  ns <- NS(id)

  # get default config
  config <- get_config()

  # NOTE: enter plottype here
  plottype <- 'enrichment_map'

  if(panel == 'sidebar'){
    tagList(

      fluidRow(
        column(6, h5('# of terms')),
        column(6,
          numericInput(ns('numcat'), label=NULL,
            value=config$ui$functional_enrichment$plots[[ plottype ]]$numcat
          ) # numericInput
        ) # column
      ), # fluidRow

      fluidRow(align='center',
        actionButton(ns('plot_do'),
                     label='Refresh plot',
                     icon=icon('arrows-rotate'),
                     class='btn-primary',
                     style='margin-bottom: 10px;')
      ) # fluidRow
    ) # tagList
  } else if(panel == 'main_btns'){
    tagList(
      fluidRow(
        column(12, align='right',
          helpButtonUI(ns(paste0(plottype, '_help')))
        ) # column
      ) # fluidRow
    ) # tagList
  } else if(panel == 'main'){
    tagList(
      withSpinner(
        visNetworkOutput(ns('plot_out'), height='600px', width='100%')
      ) # withSpinner
    ) # tagList
  }
}

#' Enrichment map plot server function
#'
#' @param id ID string used to match the ID used to call the module UI function
#' @param obj reactiveValues object containing GeneTonic object
#' @param res_obj reactive, dataframe containing enrichment results
#' @param config reactive list with config settings
#'
#' @export
enrichmapServer <- function(id, obj, res_obj, config){

  moduleServer(
    id,

    function(input, output, session){

      ns <- NS(id)

      # NOTE: enter plottype here
      plottype <- 'enrichment_map'

      # update from reactive config
      observeEvent(config(), {
        updateNumericInput(session, 'numcat',
                           value=config()$ui$functional_enrichment$plots[[ plottype ]]$numcat)
      })

      enrichplot <- eventReactive(c(obj(),
                                    input$plot_do), {
        # get defaults from reactive config
        defaults <- config()$ui$functional_enrichment$plots[[ plottype ]]

        l_gs <- obj()$l_gs
        anno_df <- obj()$anno_df
        res <- res_obj()

        validate(
          need(nrow(l_gs) > 0, '')
        )

        # this is needed to handle reactive ui
        if(is.null(input$numcat)){
          numcat <- defaults$numcat
        } else {
          numcat <- input$numcat
        }

        validate(
            need(!is.null(numcat) & numcat > 0,
                 'Number of terms must be > 0')
        )
        n_gs <- min(nrow(l_gs), numcat)

        # TODO: expose more params
        p <- enrichment_map(l_gs, res, anno_df, n_gs=n_gs)

        p
      })

      output$plot_out <- renderVisNetwork({
        enrichplot() %>%
          visIgraph() %>%
          visOptions(highlightNearest = list(enabled = TRUE,
                                             degree = 1,
                                             hover = TRUE),
                     nodesIdSelection = TRUE) %>%
          visNodes(font=list(size=25))
      })

      helpButtonServer(paste0(plottype, '_help'), size='l')
      downloadButtonServer('enrichplot_download', enrichplot,
                           plottype)
    } # function
  ) # moduleServer
} # function

######################### Cnetplot #############################################

#' Cnetplot UI
#'
#' @param id ID string used to match the ID used to call the module server function
#' @param panel string, can be 'sidebar' or 'main'
#'
#' @export
cnetPlotUI <- function(id, panel){
  ns <- NS(id)

  # get default config
  config <- get_config()

  plottype <- 'cnetplot'
  defaults <- config$ui$functional_enrichment$plots[[ plottype ]]

  if(panel == 'sidebar'){

    tagList(
      fluidRow(
        column(6, h5('# of terms')),
        column(6,
          numericInput(ns('numcat'), label=NULL,
            value=defaults$numcat
          ) # numericInput
        ) # column
      ), # fluidRow

      fluidRow(
        column(6, h5('node label')),
        column(6,
          selectInput(ns('node_label'), label=NULL,
            choices=c('category','gene','all','none'),
            selected=defaults$node_label
          ) # selectInput
        ) # column
      ), # fluidRow

      fluidRow(
        column(6, h5('max name length')),
        column(6,
          numericInput(ns('catlen'), label=NULL,
            value=defaults$catlen
          ) # numericInput
        ) # column
      ), # fluidRow

      fluidRow(
        column(7, h5('color edge by terms?')),
        column(5, align='left',
          checkboxInput(ns('colorEdge'), label=NULL,
            value=defaults$color_edge
          ) # checkboxInput
        ) # column
      ), # fluidRow

      fluidRow(
        column(7, h5('circular layout')),
        column(5, align='left',
          checkboxInput(ns('circular'), label=NULL,
            value=defaults$circular
          ) # checkboxInput
        ) # column
      ), # fluidRow

      fluidRow(align='center',
        actionButton(ns('plot_do'),
                     label='Refresh plot',
                     icon=icon('arrows-rotate'),
                     class='btn-primary',
                     style='margin-bottom: 10px;')
      ) # fluidRow

    ) # tagList
  } else if(panel == 'main_btns'){
    tagList(
      fluidRow(
        column(6, align='right',
          helpButtonUI(ns(paste0(plottype, '_help')))
        ), # column
        column(6,
          downloadButtonUI(ns('enrichplot_download'))
        ) # column
      ) # fluidRow
    ) # tagList
  } else if(panel == 'main'){
    tagList(
      withSpinner(
        plotOutput(ns('plot_out'), height='600px', width='100%')
      ) # withSpinner
    ) # tagList
  }
}

#' Cnetplot server function
#'
#' @param id ID string used to match the ID used to call the module UI function
#' @param obj reactive, dataframe with enrichment results
#' @param config reactive list with config settings
#'
#' @export
cnetPlotServer <- function(id, obj, config){

  moduleServer(
    id,

    function(input, output, session){

      ns <- NS(id)

      # NOTE: enter plottype here
      plottype <- 'cnetplot'

      # update from reactive config
      observeEvent(config(), {
        updateNumericInput(session, 'numcat',
                           value=config()$ui$functional_enrichment$plots[[ plottype ]]$numcat)
        updateNumericInput(session, 'catlen',
                           value=config()$ui$functional_enrichment$plots[[ plottype ]]$catlen)
        updateSelectInput(session, 'node_label',
                          choices=config()$ui$functional_enrichment$plots[[ plottype ]]$node_label)
        updateCheckboxInput(session, 'colorEdge',
                            value=config()$ui$functional_enrichment$plots[[ plottype ]]$colorEdge)
        updateCheckboxInput(session, 'circular',
                            value=config()$ui$functional_enrichment$plots[[ plottype ]]$circular)
      })

      enrichplot <- eventReactive(c(obj(),
                                    input$plot_do), {
        l <- obj()

        # get defaults from reactive config
        defaults <- config()$ui$functional_enrichment$plots[[ plottype ]]

        # this is needed to handle reactive ui (run first time)
        if(is.null(input$numcat)){
          numcat <- defaults$numcat
          catlen <- defaults$catlen
          node_label <- defaults$node_label
          colorEdge <- defaults$color_edge
          circular <- defaults$circular
        } else {
          numcat <- input$numcat
          catlen <- input$catlen
          node_label <- input$node_label
          colorEdge <- input$colorEdge
          circular <- input$circular
        }
        validate(
            need(!is.null(catlen) & catlen > 0,
                 'Max term name length must be > 0')
        )
        # convert df to enrichResult
        if(is.data.frame(l)){
          if('core_enrichment' %in% colnames(l)) l <- makeEnrichResult(l, type='gseaResult')
          else l <- makeEnrichResult(l, type='enrichResult')
        }

        l@result$Description <- substr(l@result$Description, 1, catlen)

        text_size <- config()$server$functional_enrichment$plot[[ plottype ]]$text_size
        p <- cnetplot(l,
                      showCategory = numcat,
                      color.params=list(edge=colorEdge),
                      circular = circular,
                      node_label = node_label) + theme(text=element_text(size=text_size))
        p
      })

      output$plot_out <- renderPlot({
        enrichplot()
      })

      helpButtonServer(paste0(plottype, '_help'), size='l')
      downloadButtonServer('enrichplot_download', enrichplot,
                           plottype)
    } # function
  ) # moduleServer
} # function

######################### Radar ################################################

#' Radar plot UI
#'
#' @param id ID string used to match the ID used to call the module server function
#' @param panel string, can be 'sidebar' or 'main'
#' @param type string, if 'comp' then show the comparison view
#'
#' @export
radarUI <- function(id, panel, type=''){
  ns <- NS(id)

  # get default config
  config <- get_config()

  # NOTE: enter plottype here
  plottype <- 'radar'

  defaults <- config$ui$functional_enrichment$plots[[ plottype ]]

  if(panel == 'sidebar'){
    tagList(

      fluidRow(
        column(6, h5('# of terms')),
        column(6,
          numericInput(ns('numcat'), label=NULL,
            value=defaults$numcat
          ) # numericInput
        ) # column
      ), # fluidRow

      fluidRow(
        column(6, h5('p-value column')),
        column(6,
          selectInput(ns('pval'), label=NULL,
            choices=unlist(defaults$pval)
          ) # selectInput
        ) # column
      ), # fluidRow

      fluidRow(
        column(6, h5('max name length')),
        column(6,
          numericInput(ns('catlen'), label=NULL,
            value=defaults$catlen
          ) # numericInput
        ) # column
      ), # fluidRow

      fluidRow(align='center',
        actionButton(ns('plot_do'),
                     label='Refresh plot',
                     icon=icon('arrows-rotate'),
                     class='btn-primary',
                     style='margin-bottom: 10px;')
      ) # fluidRow
    ) # tagList
  } else if(panel == 'main_btns'){
    if(type == 'comp'){
      hlp <- paste0('comp_', plottype, '_help')
    } else {
      hlp <- paste0(plottype, '_help')
    }
    tagList(
      fluidRow(
        column(12, align='right',
          helpButtonUI(ns(hlp))
        ) # column
      ) # fluidRow
    ) # tagList
  } else if(panel == 'main'){
    tagList(
      withSpinner(
        plotlyOutput(ns('plot_out'), height='600px', width='100%')
      ) # withSpinner
    ) # tagList
  }
}

#' Radar plot server function
#'
#' @param id ID string used to match the ID used to call the module UI function
#' @param obj reactive containing GeneTonic object
#' @param type string, if 'comp' then show the comparison view
#' @param config reactive list with config settings
#'
#' @export
radarServer <- function(id, obj, config, type=''){

  moduleServer(
    id,

    function(input, output, session){

      ns <- NS(id)

      # NOTE: enter plottype here
      plottype <- 'radar'

      # update from reactive config
      observeEvent(config(), {
        updateNumericInput(session, 'numcat',
                           value=config()$ui$functional_enrichment$plots[[ plottype ]]$numcat)
        updateNumericInput(session, 'catlen',
                           value=config()$ui$functional_enrichment$plots[[ plottype ]]$catlen)
        updateSelectInput(session, 'pval',
                          choices=config()$ui$functional_enrichment$plots[[ plottype ]]$pval)
      })

      enrichplot <- eventReactive(c(obj(),
                                    input$plot_do), {
        # get defaults from reactive config
        defaults <- config()$ui$functional_enrichment$plots[[plottype]]

        if(is.null(input$numcat)){
          numcat <- defaults$numcat
          pval <- unlist(defaults$pval)[1]
          catlen <- defaults$catlen
        } else {
          numcat <- input$numcat
          pval <- input$pval
          catlen <- input$catlen
        }

        validate(
            need(!is.null(numcat) & numcat > 0,
                 'Number of terms must be > 0')
        )

        validate(
            need(!is.null(catlen) & catlen > 0,
                 'Max term name length must be > 0')
        )

        # when single object is passed
        if('l_gs' %in% names(obj())){
          l_gs <- obj()$l_gs
          label <- obj()$label

          validate(
            need(nrow(l_gs) > 0, '')
          )
          n_gs <- min(nrow(l_gs), numcat)

          l_gs$gs_description <- substr(l_gs$gs_description, 1, catlen)

          validate(
              need(sum(duplicated(l_gs$gs_description[1:n_gs])) == 0,
                   'Duplicate category names after truncation. Please increase max name length')
          )

          # find p-value column by matching
          pval_idx <- grep(pval, colnames(l_gs))
          # make sure pvalue column exists
          if(length(pval_idx) == 0){
            showNotification(
              paste0('P-value column: ', pval, ' does not exist in enrichment table!')
            )

            validate(
              need(pval %in% colnames(l_gs), '')
            )
          }
          pval <- colnames(l_gs)[pval_idx]

          # remove NAs from pvalue column
          na.idx <- is.na(l_gs[, pval])
          l_gs[na.idx, pval] <- 1

          # TODO: expose more params
          p <- gs_radar(l_gs, n_gs = n_gs,
                        label1 = label,
                         p_value_column = pval)
        } else {
          l_gs <- obj()$obj1$l_gs
          label1 <- obj()$obj1$label

          n_gs1 <- min(nrow(l_gs), numcat)

          l_gs$gs_description <- substr(l_gs$gs_description, 1, catlen)

          validate(
              need(sum(duplicated(l_gs$gs_description[1:n_gs1])) == 0,
                   'Duplicate category names for comparison 1 after truncation. Please increase max name length')
          )

          l_gs2 <- obj()$obj2$l_gs
          label2 <- obj()$obj2$label

          n_gs2 <- min(nrow(l_gs2), numcat)

          l_gs2$gs_description <- substr(l_gs2$gs_description, 1, catlen)

          validate(
              need(sum(duplicated(l_gs2$gs_description[1:n_gs2])) == 0,
                   'Duplicate category names for comparison 2 after truncation. Please increase max name length')
          )

          # when multiple objects are passed, run radar for pair
          n_gs <- min(nrow(l_gs), nrow(l_gs2), numcat)

          # find p-value column by matching
          pval_idx <- grep(pval, colnames(l_gs))
          # make sure pvalue column exists
          if(length(pval_idx) == 0){
            showNotification(
              paste0('P-value column: ', pval, ' does not exist in enrichment table!')
            )

            validate(
              need(pval %in% colnames(l_gs), '')
            )
          }
          pval <- colnames(l_gs)[pval_idx]

          # make sure pvalue column also exists for comparison 2
          if(!pval %in% colnames(l_gs2)){
            showNotification(
              paste0('P-value column: ', pval, ' does not exist for comparison 2!')
            )

            validate(
              need(pval %in% colnames(l_gs2), '')
            )
          }

          # remove NAs from pvalue column
          na.idx <- is.na(l_gs[,pval])
          l_gs[na.idx, pval] <- 1

          na.idx2 <- is.na(l_gs2[,pval])
          l_gs2[na.idx2, pval] <- 1

          # adjust labels if comparing within same contrast
          if(label1 == label2){
              label1 <- paste(label1, obj()$obj1$geneset)
              label2 <- paste(label2, obj()$obj2$geneset)
          }

          # TODO: expose more params
          p <- gs_radar(l_gs, l_gs2, n_gs = n_gs,
                        label1 = label1,
                        label2 = label2,
                        p_value_column = pval)
        }
        p
      })

      output$plot_out <- renderPlotly({
        enrichplot()
      })

      if(type == 'comp'){
        helpButtonServer(paste0('comp_', plottype, '_help'), size='l')
      } else {
        helpButtonServer(paste0(plottype, '_help'), size='l')
      }
      downloadButtonServer('enrichplot_download', enrichplot,
                           plottype)
    } # function
  ) # moduleServer
} # function

######################## Alluvial ###############################################

#' Alluvial plot UI
#'
#' @param id ID string used to match the ID used to call the module server function
#' @param panel string, can be 'sidebar' or 'main'
#'
#' @export
alluvialUI <- function(id, panel){
  ns <- NS(id)

  # get default config
  config <- get_config()

  # NOTE: enter plottype here
  plottype <- 'alluvial'

  defaults <- config$ui$functional_enrichment$plots[[ plottype ]]

  if(panel == 'sidebar'){
    tagList(

      fluidRow(
        column(6, h5('# of terms')),
        column(6,
          numericInput(ns('numcat'), label=NULL,
            value=defaults$numcat
          ) # numericInput
        ) # column
      ), # fluidRow

      fluidRow(align='center',
        actionButton(ns('plot_do'),
                     label='Refresh plot',
                     icon=icon('arrows-rotate'),
                     class='btn-primary',
                     style='margin-bottom: 10px;')
      ) # fluidRow
    ) # tagList
  } else if(panel == 'main_btns'){
    tagList(
      fluidRow(
        column(12, align='right',
          helpButtonUI(ns(paste0(plottype, '_help')))
        ) # column
      ) # fluidRow
    ) # tagList
  } else if(panel == 'main'){
    tagList(
      withSpinner(
        plotlyOutput(ns('plot_out'), height='600px', width='100%')
      ) # withSpinner
    ) # tagList
  }
}

#' Alluvial plot server function
#'
#' @param id ID string used to match the ID used to call the module UI function
#' @param obj reactive containing GeneTonic object
#' @param res_obj reactive, dataframe containing enrichment results
#' @param config reactive list with config settings
#'
#' @export
alluvialServer <- function(id, obj, res_obj, config){

  moduleServer(
    id,

    function(input, output, session){

      ns <- NS(id)

      # NOTE: enter plottype here
      plottype <- 'alluvial'

      observeEvent(config(), {
        updateNumericInput(session, 'numcat',
                           value=config()$ui$functional_enrichment$plots[[ plottype ]]$numcat)
      })

      enrichplot <- eventReactive(c(obj(),
                                    input$plot_do), {
        l_gs <- obj()$l_gs
        anno_df <- obj()$anno_df
        res <- res_obj()

        defaults <- config()$ui$functional_enrichment$plots[[plottype]]

        if(is.null(input$numcat)){
          numcat <- defaults$numcat
        } else {
          numcat <- input$numcat
        }

        validate(
            need(!is.null(numcat) & numcat > 0,
                 'Number of terms must be > 0')
        )

        n_gs <- min(nrow(l_gs), numcat)

        p <- gs_alluvial(l_gs, res, anno_df, n_gs=n_gs)

        p
      })

      output$plot_out <- renderPlotly({
        enrichplot()
      })

      helpButtonServer(paste0(plottype, '_help'), size='l')
      downloadButtonServer('enrichplot_download', enrichplot,
                           plottype)
    } # function
  ) # moduleServer
} # function

#################### Dendrogram ################################################

#' Dendrogram plot UI
#'
#' @param id ID string used to match the ID used to call the module server function
#' @param panel string, can be 'sidebar' or 'main'
#'
#' @export
dendrogramUI <- function(id, panel){
  ns <- NS(id)

  # get default config
  config <- get_config()

  # NOTE: enter plottype here
  plottype <- 'dendrogram'

  defaults <- config$ui$functional_enrichment$plots[[ plottype ]]

  if(panel == 'sidebar'){
    tagList(

      fluidRow(
        column(6, h5('# of terms')),
        column(6,
          numericInput(ns('numcat'), label=NULL,
            value=defaults$numcat
          ) # numericInput
        ) # column
      ), # fluidRow

      fluidRow(
        column(6, h5('max name length')),
        column(6,
          numericInput(ns('catlen'), label=NULL,
            value=defaults$catlen
          ) # numericInput
        ) # column
      ), # fluidRow

      fluidRow(align='center',
        actionButton(ns('plot_do'),
                     label='Refresh plot',
                     icon=icon('arrows-rotate'),
                     class='btn-primary',
                     style='margin-bottom: 10px;')
      ) # fluidRow
    ) # tagList
  } else if(panel == 'main_btns'){
    tagList(
      fluidRow(
        column(6, align='right',
          helpButtonUI(ns(paste0(plottype, '_help')))
        ), # column
        column(6,
          downloadButtonUI(ns('enrichplot_download'))
        ) # column
      ) # fluidRow
    ) # tagList
  } else if(panel == 'main'){
    tagList(
      withSpinner(
        plotOutput(ns('plot_out'), height='600px', width='100%')
      ) # withSpinner
    ) # tagList
  }
}

#' Dendrogram plot server function
#'
#' @param id ID string used to match the ID used to call the module UI function
#' @param obj reactive containing GeneTonic object
#' @param config reactive list with config settings
#'
#' @export
dendrogramServer <- function(id, obj, config){

  moduleServer(
    id,

    function(input, output, session){

      ns <- NS(id)

      # NOTE: enter plottype here
      plottype <- 'dendrogram'

      # update from reactive config
      observeEvent(config(), {
        updateNumericInput(session, 'numcat',
                           value=config()$ui$functional_enrichment$plots[[ plottype ]]$numcat)
        updateNumericInput(session, 'catlen',
                           value=config()$ui$functional_enrichment$plots[[ plottype ]]$catlen)
      })

      enrichplot <- eventReactive(c(obj(),
                                    input$plot_do), {
        l_gs <- obj()$l_gs

        validate(
          need(nrow(l_gs) > 0, '')
        )

        # get defaults from reactive config
        defaults <- config()$ui$functional_enrichment$plots[[plottype]]

        # this is needed to handle reactive ui
        if(is.null(input$numcat)){
          numcat <- defaults$numcat
          catlen <- defaults$catlen
        } else {
          numcat <- input$numcat
          catlen <- input$catlen
        }
        validate(
            need(!is.null(numcat) & numcat > 0,
                 'Number of terms must be > 0')
        )
        n_gs <- min(nrow(l_gs), numcat)

        validate(
            need(!is.null(catlen) & catlen > 0,
                 'Max term name length must be > 0')
        )
        l_gs$gs_description <- substr(l_gs$gs_description, 1, catlen)

        # TODO: expose more params
        p <- gs_dendro(l_gs,
                       n_gs = n_gs,
                       gs_dist_type = "kappa",
                       clust_method = "ward.D2",
                       color_leaves_by = "z_score",
                       size_leaves_by = "gs_pvalue",
                       color_branches_by = "clusters",
                       create_plot = FALSE)

        p
      })

      output$plot_out <- renderPlot({
        p <- enrichplot()

        validate(
          need(!is.null(p), '')
        )

        old.par <- par(mar = c(0, 0, 1, 25), cex=1.25)
        plot(p, horiz=TRUE)
        par(old.par)
        # p %>% as.ggdend() %>%
        #     ggplot(horiz=TRUE, offset_labels=-0.1)
      })

      helpButtonServer(paste0(plottype, '_help'), size='l')
      downloadButtonServer('enrichplot_download', enrichplot,
                           plottype)
    } # function
  ) # moduleServer
} # function


######################### Horizon #########################

#' Horizon plot UI
#'
#' @param id ID string used to match the ID used to call the module server function
#' @param panel string, can be 'sidebar' or 'main'
#'
#' @export
horizonUI <- function(id, panel){
  ns <- NS(id)

  # get default config
  config <- get_config()

  # NOTE: enter plottype here
  plottype <- 'horizon'

  defaults <- config$ui$functional_enrichment$plots[[ plottype ]]

  if(panel == 'sidebar'){
    tagList(

      fluidRow(
        column(6, h5('# of terms')),
        column(6,
          numericInput(ns('numcat'), label=NULL,
            value=defaults$numcat
          ) # numericInput
        ) # column
      ), # fluidRow

      fluidRow(
        column(6, h5('sort by')),
        column(6,
          selectInput(ns('sort_by'), label=NULL,
            choices=defaults$sort_by
          ) # selectInput
        ) # column
      ), # fluidRow

      fluidRow(
        column(6, h5('max name length')),
        column(6,
          numericInput(ns('catlen'), label=NULL,
            value=defaults$catlen
          ) # numericInput
        ) # column
      ), # fluidRow

      fluidRow(align='center',
        actionButton(ns('plot_do'),
                     label='Refresh plot',
                     icon=icon('arrows-rotate'),
                     class='btn-primary',
                     style='margin-bottom: 10px;')
      ) # fluidRow
    ) # tagList
  } else if(panel == 'main_btns'){
    tagList(
      fluidRow(
        column(6, align='right',
          helpButtonUI(ns(paste0('comp_', plottype, '_help')))
        ), # column
        column(6,
          downloadButtonUI(ns('compenrich_download'))
        ) # column
      ) # fluidRow
    ) # tagList
  } else if(panel == 'main'){
    tagList(
      withSpinner(
        plotOutput(ns('plot_out'), height='600px', width='100%')
      ) # withSpinner
    )
  }
}

#' Horizon plot server function
#'
#' @param id ID string used to match the ID used to call the module UI function
#' @param obj reactive containing GeneTonic object
#' @param config reactive list with config settings
#'
#' @export
horizonServer <- function(id, obj, config){

  moduleServer(
    id,

    function(input, output, session){

      ns <- NS(id)

      # NOTE: enter plottype here
      plottype <- 'horizon'

      plot_args <- reactive({
        # NOTE: list containing plot args
        # - each element should correspond to inputs listed above
        list(
          numcat=input$numcat,
          sort_by=input$sort_by,
          catlen=input$catlen
        )
      })

      # update from reactive config
      observeEvent(config(), {
        updateNumericInput(session, 'numcat',
                           value=config()$ui$functional_enrichment$plots[[ plottype ]]$numcat)
        updateNumericInput(session, 'catlen',
                           value=config()$ui$functional_enrichment$plots[[ plottype ]]$catlen)
        updateSelectInput(session, 'sort_by',
                          choices=config()$ui$functional_enrichment$plots[[ plottype ]]$sort_by)
      })

      enrichplot <- eventReactive(c(obj(),
                                    input$plot_do), {

        # get defaults from reactive config
        defaults <- config()$ui$functional_enrichment$plots[[ plottype ]]

        # check if inputs exist (run first time/before plot options is uncollapsed)
        if(is.null(input$numcat)){
          numcat <- defaults$numcat
          sort_by <- unlist(defaults$sort_by)[1]
          catlen <- defaults$catlen
        } else {
          numcat <- input$numcat
          sort_by <- input$sort_by
          catlen <- input$catlen
        }
        validate(
          need(!is.null(numcat) & numcat > 0,
               'Number of terms must be > 0')
        )

        l_gs1 <- obj()$obj1$l_gs
        label1 <- obj()$obj1$label

        n_gs1 <- min(nrow(l_gs1), numcat)

        l_gs1$gs_description <- substr(l_gs1$gs_description, 1, catlen)

        validate(
            need(sum(duplicated(l_gs1$gs_description[1:n_gs1])) == 0,
                 'Duplicate category names for comparison 1 after truncation. Please increase max name length')
        )

        l_gs2 <- obj()$obj2$l_gs
        label2 <- obj()$obj2$label

        n_gs2 <- min(nrow(l_gs2), numcat)

        l_gs2$gs_description <- substr(l_gs2$gs_description, 1, catlen)

        validate(
            need(sum(duplicated(l_gs2$gs_description[1:n_gs2])) == 0,
                 'Duplicate category names for comparison 2 after truncation. Please increase max name length')
        )

        text_size <- config()$server$functional_enrichment$plot[[ plottype ]]$text_size
        comp_list <- list(l_gs2)

        # adjust group names if comparing within same contrast
        if(label2 == label1){
            names(comp_list) <- paste(label2, obj()$obj2$geneset)
            ref_name <- paste(label1, obj()$obj1$geneset)
        } else {
            names(comp_list) <- label2
            ref_name <- label1
        }

        n_gs <- min(nrow(l_gs1), nrow(l_gs2), numcat)

        # TODO: expose more params
        p <- gs_horizon(l_gs1,
           compared_res_enrich_list = comp_list,
           ref_name = ref_name,
           n_gs = n_gs,
           sort_by = sort_by) + theme(text=element_text(size=text_size))

        p
      })

      output$plot_out <- renderPlot({
        enrichplot()
      })

      helpButtonServer(paste0('comp_', plottype, '_help'), size='l')
      downloadButtonServer('enrichplot_download', enrichplot,
                           plottype)
    } # function
  ) # moduleServer
} # function


######################## Distill Enrichment ####################################

#' Distill enrichment plot UI
#'
#' @param id ID string used to match the ID used to call the module server function
#' @param panel string, can be 'sidebar' or 'main'
#'
#' @export
distillPlotUI <- function(id, panel){
  ns <- NS(id)

  config <- get_config()

  # NOTE: enter plottype here
  plottype <- 'emap_distill'

  if(panel == 'sidebar'){
    tagList(

      fluidRow(
        column(5, h5('# of terms')),
        column(7,
          numericInput(ns('numcat'), label=NULL,
                       value=NA
          ) # numericInput
        ) # column
      ), # fluidRow

      fluidRow(align='center',
        actionButton(ns('plot_do'),
                     label='Refresh plot',
                     icon=icon('arrows-rotate'),
                     class='btn-primary',
                     style='margin-bottom: 10px;')
      ) # fluidRow
    ) # tagList
  } else if(panel == 'main_btns'){
    tagList(
      fluidRow(
        column(12, align='right',
          helpButtonUI(ns(paste0(plottype, '_help')))
        ) # column
      ) # fluidRow
    ) # tagList
  } else if(panel == 'main'){
    tagList(
      withSpinner(
        visNetworkOutput(ns('plot_out'), height='600px', width='100%')
      ) # withSpinner
    ) # tagList
  }
}

#' Distill enrichment plot server function
#'
#' @param id ID string used to match the ID used to call the module UI function
#' @param obj reactive containing 'distilled' enrichment results
#' @param args reactive, list with plot arguments, 'numcat' (number of categories to plot)
#'
#' @export
distillPlotServer <- function(id, obj, args){

  moduleServer(
    id,

    function(input, output, session){

      ns <- NS(id)

      config <- get_config()

      # NOTE: enter plottype here
      plottype <- 'emap_distill'

      defaults <- config$ui$functional_enrichment$plots[[ plottype ]]

      #plot_args <- reactive({
      #  # NOTE: list containing plot args
      #  # - each element should correspond to plot options listed above
      #  list(
      #    numcat=curr_args$numcat
      #  )
      #})

      curr_args <- reactiveValues(numcat=NULL)

      observeEvent(args()$numcat, {
        numcat <- args()$numcat

        updt <- TRUE
        if(is.null(numcat)){
          curr_args$numcat <- defaults$numcat
          updt <- TRUE
        } else if(is.null(curr_args$numcat)){
          curr_args$numcat <- numcat
          updt <- TRUE
        } else if(numcat != curr_args$numcat){
          curr_args$numcat <- numcat
          updt <- TRUE
        }

        if(updt){
          updateNumericInput(session, 'numcat',
                             value=curr_args$numcat)
        }
      })

      observeEvent(input$numcat, {
        if(is.null(curr_args$numcat) & !is.na(input$numcat)){
          curr_args$numcat <- input$numcat
        } else if(!is.na(input$numcat) & curr_args$numcat != input$numcat){
          curr_args$numcat <- input$numcat
        } else {
          isolate({
            updateNumericInput(session, 'numcat',
                               value=curr_args$numcat)
          })
        }
      })

      enrichplot <- eventReactive(c(obj(), input$plot_do), {
        p <- obj()$distilled_em

        # defining a color palette for nicer display
        colpal <- rainbow_hcl(length(unique(igraph::V(p)$color)))[igraph::V(p)$color]
        igraph::V(p)$color.background <- alpha(colpal, alpha = 0.8)
        igraph::V(p)$color.highlight <- alpha(colpal, alpha = 1)
        igraph::V(p)$color.hover <- alpha(colpal, alpha = 0.5)
        igraph::V(p)$color.border <- "black"

        p
      })

      output$plot_out <- renderVisNetwork({
        enrichplot() %>%
          visIgraph() %>%
          visOptions(highlightNearest = list(enabled = TRUE,
                                             degree = 1,
                                             hover = TRUE),
                     nodesIdSelection = TRUE, selectedBy = 'membership') %>%
          visNodes(font=list(size=25))
      })

      helpButtonServer(paste0(plottype, '_help'), size='l')


      plot_data <- eventReactive(input$plot_do, {
        curr_args$numcat
      })

      return(
        plot_data
      )
    } # function
  ) # moduleServer
} # function

######################## Fuzzy Enrichment ########################

#' Fuzzy enrichment plot UI
#'
#' @param id ID string used to match the ID used to call the module server function
#' @param panel string, can be 'sidebar' or 'main'
#'
#' @export
fuzzyPlotUI <- function(id, panel){
  ns <- NS(id)

  config <- get_config()

  # NOTE: enter plottype here
  plottype <- 'emap_fuzzy'

  if(panel == 'sidebar'){
    tagList(

      fluidRow(
        column(5, h5('# of terms')),
        column(7,
          numericInput(ns('numcat'), label=NULL,
                       value=NA
          ) # numericInput
        ) # column
      ), # fluidRow

      fluidRow(align='center',
        actionButton(ns('plot_do'),
                     label='Refresh plot',
                     icon=icon('arrows-rotate'),
                     class='btn-primary',
                     style='margin-bottom: 10px;')
      ) # fluidRow
    ) # tagList
  } else if(panel == 'main_btns'){
    tagList(
      fluidRow(
        column(12, align='right',
          helpButtonUI(ns(paste0(plottype, '_help')))
        ) # column
      ) # fluidRow
    ) # tagList
  } else if(panel == 'main'){
    tagList(
      withSpinner(
        visNetworkOutput(ns('plot_out'), height='600px', width='100%')
      ) # withSpinner
    ) # tagList
  }
}

#' Fuzzy enrichment plot server function
#'
#' @param id ID string used to match the ID used to call the module UI function
#' @param obj reactive containing fuzzy enrichment object
#' @param args reactive, list with plot arguments, 'numcat' (number of categories to plot)
#'
#' @export
fuzzyPlotServer <- function(id, obj, args){

  moduleServer(
    id,

    function(input, output, session){

      ns <- NS(id)

      config <- get_config()

      # NOTE: enter plottype here
      plottype <- 'emap_fuzzy'

      defaults <- config$ui$functional_enrichment$plots[[ plottype ]]

      #plot_args <- reactive({
      #  # NOTE: list containing plot args
      #  # - each element should correspond to plot options listed above
      #  list(
      #    numcat=curr_args$numcat
      #  )
      #})

      curr_args <- reactiveValues(numcat=NULL)

      observeEvent(args()$numcat, {
        numcat <- args()$numcat

        updt <- TRUE
        if(is.null(numcat)){
          curr_args$numcat <- defaults$numcat
          updt <- TRUE
        } else if(is.null(curr_args$numcat)){
          curr_args$numcat <- numcat
          updt <- TRUE
        } else if(numcat != curr_args$numcat){
          curr_args$numcat <- numcat
          updt <- TRUE
        }

        if(updt){
          updateNumericInput(session, 'numcat',
                             value=curr_args$numcat)
        }
      })

      observeEvent(input$numcat, {
        if(is.null(curr_args$numcat) & !is.na(input$numcat)){
          curr_args$numcat <- input$numcat
        } else if(!is.na(input$numcat) & curr_args$numcat != input$numcat){
          curr_args$numcat <- input$numcat
        } else {
          isolate({
            updateNumericInput(session, 'numcat',
                               value=curr_args$numcat)
          })
        }
      })

      enrichplot <- eventReactive(c(obj(), input$plot_do), {
        l <- obj()

        # find p.adjust column and stop if not present
        padj_col <- grep('p.adjust', colnames(l$clusters))
        if(length(padj_col) == 0){
          showNotification(
            'Adjusted p-value column not found in GeneTonic object!'
          )

          validate(
            need(length(padj_col) > 0, '')
          )
        }

        # check for other columns
        cols_to_check <- c('gs_description', 'gs_fuzzycluster', 'gs_cluster_status')
        if(any(!cols_to_check %in% colnames(l$clusters))){
          missing <- setdiff(cols_to_check, colnames(l$clusters))
          showNotification(
            paste0('Fuzzy table missing required columns:',
                   paste(missing, collapse=','))
          )

          validate(
            need(length(missing) == 0, '')
          )
        }

        # TODO: expose more params
        p <- enrichment_map(l$clusters[order(l$clusters[[ padj_col ]]),],
                            l$res, l$anno_df,
                            n_gs = nrow(l$clusters),
                            color_by = colnames(l$clusters)[padj_col])

        # add info from data frame to graph
        term.order <- match(igraph::V(p)$name, l$clusters$gs_description)
        igraph::V(p)$color <- l$clusters$gs_fuzzycluster[term.order]
        igraph::V(p)$membership <- l$clusters$gs_fuzzycluster[term.order]
        igraph::V(p)$status <- l$clusters$gs_cluster_status[term.order]


        # defining a color palette for nicer display
        colpal <- qualitative_hcl(length(unique(igraph::V(p)$color)))[igraph::V(p)$color]
        igraph::V(p)$color.background <- alpha(colpal, alpha = 0.8)
        igraph::V(p)$color.highlight <- alpha(colpal, alpha = 1)
        igraph::V(p)$color.hover <- alpha(colpal, alpha = 0.5)
        igraph::V(p)$color.border <- ifelse(igraph::V(p)$status == "Representative",
                                    "black", "darkgrey")

        p
      })

      output$plot_out <- renderVisNetwork({
        enrichplot() %>%
          visIgraph() %>%
          visOptions(highlightNearest = list(enabled = TRUE,
                                             degree = 1,
                                             hover = TRUE),
                     nodesIdSelection = TRUE, selectedBy = 'membership') %>%
          visNodes(font=list(size=25))
      })

      helpButtonServer(paste0(plottype, '_help'), size='l')

      plot_data <- eventReactive(input$plot_do, {
        curr_args$numcat
      })

      return(
        plot_data
      )
    } # function
  ) # moduleServer
} # function

