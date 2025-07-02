#' MA plot module UI
#'
#' @param id ID string used to match the ID used to call the module server function
#' @param panel string, can be 'sidebar' or 'main'
#'
#' @export
maPlotUI <- function(id, panel){
  ns <- NS(id)

  # get default config
  config <- get_config()

  if(panel == 'sidebar'){
    tagList(
      fluidRow(
        column(6, align='left',
          tags$label(class='control-label',
                     'Comparison')
        ),
        column(6, align='right',
          helpButtonUI(ns('ma_controls_help'))
        ) # column
      ), # fluidRow

      selectizeInput(ns('comp_all'),
                    label=NULL,
                    choices=NULL,
                    selected=NULL
      ), # selectizeInput

      fluidRow(
        column(4, h5('Interactive?')),
        column(8,
          selectInput(ns("plot_interactive"), label=NULL,
                      choices=c('yes', 'no')
          ) # selectInput
        ) # column
      ), # fluidRow

      bsCollapse(
        bsCollapsePanel('Plot options',

          tags$label(class='control-label',
                     'y-axis limits'),

          fluidRow(
            column(4, h5('max')),
            column(8,
              numericInput(ns('ma_ymax'), label=NULL,
                value=config$ui$de_analysis$ma_plot$log2fc_limits$max)
            ) # column
          ), # fluidRow

          fluidRow(
            column(4, h5('min')),
            column(8,
              numericInput(ns('ma_ymin'), label=NULL,
                value=config$ui$de_analysis$ma_plot$log2fc_limits$min)
            ) # column
          ), # fluidRow

          fluidRow(
            column(4, align='left', style='margin-bottom: 10px;',
              actionButton(ns('ma_y_auto'), label='Autoscale')
            ) # column
          ) # fluidRow

        ) # bsCollapsePanel
      ) # bsCollapse

    ) # tagList
  } else if(panel == 'main'){
    tagList(
      fluidRow(
        column(6, align='left',
          helpButtonUI(ns('de_ma_help'))
        ), # column
        column(6, align='right',
          downloadButtonUI(ns('maplot_download'))
        ) # column
      ), # fluidRow

      withSpinner(
        uiOutput(ns('maplot_out'))
      ) # withSpinner
    ) # tagList
  }
}

#' MA plot module server function
#'
#' @param id ID string used to match the ID used to call the module UI function
#' @param obj reactiveValues object containing carnation object
#' @param plot_args reactive containing 'fdr.thres' (padj threshold), 'fc.thres' (log2FC threshold)
#' & 'gene.to.plot' (genes selected in scratchpad)
#' @param config reactive list with config settings
#'
#' @export
maPlotServer <- function(id, obj, plot_args, config){

  moduleServer(
    id,

    function(input, output, session){

      ns <- NS(id)

      app_object <- reactive({
          list(res=obj$res)
      })

      observeEvent(app_object()$res, {
        validate(
          need(!is.null(app_object()$res), 'Waiting for data')
        )

        updateSelectInput(session, 'comp_all',
                          choices=names(app_object()$res))
      })

      curr_thres <- reactiveValues(fdr.thres=0.1,
                                   fc.thres=0)

      # update from reactive config
      observeEvent(config(), {
        curr_thres$fdr.thres <- config()$ui$de_analysis$filters$fdr_threshold
        curr_thres$fc.thres <- config()$ui$de_analysis$filters$log2fc_threshold

        updateNumericInput(session, 'ma_ymax',
                           value=config()$ui$de_analysis$ma_plot$log2fc_limits$max)
        updateNumericInput(session, 'ma_ymin',
                           value=config()$ui$de_analysis$ma_plot$log2fc_limits$min)
      })

      observeEvent(c(plot_args()$fdr.thres, plot_args()$fc.thres), {
        fc.thres <- ifelse(plot_args()$fc.thres == '' | is.na(plot_args()$fc.thres),
                           config()$ui$de_analysis$filters$log2fc_threshold,
                           plot_args()$fc.thres)
        fdr.thres <- ifelse(plot_args()$fdr.thres == '' | is.na(plot_args()$fdr.thres),
                            config()$ui$de_analysis$filters$fdr_threshold,
                            plot_args()$fdr.thres)

        curr_thres$fdr.thres <- fdr.thres
        curr_thres$fc.thres <- fc.thres
      })

      # reactive to generate labeled MA plot
      # NOTE: this is used only for the downloaded plot
      maplot <- eventReactive(c(input$comp_all, plot_args()$gene.to.plot,
                                curr_thres$fdr.thres, curr_thres$fc.thres,
                                input$ma_ymax, input$ma_ymin), {
        validate(
          need(!is.null(app_object()$res) & !is.null(input$comp_all) & input$comp_all != '',
               'Waiting for selection')
        )

        validate(
          need(input$comp_all %in% names(app_object()$res), 'Waiting for selection')
        )

        if(is.null(plot_args()$gene.to.plot) | all(plot_args()$gene.to.plot %in% '')){
            lab.genes <- NULL
        } else {
            lab.genes <- plot_args()$gene.to.plot
        }

        validate(
          need(input$ma_ymin != '' & input$ma_ymax != '',
               'y-axis limits missing')
        )

        validate(
          need(input$ma_ymin < input$ma_ymax,
               'y-axis min must be < y-axis max')
        )
        plotMA.label(app_object()$res[[input$comp_all]],
                     fdr.thres=curr_thres$fdr.thres,
                     fc.thres=curr_thres$fc.thres,
                     fc.lim=c(input$ma_ymin, input$ma_ymax),
                     lab.genes=lab.genes)
      }) # eventReactive maplot

      # observer for maplot ylim autoscale btn
      observeEvent(input$ma_y_auto, {
        showNotification('Autoscaling y-axis limits')

        df <- app_object()$res[[input$comp_all]]

        # add a 5% buffer and round to 3 significant digits
        df.max <- round(max(df$log2FoldChange, na.rm=T)*1.05,
                        digits=3)
        df.min <- round(min(df$log2FoldChange, na.rm=T)*1.05,
                        digits=3)

        updateNumericInput(session, 'ma_ymin', value=df.min)
        updateNumericInput(session, 'ma_ymax', value=df.max)
      })


      # this is the interactive plot_ly version
      maplot_ly <- eventReactive(c(app_object()$res, input$comp_all, plot_args()$gene.to.plot,
                                curr_thres$fdr.thres, curr_thres$fc.thres,
                                input$ma_ymin, input$ma_ymax), {
        validate(
          need(!is.null(app_object()$res) & !is.null(input$comp_all) & input$comp_all != '',
               'Waiting for selection')
        )

        validate(
          need(input$comp_all %in% names(app_object()$res), 'Waiting for selection')
        )

        g <- plot_args()$gene.to.plot
        if(is.null(g) || all(g %in% '')){
            lab.genes <- NULL
        } else {
            lab.genes <- g
        }

        validate(
          need(input$ma_ymin != '' & input$ma_ymax != '',
               'y-axis limits missing')
        )

        validate(
          need(input$ma_ymin < input$ma_ymax,
               'y-axis min must be < y-axis max')
        )

        plotMA.label_ly(app_object()$res[[input$comp_all]],
                 fdr.thres=curr_thres$fdr.thres,
                 fc.thres=curr_thres$fc.thres,
                 fc.lim=c(input$ma_ymin, input$ma_ymax),
                 lab.genes=lab.genes)
      }) # eventReactive maplot_ly

      output$maplot_out <- renderUI({
        if(input$plot_interactive == 'yes'){
          p <- maplot_ly() %>% toWebGL()

          output$plot1 <- renderPlotly({ p })

          withSpinner(
            plotlyOutput(ns('plot1'), height='600px')
          )
        } else if(input$plot_interactive == 'no'){
          p <- maplot() + theme(text=element_text(size=18))

          output$plot2 <- renderPlot({ p })

          withSpinner(
            plotOutput(ns('plot2'), height='600px')
          )
        }
      })

      helpButtonServer('ma_controls_help')
      helpButtonServer('de_ma_help', size='l')
      downloadButtonServer('maplot_download', maplot, 'maplot')
    }
  ) # moduleServer
} # maPlotServer
