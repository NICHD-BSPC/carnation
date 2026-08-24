volcanoPlotUI <- function(id, panel) {
  ns <- NS(id)
  config <- get_config()

  ########### Sidebar ################
  if (panel == 'sidebar') {
    tagList(
      fluidRow(
        column(
          6,
          align = 'left',
          tags$label(class = 'control-label', 'Comparison')
        ),
        column(6, align = 'right', helpButtonUI(ns('volcano_controls_help')))
      ), # fluidRow

      selectizeInput(
        ns('comp_all'),
        label = NULL,
        choices = NULL,
        selected = NULL
      ), # selectizeInput
      fluidRow(
        column(4, h5('Interactive?')),
        column(
          8,
          selectInput(
            ns("plot_interactive"),
            label = NULL,
            choices = c('yes', 'no')
          )
        ) # column
      ), # fluid row

      ############## Plot Options Menu ###############

      bsCollapse(
        id = ns('plot_opts'),
        bsCollapsePanel(
          'Plot options',
          fluidRow(
            column(4, h5('Color by')),
            column(
              8,
              selectInput(
                ns('color_by'),
                label = NULL,
                choices = c('baseMean', 'significance'),
                selected = 'baseMean'
              )
            ) # column
          ), # fluidRow

          ## alpha-slider ##########################
          fluidRow(
            column(
              8,
              sliderInput(ns('volcano_alpha'), 'Opacity/Alpha value', 0, 1, 0.6)
            )
          ),

          ## x-axis limits ##########################
          tags$label(class = 'control-label', 'x-axis limits'),
          fluidRow(
            column(4, h5('max')),
            column(
              8,
              numericInput(
                ns('volcano_xmax'),
                label = NULL,
                value = config$ui$de_analysis$volcano_plot$log2fc_limits$max
              )
            ) # column
          ), # fluidRow
          fluidRow(
            column(4, h5('min')),
            column(
              8,
              numericInput(
                ns('volcano_xmin'),
                label = NULL,
                value = config$ui$de_analysis$volcano_plot$log2fc_limits$min
              )
            ) # column
          ), # fluidRow

          ## y-axis limits ########################
          tags$label(class = 'control-label', 'y-axis limits'),
          fluidRow(
            column(4, h5('max')),
            column(
              8,
              numericInput(
                ns('volcano_ymax'),
                label = NULL,
                value = config$ui$de_analysis$volcano_plot$neg_log_padj_limits$max
              )
            ) # column
          ), # fluidRow
          fluidRow(
            column(4, h5('min')),
            column(
              8,
              numericInput(
                ns('volcano_ymin'),
                label = NULL,
                value = config$ui$de_analysis$volcano_plot$neg_log_padj_limits$min
              )
            ) # column
          ), # fluidRow

          ## autoscaler #########################
          fluidRow(
            column(
              4,
              align = 'left',
              style = 'margin-bottom: 10px;',
              actionButton(ns('volcano_auto'), label = 'Autoscale')
            ) # column
          ) # fluidRow
        ) # bsCollapsePanel
      ) # bsCollapse
    )

    ################# MAIN ##########################
  } else if (panel == 'main') {
    tagList(
      fluidRow(
        column(6, align = 'left', helpButtonUI(ns('de_volcano_help'))),
        column(
          6,
          align = 'right',
          downloadButtonUI(ns('volcano_plot_download'))
        ),
      ), # fluidRow
      withSpinner(
        uiOutput(ns('volcano_plot_out'))
      ) # withSpinner
    ) # tagList
  }
}

################# SERVER ######################
volcanoPlotServer <- function(id, obj, plot_args, config) {
  moduleServer(
    id,
    function(input, output, session) {
      ns <- NS(id)

      # loads in the app obj containing the deseq2 processed data
      app_object <- reactive({
        list(res = obj$res)
      })

      helpButtonServer('volcano_controls_help')
      helpButtonServer('de_volcano_help', size = 'l')

      # Watches the app_object to update the drop down menu with the correct
      # options
      observeEvent(app_object()$res, {
        validate(need(!is.null(app_object()$res), 'waiting for data'))
        updateSelectizeInput(
          session,
          'comp_all',
          choices = names(app_object()$res)
        )
      })

      # Instantiates the reactive values used in making the plots
      curr_thres <- reactiveValues(
        fdr.thres = 0.1,
        fc.thres = 0.0,
        colorscale = 'viridis'
      )

      # Loads in the settings from the config
      observeEvent(config(), {
        curr_thres$fdr.thres <- config()$ui$de_analysis$filters$fdr_threshold
        curr_thres$fc.thres <- config()$ui$de_analysis$filters$log2fc_threshold

        # gets colorscale from config and falls back to viridis if it is not found
        cs <- config()$ui$de_analysis$volcano_plot$colorscale
        curr_thres$colorscale <- if (!is.null(cs)) cs else 'viridis'

        updateNumericInput(
          session,
          'volcano_xmax',
          value = config()$ui$de_analysis$volcano_plot$log2fc_limits$max
        )
        updateNumericInput(
          session,
          'volcano_xmin',
          value = config()$ui$de_analysis$volcano_plot$log2fc_limits$min
        )
        updateNumericInput(
          session,
          'volcano_ymax',
          value = config()$ui$de_analysis$volcano_plot$neg_log_padj_limits$max
        )
        updateNumericInput(
          session,
          'volcano_ymin',
          value = config()$ui$de_analysis$volcano_plot$neg_log_padj_limits$min
        )
      })

      # Syncs thresholds when plot_args() changes
      observeEvent(c(plot_args()$fdr.thres, plot_args()$fc.thres), {
        fc.thres <- ifelse(
          plot_args()$fc.thres == '' | is.na(plot_args()$fc.thres),
          config()$ui$de_analysis$filters$log2fc_threshold,
          plot_args()$fc.thres
        )
        fdr.thres <- ifelse(
          plot_args()$fdr.thres == '' | is.na(plot_args()$fdr.thres),
          config()$ui$de_analysis$filters$fdr_threshold,
          plot_args()$fdr.thres
        )

        curr_thres$fdr.thres <- fdr.thres
        curr_thres$fc.thres <- fc.thres
      })

      # Server logic for the autoscaler button
      observeEvent(input$volcano_auto, {
        showNotification('autoscaling axis limits')
        df <- app_object()$res[[input$comp_all]]
        filtered_lfc <- df$log2FoldChange[is.finite(df$log2FoldChange)]
        df.x_max <- round(max(filtered_lfc, na.rm = TRUE) * 1.05, digits = 3)
        df.x_min <- round(min(filtered_lfc, na.rm = TRUE) * 1.05, digits = 3)

        # Computes the values the plot will use for the
        # y-axis and updated the fields.
        log_padj <- -log10(df$padj)
        log_padj <- log_padj[is.finite(log_padj)]
        df.y_max <- round(max(log_padj, na.rm = TRUE) * 1.1, digits = 3)
        df.y_min <- round(min(log_padj, na.rm = TRUE) * 1.05, digits = 3)
        updateNumericInput(session, 'volcano_xmin', value = df.x_min)
        updateNumericInput(session, 'volcano_xmax', value = df.x_max)
        updateNumericInput(session, 'volcano_ymin', value = df.y_min)
        updateNumericInput(session, 'volcano_ymax', value = df.y_max)
      })

      # eventreactive for the static volcano plot
      volcano_plot <- eventReactive(
        c(
          input$comp_all,
          curr_thres$fdr.thres,
          curr_thres$fc.thres,
          curr_thres$colorscale,
          input$volcano_xmin,
          input$volcano_xmax,
          input$volcano_ymin,
          input$volcano_ymax,
          input$color_by,
          plot_args()$gene.to.plot,
          input$volcano_alpha
        ),
        {
          # Checks that the required inputs exist and are valid
          validate(
            need(
              !is.null(app_object()$res) &
                !is.null(input$comp_all) &
                input$comp_all != '',
              'waiting for selection'
            ),
            need(
              input$comp_all %in% names(app_object()$res),
              'selection not found in data'
            ),
            need(
              input$volcano_xmin != '' & input$volcano_xmax != '',
              'x-axis limits missing'
            ),
            need(
              input$volcano_xmin < input$volcano_xmax,
              'x-axis min must be less than max'
            ),
            need(
              input$volcano_ymin != '' & input$volcano_ymax != '',
              'y-axis limits missing'
            ),
            need(
              input$volcano_ymin < input$volcano_ymax,
              'y-axis min must be less than max'
            )
          )

          plotVolcano.label(
            app_object()$res[[input$comp_all]],
            fc.thres = curr_thres$fc.thres,
            fdr.thres = curr_thres$fdr.thres,
            neg_log_padj.lim = c(input$volcano_ymin, input$volcano_ymax),
            fc.lim = c(input$volcano_xmin, input$volcano_xmax),
            color_by = input$color_by,
            lab.genes = plot_args()$gene.to.plot,
            alpha = input$volcano_alpha
          )
        }
      )

      # eventreactive for the plot_ly
      volcano_plot_ly <- eventReactive(
        c(
          app_object()$res,
          input$comp_all,
          curr_thres$fdr.thres,
          curr_thres$fc.thres,
          curr_thres$colorscale,
          input$volcano_xmin,
          input$volcano_xmax,
          input$volcano_ymin,
          input$volcano_ymax,
          input$color_by,
          plot_args()$gene.to.plot,
          input$volcano_alpha
        ),
        {
          # Checks that the required inputs exist and are valid
          validate(
            need(
              !is.null(app_object()$res) &
                !is.null(input$comp_all) &
                input$comp_all != '',
              'waiting for selection'
            ),
            need(
              input$comp_all %in% names(app_object()$res),
              'selection not found in data'
            ),
            need(
              input$volcano_xmin != '' & input$volcano_xmax != '',
              'x-axis limits missing'
            ),
            need(
              input$volcano_xmin < input$volcano_xmax,
              'x-axis min must be < x-axis max'
            ),
            need(
              input$volcano_ymin != '' & input$volcano_ymax != '',
              'y-axis limits missing'
            ),
            need(
              input$volcano_ymin < input$volcano_ymax,
              'y-axis min must be < y-axis max'
            )
          )

          plotVolcano.label_ly(
            app_object()$res[[input$comp_all]],
            fc.thres = curr_thres$fc.thres,
            fdr.thres = curr_thres$fdr.thres,
            colorscale = curr_thres$colorscale,
            fc.lim = c(input$volcano_xmin, input$volcano_xmax),
            neg_log_padj.lim = c(input$volcano_ymin, input$volcano_ymax),
            color_by = input$color_by,
            lab.genes = plot_args()$gene.to.plot,
            alpha = input$volcano_alpha
          )
        }
      )

      output$volcano_plot_out <- renderUI({
        # Renders interactive plot
        if (input$plot_interactive == 'yes') {
          p <- volcano_plot_ly() %>% toWebGL()

          p <- plotly::plotly_build(p)

          output$plot1 <- renderPlotly({
            p
          })

          withSpinner(
            plotlyOutput(ns('plot1'), height = '600px')
          )
        } else if (input$plot_interactive == 'no') {
          # Renders non-interactive plot
          p <- volcano_plot() + theme(text = element_text(size = 18))
          output$plot2 <- renderPlot({
            p
          })
          withSpinner(
            plotOutput(ns('plot2'), height = '600px')
          )
        }
      })
      downloadButtonServer(
        'volcano_plot_download',
        volcano_plot,
        'volcano_plot'
      )
    }
  )
}
