#' Scatterplot module UI
#'
#' @param id ID string used to match the ID used to call the module server function
#' @param panel string, can be 'sidebar' or 'main'
#'
#' @export
scatterPlotUI <- function(id, panel){
  ns <- NS(id)

  if(panel == 'sidebar'){
    tagList(
      fluidRow(
        column(12, align='right',
          helpButtonUI(ns('de_cmp_scatter_help'))
        ) # column
      ), # fluidRow


      selectizeInput(ns('x_axis_comp'),
                    label='Comparison 1 (x-axis)',
                    choices=NULL,
                    selected=NULL
      ), # selectizeInput

      selectizeInput(ns('y_axis_comp'),
                    label='Comparison 2 (y-axis)',
                    choices=NULL,
                    selected=NULL
      ), # selectizeInput

      fluidRow(
        column(12, align='center',
          actionButton(ns('swap_comp'), label='Swap comparisons',
            icon=icon('arrows-rotate'),
              status='info',
              style='margin-bottom: 20px;')
        ) # column
      ), # fluidRow

      fluidRow(
        column(6, h5('Values to use')),
        column(6,
          selectizeInput(ns('compare'),
                        label=NULL,
                        choices=c('LFC'='log2FoldChange', 'P-adj'='padj')
          ) # selectizeInput
        ) # column
      ), # fluidRow

      fluidRow(
        column(6, h5('Interactive?')),
        column(6,
          selectInput(ns("plot_interactive"), label=NULL,
                      choices=c('yes', 'no'),
                      selected='yes'
          ) # selectInput
        ) # column
      ), # fluidRow

      fluidRow(
        column(6, h5('Show table?')),
        column(6,
          selectInput(ns("show_table"), label=NULL,
                      choices=c('yes', 'no'),
                      selected='yes'
          ) # selectInput
        ) # column
      ), # fluidRow

      wellPanel(
        style='background: white',

        h4('Plot settings', style='margin-bottom: 10px; font-weight: bold;'),

        bsCollapse(
          bsCollapsePanel('Axes limits',

            tags$label(class='control-label',
                       'x-axis limits'),

            fluidRow(
              column(6, h5('max')),
              column(6,
                numericInput(ns('scatter_xmax'), label=NULL,
                  value=NULL)
              ) # column
            ), # fluidRow

            fluidRow(
              column(6, h5('min')),
              column(6,
                numericInput(ns('scatter_xmin'), label=NULL,
                  value=NULL)
              ) # column
            ), # fluidRow

            fluidRow(
              column(4, align='left', style='margin-bottom: 10px;',
                actionButton(ns('scatter_x_auto'), label='Autoscale')
              ) # column
            ), # fluidRow

            tags$label(class='control-label',
                       'y-axis limits'),

            fluidRow(
              column(6, h5('max')),
              column(6,
                numericInput(ns('scatter_ymax'), label=NULL,
                  value=NULL)
              ) # column
            ), # fluidRow

            fluidRow(
              column(6, h5('min')),
              column(6,
                numericInput(ns('scatter_ymin'), label=NULL,
                  value=NULL)
              ) # column
            ), # fluidRow

            fluidRow(
              column(4, align='left', style='margin-bottom: 10px;',
                actionButton(ns('scatter_y_auto'), label='Autoscale')
              ) # column
            ) # fluidRow
          ) # bsCollapsePanel
        ), # bsCollapse

        bsCollapse(
          bsCollapsePanel('Point aesthetics',

            fluidRow(
              column(6, h5('Color palette')),
              column(6,
                selectInput(ns('color.palette'), label=NULL,
                            choices=c('Accent', 'Dark2', 'Paired', 'Pastel1', 'Pastel2', 'Set1', 'Set2', 'Set3'),
                            selected='Set2')
              ) # column
            ), # fluidRow

            fluidRow(
              column(6, h5('Marker opacity')),
              column(6,
                numericInput(ns("alpha"), label=NULL,
                  value=0.7,
                  min=0,
                  max=1,
                  step=0.1
                ) # numericInput
              ) # column
            ), # fluidRow

            fluidRow(
              column(6, h5('Marker size')),
              column(6,
                numericInput(ns("size"), label=NULL,
                  value=4,
                  min=0,
                  max=10,
                  step=0.1
                ) # numericInput
              ) # column
            ) # fluidRow

          ) # bsCollapsePanel
        ), # bsCollapse

        bsCollapse(
          bsCollapsePanel('Grid lines',

            fluidRow(
              column(6, h5('Show x=0?')),
              column(6,
                selectInput(ns("vline"), label=NULL,
                  choices=c('yes', 'no'),
                  selected='yes'
                ) # selectInput
              ) # column
            ), # fluidRow

            fluidRow(
              column(6, h5('Show y=0?')),
              column(6,
                selectInput(ns("hline"), label=NULL,
                            choices=c('yes', 'no'),
                            selected='yes'
                ) # selectInput
              ) # column
            ), # fluidRow

            fluidRow(
              column(6, h5('Show grid?')),
              column(6,
                selectInput(ns("show_grid"), label=NULL,
                            choices=c('yes', 'no'),
                            selected='yes'
                ) # selectInput
              ) # column
            ), # fluidRow

            fluidRow(
              column(6, h5('Show diagonal?')),
              column(6,
                selectInput(ns("dline"), label=NULL,
                            choices=c('yes', 'no'),
                            selected='yes'
                ) # selectInput
              ) # column
            ) # fluidRow

          ) # bsCollapsePanel
        ), # bsCollapse

        fluidRow(align='center',
          actionButton(ns('refresh'),
                       label='Refresh plot',
                       icon=icon('arrows-rotate'),
                       class='btn-primary',
                       style='margin-bottom: 10px;')
        ) # fluidRow
      ) # wellPanel

    ) # tagList

  } else if (panel == 'main') {
    tagList(
      fluidRow(
        column(6, align='left',
          helpButtonUI(ns('de_scatter_help'))
        ), # column
        column(6, align='right',
          downloadButtonUI(ns('scatterplot_download'))
        ) # column
      ), # fluidRow

      withSpinner(
        uiOutput(ns('scatterplot_out'))
      ), # withSpinner

      withSpinner(
        DTOutput(ns('scatter_datatable_out'))
      ) # withSpinner

     # conditionalPanel("input.show_table == 'yes'",
     # ) # conditionalPanel
    ) # tagList
  } # else if panel='main'
} # scatterPlotUI

#' Scatterplot module server function
#'
#' @param id ID string used to match the ID used to call the module UI function
#' @param obj reactiveValues object containing carnation object
#' @param plot_args reactive containing 'fdr.thres' (padj threshold), 'fc.thres' (log2FC) &
#' 'gene.to.plot' (genes to be labeled)
#' @param config reactive list with config settings
#'
#' @export
scatterPlotServer <- function(id, obj, plot_args, config){

  moduleServer(
    id,

    function(input, output, session){
      # -- Set ns, Get config, Load data, Trigger data laoded flag -- #
      ns <- NS(id)

      app_object <- reactive({
        list(res=obj$res)
      })

      # Reactive list of all available comparisons
      comp_all <- reactive({
        names(app_object()$res)
      })

      # reactive values used to track/store different settings
      flags <- reactiveValues(data_loaded=0)

      # current FDR/fold-change thresholds
      curr_thres <- reactiveValues(fdr.thres=0.1,
                                   fc.thres=0)

      # Reactive df that will be used for autoscaling and plotting
      # df_react: data frame used for plotting
      # df_full: full data frame
      df_react <- reactiveVal(NULL)
      df_full <- reactiveVal(NULL)

      # reactive values to keep track of axis limits
      axis_limits <- reactiveValues(lim.x=NULL, lim.y=NULL)

      # Initialize comparison selection dropdowns when data is available
      observeEvent(comp_all(), {
        validate(
          need(!is.null(app_object()$res), 'Waiting for data')
        )
        # Populate the x and y-axis dropdowns
        updateSelectInput(session, 'x_axis_comp', choices = comp_all(), selected = comp_all()[1])
        if (length(comp_all()) == 1) {
          available_y <- comp_all()[1]
        } else {
          available_y <- comp_all()[2]
        }
        updateSelectizeInput(session, 'y_axis_comp', choices = comp_all(), selected = available_y)

        df_react(NULL)
        df_full(NULL)
        })
      # -------------------------------------------------------------- #

      # --------------- Set FDR and FC thresholds ---------------- #

      # update from reactive config
      observeEvent(config(), {
        curr_thres$fdr.thres <- config()$ui$de_analysis$filters$fdr_threshold
        curr_thres$fc.thres <- config()$ui$de_analysis$filters$log2fc_threshold
      })

      observeEvent(c(plot_args()$fdr.thres, plot_args()$fc.thres), {
        fc.thres <- ifelse(plot_args()$fc.thres == '' | is.na(plot_args()$fc.thres), 0, plot_args()$fc.thres)
        fdr.thres <- ifelse(plot_args()$fdr.thres == '' | is.na(plot_args()$fdr.thres), 0.1, plot_args()$fdr.thres)

        curr_thres$fdr.thres <- fdr.thres
        curr_thres$fc.thres <- fc.thres
      })
      # ------------------------------------------------------------#

      # ------------- helper functions -----------------------------#

      # function to get range of column in df after handling NAs
      get_range <- function(df, column, lim){
        lim[1] <- min(df[[ column ]], na.rm=TRUE)
        lim[2] <- max(df[[ column ]], na.rm=TRUE)

        diff <- diff(lim) * 0.05

        lim[1] <- round(lim[1] - diff, digits=2)
        lim[2] <- round(lim[2] + diff, digits=2)

        return(lim)
      }

      update_geneid <- function(df) {
        na_idx <- is.na(df$symbol)
        tbl <- table(df$symbol)
        dups <- names(tbl)[tbl>1]
        duplicates_idx <- df$symbol %in% dups
        symbols_to_fix <- na_idx | duplicates_idx
        symbols_ok <- !symbols_to_fix
        df$geneid[symbols_ok] <- df$symbol[symbols_ok]
        return(df)
      }

      # Function to get autoscale limits for both plot types and both axes
      autoscale <- function(df, compare, lim.x, lim.y) {
        # Determine columns based on compare parameter
        x_column <- paste0(compare, '.x')
        y_column <- paste0(compare, '.y')
        df_temp <- df
        if (compare == 'padj') {
          df_temp[[x_column]] <- -log10(df_temp[[x_column]])
          df_temp[[y_column]] <- -log10(df_temp[[y_column]])
        }

        # Set axis limits
        lim.x <- get_range(df_temp, x_column, lim.x)
        lim.y <- get_range(df_temp, y_column, lim.y)

        return(list(lim.x = lim.x, lim.y = lim.y))
      }
      # ---------------------------------------------------------- #


      # --------- Generate universal scatter data frame ---------- #
 
      observeEvent({
        list(input$compare,
             input$x_axis_comp,
             input$y_axis_comp,
             curr_thres$fc.thres,
             curr_thres$fdr.thres)
      }, {

        req(app_object()$res)

        res_i <- as.data.frame(app_object()$res[[input$x_axis_comp]])
        res_j <- as.data.frame(app_object()$res[[input$y_axis_comp]])

        # Check for matching genes in res_i and res_j
        diff.genes <- unique(unlist(c(
            setdiff(rownames(res_i), rownames(res_j)),
            setdiff(rownames(res_j), rownames(res_i))
        )))

        if (length(diff.genes) > 0) {
            warning(paste0(length(diff.genes), ' genes were discarded because found in one res but not the other'))
        }

        # Lowercase specified column names in both dataframes
        name.col <- 'SYMBOL'
        colnames(res_i)[colnames(res_i) %in% name.col] <- tolower(name.col)
        colnames(res_j)[colnames(res_j) %in% name.col] <- tolower(name.col)

        # Initialize geneid column to join by
        res_i$geneid <- rownames(res_i)
        res_j$geneid <- rownames(res_j)

        # Prepare geneid column in both dfs to prepare for joining
        res_i <- update_geneid(res_i)
        res_j <- update_geneid(res_j)

        # Make a temp df that will be used to determine significance in the single column df
          # required for df_react() to react to changed in input$compare
          cols.sub <- c('log2FoldChange', 'padj', 'geneid')
          df_full <- dplyr::inner_join(
            dplyr::select(res_i, all_of(cols.sub)),
            dplyr::select(res_j, all_of(cols.sub)),
            by = 'geneid',
            suffix = c('.x', '.y')
          )

        if (input$compare == 'LFC') {
          # Join to make df for df_react(). df_react() must change in response
          # to input$compare so we make a distinct df for each input$compare choice
          cols.sub <- c('log2FoldChange', 'geneid')
          df <- dplyr::inner_join(
            dplyr::select(res_i, all_of(cols.sub)),
            dplyr::select(res_j, all_of(cols.sub)),
            by = 'geneid',
            suffix = c('.x', '.y')
          )
        } else if (input$compare == 'P-adj') {
          cols.sub <- c('padj', 'geneid')
          df <- dplyr::inner_join(
            dplyr::select(res_i, all_of(cols.sub)),
            dplyr::select(res_j, all_of(cols.sub)),
            by = 'geneid',
            suffix = c('.x', '.y')
          )
        }

        compare <- input$compare

        # Determine columns based on compare parameter
        if (compare == 'LFC') {
          x_column <- 'log2FoldChange.x'
          y_column <- 'log2FoldChange.y'
          df_temp <- df
        } else if (compare == 'P-adj') {
          x_column <- 'padj.x'
          y_column <- 'padj.y'
          # Need to -log10 transform padj.x and padj.y to get proper limits
          df_temp <- df
          df_temp[[x_column]] <- -log10(df_temp[[x_column]])
          df_temp[[y_column]] <- -log10(df_temp[[y_column]])
        }

        # Fetch user inputs
        lim.x <- c(input$scatter_xmin, input$scatter_xmax)
        lim.y <- c(input$scatter_ymin, input$scatter_ymax)

        # update & handle NAs
        x_lfc_allna <- all(is.na(df_temp[, x_column]))
        if(x_lfc_allna){
          showNotification(
            paste('Scatter-plot warning:', input$x_axis_comp,
                  'LFC column has all NAs!',
                  'Please choose different contrast for x-axis'),
            duration=15
          )
        }

        validate(
          need(!x_lfc_allna,
               paste(input$x_axis_comp, 'log2FoldChange column has all NAs!',
                     'Please choose different contrast for x-axis'))
        )

        y_lfc_allna <- all(is.na(df_temp[, y_column]))
        if(y_lfc_allna){
          showNotification(
            paste('Scatter-plot warning:', input$y_axis_comp,
                  'log2FoldChange column has all NAs!',
                  'Please choose different contrast for x-axis'),
            duration=15
          )
        }

        validate(
          need(!y_lfc_allna,
               paste(input$y_axis_comp, 'log2FoldChange column has all NAs!',
                     'Please choose different contrast for x-axis'))
        )
        lim.x <- get_range(df_temp, x_column, lim.x)
        lim.y <- get_range(df_temp, y_column, lim.y)

        # save to reactive values & updt inputs
        axis_limits$lim.x <- lim.x
        axis_limits$lim.y <- lim.y

        updateNumericInput(session, 'scatter_xmin', value=lim.x[1])
        updateNumericInput(session, 'scatter_xmax', value=lim.x[2])
        updateNumericInput(session, 'scatter_ymin', value=lim.y[1])
        updateNumericInput(session, 'scatter_ymax', value=lim.y[2])

        # Fetch  remaining user inputs
        label_i <- input$x_axis_comp
        label_j <- input$y_axis_comp
        fdr.thres <- curr_thres$fdr.thres
        fc.thres <- curr_thres$fc.thres

        # Add significance column in df using df_color values
        df <- df %>%
          mutate(significance = case_when(
            df_full$padj.x <= fdr.thres &
            df_full$padj.y <= fdr.thres &
            (df_full$log2FoldChange.x * df_full$log2FoldChange.y >= 0) &
            abs(df_full$log2FoldChange.x) >= fc.thres &
            abs(df_full$log2FoldChange.y) >= fc.thres ~ 'Both - same LFC sign',

            df_full$padj.x <= fdr.thres &
            df_full$padj.y <= fdr.thres &
            (df_full$log2FoldChange.x * df_full$log2FoldChange.y < 0) &
            abs(df_full$log2FoldChange.x) >= fc.thres &
            abs(df_full$log2FoldChange.y) >= fc.thres  ~ 'Both - opposite LFC sign',

            (df_full$padj.x <= fdr.thres & abs(df_full$log2FoldChange.x) >= fc.thres) &
            (df_full$padj.y > fdr.thres | abs(df_full$log2FoldChange.y) < fc.thres) ~ label_i,

            (df_full$padj.y <= fdr.thres & abs(df_full$log2FoldChange.y) >= fc.thres) &
            (df_full$padj.x > fdr.thres | abs(df_full$log2FoldChange.x) < fc.thres) ~ label_j,

            TRUE ~ 'None'))

        if (label_i == label_j) {
          label_i <- paste0(label_i, '_x')
          label_j <- paste0(label_j, '_y')
        }

        # Significance as factor, to reorder in the graph
        df$significance <- factor(df$significance, levels = c('None', label_i, label_j, 'Both - opposite LFC sign', 'Both - same LFC sign'))

        # Store the dataframe in the df_react reactiveVal
        df_react(df)
        df_full(df_full)

        flags$data_loaded <- flags$data_loaded + 1
      })
      # ---------------------------------------------------------- #

      # --------------- Swap comparisons button  ----------------- #
      observeEvent(input$swap_comp, {
        updateSelectInput(session, 'y_axis_comp',
                          selected=input$x_axis_comp)
        updateSelectInput(session, 'x_axis_comp',
                          selected=input$y_axis_comp)
      })
      # ---------------------------------------------------------- #

      # ------------------ Observer: Axis limits ----------------- #
      observeEvent(input$scatter_xmin, {
        axis_limits$lim.x[1] <- input$scatter_xmin
      })

      observeEvent(input$scatter_xmax, {
        axis_limits$lim.x[2] <- input$scatter_xmax
      })

      observeEvent(input$scatter_ymin, {
        axis_limits$lim.y[1] <- input$scatter_ymin
      })

      observeEvent(input$scatter_ymax, {
        axis_limits$lim.y[2] <- input$scatter_ymax
      })

      # ---------------------------------------------------------- #

      # --------------- Observer: Autoscale x button ------------- #
      observeEvent(list(input$scatter_x_auto, df_react()), {
        validate(
          need(!is.null(df_react()), "Waiting for selection")
        )

        # Get x lims
        lims <- autoscale(df=df_react(), compare=input$compare, lim.x=numeric(2), lim.y=numeric(2))
        lim.x <- lims[['lim.x']]
        updateNumericInput(session, 'scatter_xmin', value=lim.x[1])
        updateNumericInput(session, 'scatter_xmax', value=lim.x[2])
      }) # observeEvent
      # ---------------------------------------------------------- #

      # --------------- Observer: Autoscale y button ------------- #
      observeEvent(list(input$scatter_y_auto, df_react()), {
        validate(
          need(!is.null(df_react()), "Waiting for selection")
        )

        # Get y lims
        lims <- autoscale(df=df_react(), compare=input$compare, lim.x=numeric(2), lim.y=numeric(2))
        lim.y <- lims[['lim.y']]
        updateNumericInput(session, 'scatter_ymin', value=lim.y[1])
        updateNumericInput(session, 'scatter_ymax', value=lim.y[2])
      }) # observeEvent
      # --------------------------------------------------------- #

      # -------------------- Plot preparations ------------------ #
      plot_prep <- function() {
        df <- df_react()
        compare <- input$compare
        lim.x <- axis_limits$lim.x
        lim.y <- axis_limits$lim.y

        xcol <- paste0(compare, '.x')
        ycol <- paste0(compare, '.y')

        if (compare=='padj') {
          # Convert padj to -log10(padj) for x and y
          df$padj.x <- -log10(df$padj.x)
          df$padj.y <- -log10(df$padj.y)
        }

        # filter rows with NA values
        df <- df %>% filter(!is.na(.data[[ xcol ]]))
        df <- df %>% filter(!is.na(.data[[ xcol ]]))

        # Create column with plotting character based on lim.x
        # Change point values for those outside plot limits to values that are within the limits
        df$shape <- 'in'

        # count number of points out of axis limits
        num_ob_points <- sum(df[[ xcol ]] > lim.x[2] |
                             df[[ xcol ]] < lim.x[1] |
                             df[[ ycol ]] > lim.y[2] |
                             df[[ ycol ]] < lim.y[1])
        if(!is.na(num_ob_points)){
          if(num_ob_points > 0){
            showNotification(
              paste0('Warning: ', num_ob_points, ' points outside plot axes limits. ',
                     'These are being shown along the respective boundary.'),
              type='warning'
            )
          }
        }

        # replace x values outside limits with axes limits
        # build shape column to use for different plotting symbols
        df[[ xcol ]][ df[[ xcol ]] > lim.x[2] ] <- lim.x[2]
        df[[ xcol ]][ df[[ xcol ]] < lim.x[1] ] <- lim.x[1]
        df[[ 'shape' ]][ df[[ xcol ]] == lim.x[2] ] <- 'right'
        df[[ 'shape' ]][ df[[ xcol ]] == lim.x[1] ] <- 'left'

        # same for y
        df[[ ycol ]][ df[[ ycol ]] > lim.y[2] ] <- lim.y[2]
        df[[ ycol ]][ df[[ ycol ]] < lim.y[1] ] <- lim.y[1]
        df[[ 'shape' ]][ df[[ ycol ]] == lim.y[2] ] <- 'above'
        df[[ 'shape' ]][ df[[ ycol ]] == lim.y[1] ] <- 'below'
        df[[ 'shape' ]] <- as.factor(df[[ 'shape' ]])

        validate(
          if (!is.na(input$scatter_xmin) && !is.na(input$scatter_xmax)) {
            need(input$scatter_xmin < input$scatter_xmax,
                 'x-axis min must be < x-axis max')
          }
        )

        validate(
          if (!is.na(input$scatter_ymin) && !is.na(input$scatter_ymax)) {
            need(input$scatter_ymin < input$scatter_ymax,
                 'y-axis min must be < y-axis max')
          }
        )

        # Get genes to label
        genes <- plot_args()$gene.to.plot
        if(is.null(genes) || all(genes %in% '')){
            lab.genes <- NULL
        } else {
            lab.genes <- genes
        }

        color.palette <- brewer.pal(n = 5, name=input$color.palette)

        return(list(df=df, lab.genes=lab.genes, color.palette=color.palette,
                    lim.x=axis_limits$lim.x,
                    lim.y=axis_limits$lim.y))
      }
      # --------------------------------------------------------- #

      # -------------- Generate static scatter plot  ------------ #
      # And show OB DE points warning
      # reactive to generate labeled scatter plot
      scatterplot <- eventReactive(c(input$refresh, flags$data_loaded), {

        # Validation
        needs <- c(!is.null(df_react()), !is.null(df_full()))
        for (need in needs) {
          validate(need(need, "Waiting for selection"))
        }

        params <- plot_prep()
        df <- params[['df']]
        lab.genes <- params[['lab.genes']]
        color.palette <- params[['color.palette']]

        # Return scatterplot
        plotScatter.label(
          compare=input$compare,
          df=df,
          label_x=input$x_axis_comp,
          label_y=input$y_axis_comp,
          lab.genes=lab.genes,
          lim.x=params[['lim.x']],
          lim.y=params[['lim.y']],
          plot_all='yes',
          name.col='geneid',
          lines=c(input$vline, input$hline, input$dline),
          alpha=input$alpha,
          size=input$size,
          show.grid=input$show_grid,
          color.palette=color.palette
        )

      }) # eventReactive
      # ----------------------------------------------------- #

      # ------------ Generate plotly scatter plot ----------- #
      scatterplot_ly <- eventReactive(c(input$refresh, flags$data_loaded), {

        # Validation
        needs <- c(!is.null(df_react()), !is.null(df_full()))
        for (need in needs) {
          validate(need(need, "Waiting for selection"))
        }

        params <- plot_prep()
        df <- params[['df']]
        lab.genes <- params[['lab.genes']]
        color.palette <- params[['color.palette']]

        plotScatter.label_ly(
          compare=input$compare,
          df=df,
          label_x=input$x_axis_comp,
          label_y=input$y_axis_comp,
          lim.x=params[['lim.x']],
          lim.y=params[['lim.y']],
          name.col='geneid',
          lines=c(input$vline, input$hline, input$dline),
          alpha=input$alpha,
          size=input$size,
          show.grid=input$show_grid,
          color.palette=color.palette,
          lab.genes=lab.genes
        )

     }) # eventReactive scatterplot_ly
      # ----------------------------------------------------- #

      # ---------------------- renerUI ---------------------- #
      output$scatterplot_out <- renderUI({
        if (input$plot_interactive == 'yes') {

          p <- scatterplot_ly() %>% toWebGL()

          output$plot1 <- renderPlotly({ p })

          withSpinner(
            plotlyOutput(ns('plot1'), height='600px')
          )

        } else if (input$plot_interactive == 'no') {

          p <- scatterplot() + theme(text=element_text(size=18))

          output$plot2 <- renderPlot({ p })

          withSpinner(
            plotOutput(ns('plot2'), height='600px')
          )
        }
      }) # renderUI
      # ------------------------------------------------------- #

      # -------------------- renerDataTable  ------------------ #
      # Optionally display datatable for scatter plot input data
      output$scatter_datatable_out <- renderDT({
        if (input$show_table == 'yes') {
          validate(
            need(!is.null(df_full()), '')
          )
          df <- df_full()

          # move geneid to beginning to work with container
          df <- df %>% relocate('geneid')

          # Define the columns to format to 3 sig figs
          columns_to_format <- c("padj.x", "padj.y", "log2FoldChange.x", "log2FoldChange.y")
          which_cols <- which(colnames(df) %in% columns_to_format)
          border_cols <- c(1, grep('padj', colnames(df)))

          all_comps <- c(input$x_axis_comp, input$y_axis_comp)
          validate(
            need(!any(is.null(all_comps)), 'Waiting for selection')
          )

          # build container for table
          sketch <- htmltools::withTags(table(
            class = 'display',
            tags$thead(
              tags$tr(
                tags$th(rowspan=2, 'geneid'),
                lapply(all_comps,
                       function(x) tags$th(class='dt-center', colspan=2, x))
              ),
              tags$tr(
                lapply(rep(c('log2FoldChange', 'padj'), 2), tags$th)
              )
            )
          ))

          df %>%
            datatable(rownames=FALSE,
                      selection='none',
                      container=sketch,
                      options=list(autoWidth=TRUE,
                                   columnDefs=list(list(className='dt-center',
                                                        targets=1:(ncol(df)-1))))) %>%
            formatStyle(columns=border_cols,
                        'border-right'='solid 1px') %>%
            formatSignif(columns=which_cols, digits=5)
        } # if show_table == 'yes'
      }) # renderDT
      # ------------------------------------------------------- #

      # ---------------- Help and download buttons ------------ #
      helpButtonServer('de_cmp_scatter_help', size='l')
      helpButtonServer('de_scatter_help', size='l')
      downloadButtonServer('scatterplot_download', scatterplot, 'scatterplot')
      # ----------------------------------------------------- #
    } # Server function
  ) # moduleServer
} # maPlotServer
