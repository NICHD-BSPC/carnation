#' Scatterplot module
#'
#' @description
#' Module UI + server for generating scatter plots.
#'
#' @param id Module id
#' @param panel string, can be 'sidebar' or 'main' passed to UI
#' @param obj reactiveValues object containing carnation object passed to server
#' @param plot_args reactive containing 'fdr.thres' (padj threshold), 'fc.thres' (log2FC)
#' @param gene_scratchpad reactive containing gene scratchpad genes
#' @param reset_genes reactive to reset gene scratchpad selection
#' @param config reactive list with config settings passed to server
#'
#' @returns
#' UI returns tagList with scatter plot UI.
#' Server invisibly returns NULL (used for side effects).
#'
#' @examplesIf interactive()
#' library(shiny)
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
#' plot_args <- reactive({
#'   list(
#'     fdr.thres=0.1,
#'     fc.thres=0
#'   )
#' })
#'
#' gene_scratchpad <- reactive({ c('gene1', 'gene2') })
#' reset_genes <- reactiveVal()
#'
#' config <- reactiveVal(get_config())
#'
#' shinyApp(
#'   ui = fluidPage(
#'          sidebarPanel(scatterPlotUI('p', 'sidebar')),
#'          mainPanel(scatterPlotUI('p', 'sidebar'))
#'        ),
#'   server = function(input, output, session){
#'              scatter_data <- scatterPlotServer('p', obj, plot_args,
#'                                gene_scratchpad, reset_genes, config)
#'            }
#' )
#'
#' @name scattermod
#' @rdname scattermod
NULL

#' @rdname scattermod
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

      conditionalPanel(paste0('input["', ns('plot_interactive'), '"] == "yes"'),
        withSpinner(
          plotlyOutput(ns('plotly_out'), height='600px')
        )
      ), # conditionalPanel
      conditionalPanel(paste0('input["', ns('plot_interactive'), '"] == "no"'),
        withSpinner(
          plotOutput(ns('plot_out'), height='600px')
        )
      ), # conditionalPanel

      fluidRow(
        column(3,
          h4('Table filters'),
          selectizeInput(ns('filter_tbl'), h5('Filter by significance'),
                         choices=NULL, selected=NULL, multiple=TRUE),

          # For 'Select all' and 'Select none' buttons
          fluidRow(style='margin-bottom: 5px; margin-left: 2px;',
            actionButton(ns('select_all'), 'Select all', class = "btn-secondary"),
            actionButton(ns('select_none'), 'Select none', class = "btn-secondary")
          ), # fluidRow

          actionButton(ns('filter_tbl_do'),
                       style='margin-bottom: 20px;',
                       label='Apply filters',
                       icon=icon('filter'),
                       class='btn-primary'),

          h4('Table selection'),
          fluidRow(style='margin-left: 2px;',
            actionButton(ns('add_selected'), 'Add to scratchpad'),
            actionButton(ns('reset_tbl'),
                         'Reset selection',
                         class='btn-primary')
          ) # fluidRow
        ), # column
        column(9, style='margin-top: 10px;',
          withSpinner(
            DTOutput(ns('scatter_tbl'))
          ) # withSpinner
        ) # column
      ) # fluidRow

    ) # tagList
  } # else if panel='main'
} # scatterPlotUI

#' @rdname scattermod
#' @export
scatterPlotServer <- function(id, obj, plot_args, gene_scratchpad, reset_genes, config){

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

      # reactive to hold labeled genes
      genes_clicked <- reactiveValues(g=NULL)

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

      # gene scratchpad

      observeEvent(gene_scratchpad(), {
        g <- gene_scratchpad()
        if(any(g != '')){
          if(!all(g %in% genes_clicked$g))
            genes_clicked$g <- unique(c(genes_clicked$g, g))
        }
      })

      observeEvent(reset_genes(), {
        genes_clicked$g <- NULL
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
             curr_thres$fdr.thres,
             genes_clicked$g)
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
            warning(length(diff.genes), ' genes were discarded because found in one res but not the other')
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

        cols.sub <- c(input$compare, 'geneid')
        df <- dplyr::inner_join(
          dplyr::select(res_i, all_of(cols.sub)),
          dplyr::select(res_j, all_of(cols.sub)),
          by = 'geneid',
          suffix = c('.x', '.y')
        )

        compare <- input$compare

        x_column <- paste0(input$compare, '.x')
        y_column <- paste0(input$compare, '.y')
        df_temp <- df

        # Need to -log10 transform padj.x and padj.y to get proper limits
        if(input$compare == 'padj'){
          df_temp[[x_column]] <- -log10(df_temp[[x_column]])
          df_temp[[y_column]] <- -log10(df_temp[[y_column]])
        }

        # Fetch user inputs
        lim.x <- c(input$scatter_xmin, input$scatter_xmax)
        lim.y <- c(input$scatter_ymin, input$scatter_ymax)

        # update & handle NAs
        x_lfc_allna <- all(is.na(df_temp[, x_column]))

        validate(
          need(!x_lfc_allna,
               paste(input$x_axis_comp, compare, ' column has all NAs!',
                     'Please choose different contrast for x-axis'))
        )

        y_lfc_allna <- all(is.na(df_temp[, y_column]))

        validate(
          need(!y_lfc_allna,
               paste(input$y_axis_comp, compare, ' column has all NAs!',
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
        sig_levels <- c('None', label_i, label_j,
                        'Both - opposite LFC sign', 'Both - same LFC sign')
        df$significance <- factor(df$significance, levels = sig_levels)

        # Store the dataframe in the df_react reactiveVal
        updateSelectizeInput(session, 'filter_tbl', choices=sig_levels, selected=sig_levels)

        df_react(df)
        df_full(df_full)

        flags$data_loaded <- flags$data_loaded + 1
      })

      # observers for tbl filters
      observeEvent(input$select_all, {
        # Get all possible comparison options
        all_sig <- levels(df_react()$significance)
        # Update the select_none checkbox
        updateSelectizeInput(session, 'filter_tbl', selected=all_sig)
      })

      # Observer for Select none checkbox
      observeEvent(input$select_none, {
        # Update comp_all with no selected comparisons
        updateSelectizeInput(session, 'filter_tbl', selected=character(0))
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
        validate(
          need(!is.na(input$scatter_xmin), '')
        )
        axis_limits$lim.x[1] <- input$scatter_xmin
      })

      observeEvent(input$scatter_xmax, {
        validate(
          need(!is.na(input$scatter_xmax), '')
        )
        axis_limits$lim.x[2] <- input$scatter_xmax
      })

      observeEvent(input$scatter_ymin, {
        validate(
          need(!is.na(input$scatter_ymin), '')
        )
        axis_limits$lim.y[1] <- input$scatter_ymin
      })

      observeEvent(input$scatter_ymax, {
        validate(
          need(!is.na(input$scatter_ymax), '')
        )
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
        df <- df %>% filter(!is.na(.data[[ xcol ]]), !is.na(.data[[ ycol ]]))

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
        genes <- genes_clicked$g
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

      output$plotly_out <- renderPlotly({
        scatterplot_ly()
      })

      output$plot_out <- renderPlot({
        scatterplot() + theme(text=element_text(size=18))
      })

      # ------------------------------------------------------- #

      # -------------------- renerDataTable  ------------------ #
      scatter_dt <- eventReactive(c(df_react(), input$x_axis_comp, input$y_axis_comp,
                                    input$filter_tbl_do), {
        validate(
          need(!is.null(df_full()), '')
        )
        df <- df_full()

        # add significance column
        sig_df <- df_react()

        # map geneid to significance
        sigvec <- sig_df$significance
        names(sigvec) <- sig_df$geneid
        df$significance <- sigvec[df$geneid]

        # get current filters
        curr_filters <- input$filter_tbl
        idx <- df$significance %in% curr_filters
        if(sum(idx) == 0){
          showNotification(
            'No rows left in table after filtering', type='error'
          )

          validate(need(sum(idx) > 0, ''))
        }
        # move geneid & significance to beginning to work with container
        df[which(idx),] %>% relocate('significance') %>% relocate('geneid')

      })

      # Optionally display datatable for scatter plot input data
      output$scatter_tbl <- renderDT({
        validate(
          need(!is.null(scatter_dt()), '')
        )
        df <- scatter_dt()

        # Define the columns to format to 3 sig figs
        columns_to_format <- c("padj.x", "padj.y", "log2FoldChange.x", "log2FoldChange.y")
        which_cols <- which(colnames(df) %in% columns_to_format)
        border_cols <- c(1, 2, grep('padj', colnames(df)))

        all_comps <- c(input$x_axis_comp, input$y_axis_comp)
        validate(
          need(!any(is.null(all_comps)), 'Waiting for selection')
        )

        # build container for table
        sketch <- htmltools::withTags(table(
          class = 'display',
          tags$thead(
            tags$tr(
              lapply(c('geneid', 'significance'), function(x) tags$th(rowspan=2, x)),
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
                    selection=list(mode='multiple'),
                    container=sketch,
                    options=list(autoWidth=TRUE,
                                 columnDefs=list(list(className='dt-center',
                                                      targets=seq_len((ncol(df)-1)))))) %>%
          formatStyle(columns=border_cols,
                      'border-right'='solid 1px') %>%
          formatSignif(columns=which_cols, digits=5)
      }) # renderDT
      # ------------------------------------------------------- #

      # table selection handling

      scatter_proxy <- dataTableProxy('scatter_tbl')

      observeEvent(input$reset_tbl, {
        scatter_proxy %>% selectRows(NULL)
      })

      observeEvent(input$add_selected, {
        tbl <- scatter_dt()
        sel <- input$scatter_tbl_rows_selected

        # handle NAs in symbol
        s <- tbl$geneid
        #s[is.na(s)] <- tbl$gene[is.na(s)]

        if(is.null(sel)){
          showNotification(
            'Cannot add genes, no rows selected', type='warning'
          )
          validate(
            need(!is.null(sel), '')
          )
        } else if(all(s[sel] %in% genes_clicked$g)){
          showNotification(
            'Selected genes already present in scratchpad, skipping', type='warning'
          )
        } else {
          if(is.null(genes_clicked$g)) selected <- s[sel]
          else selected <- unique(c(genes_clicked$g, s[sel]))

          new_genes <- setdiff(selected, genes_clicked$g)
          showNotification(
            paste('Adding', length(new_genes),
                  'new genes to scratchpad')
          )

          genes_clicked$g <- c(genes_clicked$g, new_genes)
        }
      }) # observeEvent

      # ---------------- Help and download buttons ------------ #
      helpButtonServer('de_cmp_scatter_help', size='l')
      helpButtonServer('de_scatter_help', size='l')
      downloadButtonServer('scatterplot_download', scatterplot, 'scatterplot')
      # ----------------------------------------------------- #

      return(
        reactive({
          list(genes=genes_clicked$g)
        })
      )
    } # Server function
  ) # moduleServer
} # maPlotServer
