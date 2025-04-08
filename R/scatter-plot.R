#' Scatterplot module UI
#'
#' @param id ID string used to match the ID used to call the module server function
#' @param panel string, can be 'sidebar' or 'main'
#'
#' @export
scatterPlotUI <- function(id, panel){
  ns <- NS(id)

  config <- get_config()

  if(panel == 'sidebar'){
    tagList(
      fluidRow(
        column(12, align='right',
          helpButtonUI(ns('de_cmp_scatter_help'))
        )
      ),

      selectizeInput(ns('x_axis_comp'),
        label='Comparison 1 (x-axis)',
        choices=NULL,
        selected=NULL
      ),

      selectizeInput(ns('y_axis_comp'),
        label='Comparison 2 (y-axis)',
        choices=NULL,
        selected=NULL
      ),

      fluidRow(
        column(12, align='center',
          actionButton(ns('swap_comp'), label='Swap comparisons',
            icon=icon('arrows-rotate'),
            status='info',
            style='margin-bottom: 20px;')
        )
      ),

      fluidRow(
        column(6, h5('Values to use')),
        column(6,
          selectizeInput(ns('compare'),
            label=NULL,
            choices=c('LFC', 'P-adj'),
            selected=config$ui$de_analysis$scatter_plot$compare
          )
        )
      ),

      fluidRow(
        column(6, h5('Interactive?')),
        column(6,
          selectInput(ns("plot_interactive"), label=NULL,
            choices=c('yes', 'no'),
            selected='yes'
          )
        )
      ),

      fluidRow(
        column(6, h5('Show table?')),
        column(6,
          selectInput(ns("show_table"), label=NULL,
            choices=c('yes', 'no'),
            selected='yes'
          )
        )
      ),

      wellPanel(
        style='background: white',

        h4('Plot settings', style='margin-bottom: 10px; font-weight: bold;'),

        bsCollapse(
          bsCollapsePanel('Axes limits',

            tags$label(class='control-label', 'x-axis limits'),
            fluidRow(
              column(6, h5('max')),
              column(6,
                numericInput(ns('scatter_xmax'), label=NULL,
                  value=NULL)
              )
            ),

            fluidRow(
              column(6, h5('min')),
              column(6,
                numericInput(ns('scatter_xmin'), label=NULL,
                  value=NULL)
              )
            ),

            fluidRow(
              column(4, align='left', style='margin-bottom: 10px;',
                actionButton(ns('scatter_x_auto'), label='Autoscale')
              )
            ),

            tags$label(class='control-label', 'y-axis limits'),
            fluidRow(
              column(6, h5('max')),
              column(6,
                numericInput(ns('scatter_ymax'), label=NULL,
                  value=NULL)
              )
            ),

            fluidRow(
              column(6, h5('min')),
              column(6,
                numericInput(ns('scatter_ymin'), label=NULL,
                  value=NULL)
              )
            ),

            fluidRow(
              column(4, align='left', style='margin-bottom: 10px;',
                actionButton(ns('scatter_y_auto'), label='Autoscale')
              )
            )
          )
        ),

        bsCollapse(
          bsCollapsePanel('Point aesthetics',

            fluidRow(
              column(6, h5('Color palette')),
              column(6,
                selectInput(ns('color.palette'), label=NULL,
                  choices=c('Accent', 'Dark2', 'Paired', 'Pastel1', 'Pastel2', 'Set1', 'Set2', 'Set3'),
                  selected='Set2')
              )
            ),

            fluidRow(
              column(6, h5('Marker opacity')),
              column(6,
                numericInput(ns("alpha"), label=NULL,
                  value=1, min=0, max=1, step=0.1
                )
              )
            ),

            fluidRow(
              column(6, h5('Marker size')),
              column(6,
                numericInput(ns("size"), label=NULL,
                  value=4, min=0, max=10, step=0.1
                )
              )
            ),

            fluidRow(
              column(6, h5('Show all points?')),
              column(6,
                selectInput(ns("plot_all"), label=NULL,
                  choices=c('yes', 'no'),
                  selected='no'
                )
              )
            )
          )
        ),

        bsCollapse(
          bsCollapsePanel('Grid lines',

            fluidRow(
              column(6, h5('Show x=0?')),
              column(6,
                selectInput(ns("vline"), label=NULL,
                  choices=c('yes', 'no'),
                  selected='yes'
                )
              )
            ),

            fluidRow(
              column(6, h5('Show y=0?')),
              column(6,
                selectInput(ns("hline"), label=NULL,
                  choices=c('yes', 'no'),
                  selected='yes'
                )
              )
            ),

            fluidRow(
              column(6, h5('Show grid?')),
              column(6,
                selectInput(ns("show_grid"), label=NULL,
                  choices=c('yes', 'no'),
                  selected='yes'
                )
              )
            ),

            fluidRow(
              column(6, h5('Show diagonal?')),
              column(6,
                selectInput(ns("dline"), label=NULL,
                  choices=c('yes', 'no'),
                  selected='yes'
                )
              )
            )
          )
        ),

        fluidRow(align='center',
          actionButton(ns('refresh'),
            label='Refresh',
            icon=icon('arrows-rotate'),
            class='btn-primary',
            style='margin-bottom: 10px;')
        )
      )
    )
  } else if (panel == 'main') {
    tagList(
      fluidRow(
        column(6, align='left',
          helpButtonUI(ns('de_scatter_help'))
        ),
        column(6, align='right',
          downloadButtonUI(ns('scatterplot_download'))
        )
      ),

      withSpinner(
        uiOutput(ns('scatterplot_out'))
      ),

      withSpinner(
        DTOutput(ns('scatter_datatable_out'))
      )
    )
  }
}


#' Scatterplot module server function
#'
#' @param id ID string used to match the ID used to call the module UI function
#' @param obj reactiveValues object containing carnation object
#' @param plot_args reactive containing 'fdr.thres' (padj threshold), 'fc.thres' (log2FC) & 'gene.to.plot' (genes to be labeled)
#'
#' @export
scatterPlotServer <- function(id, obj, plot_args){
  moduleServer(
    id,
    function(input, output, session){

      ns <- NS(id)
      config <- get_config()

      app_object <- reactive({
        list(res=obj$res)
      })

      # Reactive list of all available comparisons
      comp_all <- reactive({
        names(app_object()$res)
      })

      # Populate the x and y-axis dropdowns once data is available
      observeEvent(comp_all(), {
        validate(need(!is.null(app_object()$res), 'Waiting for data'))

        updateSelectInput(session, 'x_axis_comp',
          choices = comp_all(),
          selected = comp_all()[1]
        )
        if (length(comp_all()) == 1) {
          # if only 1 contrast, y is the same
          available_y <- comp_all()[1]
        } else {
          available_y <- comp_all()[2]
        }
        updateSelectizeInput(session, 'y_axis_comp',
          choices = comp_all(),
          selected = available_y
        )
      })

      # Swap comparisons button - does NOT trigger plot or table refresh
      observeEvent(input$swap_comp, {
        current_x <- isolate(input$x_axis_comp)
        current_y <- isolate(input$y_axis_comp)
        updateSelectInput(session, 'y_axis_comp', selected=current_x)
        updateSelectInput(session, 'x_axis_comp', selected=current_y)
      })

      # We keep axis limits in a reactiveValues
      axis_limits <- reactiveValues(lim.x=NULL, lim.y=NULL)

      # ============= Auto-scale x-axis limits =============
      observeEvent(input$scatter_x_auto, {
        autoScaleX(input, app_object(), sesion)
      })

      # ============= Auto-scale y-axis limits =============
      observeEvent(input$scatter_y_auto, {
        autoScaleY(input, app_object(), session)
      })

      # Observe 'Values to use' changes -> call both autoScaleX/Y
      observeEvent(input$compare, {
        autoScaleX(input, app_object(), session)
        autoScaleY(input, app_object(), session)
      })

      # ============= Master eventReactive: only on Refresh (or initial 'Scatter plot' main panel button click) =============
      # This merges data, sets significance, determines "show_table",
      # "plot_interactive", etc. and returns them all in a list.
      merged_data <- eventReactive({
        req(!is.null(app_object()$res))
        # The "trigger expression" – if *any* of these reactives change, we rerun:
        list(
          input$refresh,       # user clicks the "Refresh" button
          names(obj$res)       # the set of available comparisons changes
        )
        }, {
        # Merge data & set up the plot:
        validate(
          need(!is.null(obj$res), "Waiting for data"),
          need(length(names(obj$res)) > 0, "No comparisons found")
        )

        # Build minimal df
        df <- buildMergedData(input, app_object())

        # -log10 transform if P-adj
        compare <- input$compare
        if(compare == 'P-adj'){
          df[["padj.x"]] <- -log10(df[["padj.x"]])
          df[["padj.y"]] <- -log10(df[["padj.y"]])
        }

        # Validate no all-NA columns
        if(compare == 'LFC'){
          x_col <- "log2FoldChange.x"
          y_col <- "log2FoldChange.y"
        } else {
          x_col <- "padj.x"
          y_col <- "padj.y"
        }
        validate(need(!all(is.na(df[[x_col]])),
          paste(input$x_axis_comp, "column has all NAs!")))
        validate(need(!all(is.na(df[[y_col]])),
          paste(input$y_axis_comp, "column has all NAs!")))

        # 3) Determine final axis limits
        lim.x <- c(input$scatter_xmin, input$scatter_xmax)
        lim.y <- c(input$scatter_ymin, input$scatter_ymax)

        local_get_range <- function(df, column, lim){
          lim[1] <- min(df[[column]], na.rm=TRUE)
          lim[2] <- max(df[[column]], na.rm=TRUE)
          diff <- diff(lim) * 0.05
          lim[1] <- round(lim[1] - diff, 2)
          lim[2] <- round(lim[2] + diff, 2)
          lim
        }
        # Autoscale if needed
        if(is.null(lim.x[1]) || is.na(lim.x[1]) ||
           is.null(lim.x[2]) || is.na(lim.x[2])){
          lim.x <- local_get_range(df, x_col, lim.x)
        }
        if(is.null(lim.y[1]) || is.na(lim.y[1]) ||
           is.null(lim.y[2]) || is.na(lim.y[2])){
          lim.y <- local_get_range(df, y_col, lim.y)
        }
        axis_limits$lim.x <- lim.x
        axis_limits$lim.y <- lim.y

        # 4) Full data for significance labeling
        res_i_full <- as.data.frame(app_object()$res[[input$x_axis_comp]])
        res_j_full <- as.data.frame(app_object()$res[[input$y_axis_comp]])
        res_i_full$geneid <- rownames(res_i_full)
        res_j_full$geneid <- rownames(res_j_full)

        df_full <- dplyr::inner_join(
          dplyr::select(res_i_full, geneid, log2FoldChange, padj),
          dplyr::select(res_j_full, geneid, log2FoldChange, padj),
          by="geneid",
          suffix=c(".x", ".y")
        )

        # thresholds
        fdr.thres <- ifelse(is.na(plot_args()$fdr.thres), 0.1, plot_args()$fdr.thres)
        fc.thres  <- ifelse(is.na(plot_args()$fc.thres), 0, plot_args()$fc.thres)

        # labeling
        label_i <- input$x_axis_comp
        label_j <- input$y_axis_comp
        if(label_i == label_j){
          label_i <- paste0(label_i, "_x")
          label_j <- paste0(label_j, "_y")
        }

        df$significance <- dplyr::case_when(
          df_full$padj.x <= fdr.thres &
            df_full$padj.y <= fdr.thres &
            (df_full$log2FoldChange.x * df_full$log2FoldChange.y >= 0) &
            abs(df_full$log2FoldChange.x) >= fc.thres &
            abs(df_full$log2FoldChange.y) >= fc.thres ~ 'Both - same LFC sign',

          df_full$padj.x <= fdr.thres &
            df_full$padj.y <= fdr.thres &
            (df_full$log2FoldChange.x * df_full$log2FoldChange.y < 0) &
            abs(df_full$log2FoldChange.x) >= fc.thres &
            abs(df_full$log2FoldChange.y) >= fc.thres ~ 'Both - opposite LFC sign',

          (df_full$padj.x <= fdr.thres & abs(df_full$log2FoldChange.x) >= fc.thres) &
            (df_full$padj.y > fdr.thres | abs(df_full$log2FoldChange.y) < fc.thres) ~ label_i,

          (df_full$padj.y <= fdr.thres & abs(df_full$log2FoldChange.y) >= fc.thres) &
            (df_full$padj.x > fdr.thres | abs(df_full$log2FoldChange.x) < fc.thres) ~ label_j,

          TRUE ~ 'None'
        )
        df$significance <- factor(df$significance,
          levels=c('None', label_i, label_j,
                   'Both - opposite LFC sign',
                   'Both - same LFC sign')
        )

        # 5) Capture user choices for show_table and plot_interactive
        show_table_choice    <- input$show_table
        interactive_choice   <- input$plot_interactive

        list(
          df              = df,
          df_full         = df_full,
          lim.x           = lim.x,
          lim.y           = lim.y,
          show_table      = show_table_choice,
          plot_interactive= interactive_choice
        )
      })

      # ============= Build static ggplot upon Refresh =============
      scatterplot <- eventReactive({
        list(input$refresh, names(obj$res)) # User clicks Refresh or data loads
      }, {
        md <- merged_data()
        df <- md$df

        genes <- plot_args()$gene.to.plot
        if(is.null(genes) || all(genes == "")){
          lab.genes <- NULL
        } else {
          lab.genes <- genes
        }

        df <- handleOutOfBounds(df, input$compare,
                                md$lim.x, md$lim.y,
                                input$plot_all)

        color.palette <- brewer.pal(5, name=input$color.palette)

        p <- plotScatter.label(
          compare=input$compare,
          df=df,
          label_x=input$x_axis_comp,
          label_y=input$y_axis_comp,
          lab.genes=lab.genes,
          lim.x=md$lim.x,
          lim.y=md$lim.y,
          plot_all=input$plot_all,
          name.col='geneid',
          lines=c(input$vline, input$hline, input$dline),
          alpha=input$alpha,
          size=input$size,
          show.grid=input$show_grid,
          color.palette=color.palette
        )
        p
      })

      # ============= Build plotly upon Refresh =============
      scatterplot_ly <- eventReactive({
        list(input$refresh, names(obj$res)) # User clicks Refresh or data loads
      }, {
        md <- merged_data()
        df <- md$df

        genes <- plot_args()$gene.to.plot
        if(is.null(genes) || all(genes == "")){
          lab.genes <- NULL
        } else {
          lab.genes <- genes
        }

        df <- handleOutOfBounds(df, input$compare,
                                md$lim.x, md$lim.y,
                                input$plot_all)

        color.palette <- brewer.pal(5, name=input$color.palette)

        p <- plotScatter.label_ly(
          compare=input$compare,
          df=df,
          label_x=input$x_axis_comp,
          label_y=input$y_axis_comp,
          lim.x=md$lim.x,
          lim.y=md$lim.y,
          name.col='geneid',
          lines=c(input$vline, input$hline, input$dline),
          alpha=input$alpha,
          size=input$size,
          show.grid=input$show_grid,
          color.palette=color.palette,
          lab.genes=lab.genes
        )
        p
      })

      # Helper to optionally clip out-of-bounds points
      handleOutOfBounds <- function(df, compare, lim.x, lim.y, plot_all){
        df$shape <- 'in'
        if(plot_all == 'yes'){
          if(compare == 'LFC'){
            df$log2FoldChange.x <- pmin(df$log2FoldChange.x, lim.x[2])
            df$log2FoldChange.x <- pmax(df$log2FoldChange.x, lim.x[1])
            df$shape[df$log2FoldChange.x == lim.x[2]] <- 'right'
            df$shape[df$log2FoldChange.x == lim.x[1]] <- 'left'

            df$log2FoldChange.y <- pmin(df$log2FoldChange.y, lim.y[2])
            df$log2FoldChange.y <- pmax(df$log2FoldChange.y, lim.y[1])
            df$shape[df$log2FoldChange.y == lim.y[2]] <- 'above'
            df$shape[df$log2FoldChange.y == lim.y[1]] <- 'below'
          } else {
            df$padj.x <- pmin(df$padj.x, lim.x[2])
            df$padj.x <- pmax(df$padj.x, lim.x[1])
            df$shape[df$padj.x == lim.x[2]] <- 'right'
            df$shape[df$padj.x == lim.x[1]] <- 'left'

            df$padj.y <- pmin(df$padj.y, lim.y[2])
            df$padj.y <- pmax(df$padj.y, lim.y[1])
            df$shape[df$padj.y == lim.y[2]] <- 'above'
            df$shape[df$padj.y == lim.y[1]] <- 'below'
          }
          df$shape <- factor(df$shape)
        }
        df
      }

      # ============= Render the plot UI only on refresh data =============
      output$scatterplot_out <- renderUI({
        # Instead of referencing input$plot_interactive directly,
        # we check the eventReactive data
        md <- merged_data()

        # If user picks interactive, show the plotly version
        if(md$plot_interactive == 'yes'){
          p <- scatterplot_ly() %>% toWebGL()
          output$plot1 <- renderPlotly(p)
          plotlyOutput(ns('plot1'), height='600px')

        } else {
          p <- scatterplot() + theme(text=element_text(size=18))
          output$plot2 <- renderPlot(p)
          plotOutput(ns('plot2'), height='600px')
        }
      })

      # ============= Render the table only if show_table == 'yes' =============
      output$scatter_datatable_out <- renderDT({
        md <- merged_data()
        req(md$show_table == 'yes')  # only proceed if 'yes' at refresh time

        df_full <- md$df_full
        validate(need(!is.null(df_full), 'No data for table'))

        df_full <- dplyr::relocate(df_full, geneid)

        columns_to_format <- c("padj.x", "padj.y", "log2FoldChange.x", "log2FoldChange.y")
        which_cols <- which(colnames(df_full) %in% columns_to_format)
        border_cols <- c(1, grep('padj', colnames(df_full)))

        # Instead of referencing input$x_axis_comp, etc. here,
        # we can store them in md if we truly need them in the header.
        x_comp_name <- isolate(input$x_axis_comp)
        y_comp_name <- isolate(input$y_axis_comp)
        all_comps   <- c(x_comp_name, y_comp_name)

        sketch <- htmltools::withTags(
          table(
            class = 'display',
            thead(
              tr(
                th(rowspan=2, 'geneid'),
                lapply(all_comps, function(x) th(class='dt-center', colspan=2, x))
              ),
              tr(
                lapply(rep(c('log2FoldChange', 'padj'), 2), th)
              )
            )
          )
        )

        datatable(df_full,
          rownames=FALSE,
          selection='none',
          container=sketch,
          options=list(
            autoWidth=TRUE,
            columnDefs=list(
              list(className='dt-center', targets=1:(ncol(df_full)-1))
            )
          )
        ) %>%
          formatStyle(
            columns=border_cols,
            'border-right'='solid 1px'
          ) %>%
          formatSignif(columns=which_cols, digits=5)
      })

      # Help and download
      helpButtonServer('de_cmp_scatter_help', size='l')
      helpButtonServer('de_scatter_help', size='l')

      # The download always returns the static ggplot
      downloadButtonServer('scatterplot_download', scatterplot, 'scatterplot')
    }
  )
}

