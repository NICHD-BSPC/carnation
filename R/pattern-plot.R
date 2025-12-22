#' Pattern plot module
#'
#' @description
#' Module UI & server to generate pattern plots.
#'
#' @param id Module id
#' @param panel string, can be 'sidebar' or 'main'
#' @param tab string, if 'plot' show plot settings, if 'table' show table settings;
#' if 'both', show settings for both.
#' @param obj reactiveValues object containing carnation object
#' @param coldata reactiveValues object containing object metadata
#' @param plot_args reactive containing 'gene_scratchpad' (genes selected in scratchpad) &
#' 'upset_data' (list containing data from upset plot module)
#' @param config reactive list with config settings
#'
#' @returns
#' UI returns tagList with module UI
#' server invisibly returns NULL (used for side effects)
#'
#' @examples
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
#' cdata <- lapply(oobj$rld, function(x) colData(x))
#'
#' coldata <- reactiveValues( all=cdata, curr=cdata )
#'
#' plot_args <- reactive({
#'   list(
#'     gene_scratchpad=c('gene1', 'gene2'),
#'     upset_data=list(genes=NULL, labels=NULL)
#'   )
#' })
#'
#' config <- reactiveVal(get_config())
#'
#' if(interactive()){
#'   shinyApp(
#'     ui = fluidPage(
#'            sidebarPanel(
#'              patternPlotUI('p', 'sidebar', 'both'),
#'              conditionalPanel(condition = "input.pattern_mode == 'Plot'",
#'                patternPlotUI('p', 'sidebar', 'plot')
#'              ),
#'              conditionalPanel(condition = "input.pattern_mode == 'Table'",
#'                patternPlotUI('p', 'sidebar', 'table')
#'              )
#'            ),
#'            mainPanel(
#'              tabsetPanel(id='pattern_mode',
#'                tabPanel('Plot',
#'                  patternPlotUI('p', 'plot')
#'                ), # tabPanel plot
#'
#'                tabPanel('Cluster membership',
#'                  patternPlotUI('p', 'table')
#'                ) # tabPanel cluster_membership
#'
#'              ) # tabsetPanel pattern_mode
#'            ) # tabPanel pattern_analysis
#'          ),
#'     server = function(input, output, session){
#'                patternPlotServer('deg_plot', obj, coldata,
#'                                  plot_args, config)
#'              }
#'   )
#' }
#'
#' @rdname degmod
#' @name degmod
NULL

#' @rdname degmod
#' @export
patternPlotUI <- function(id, panel, tab){
  ns <- NS(id)

  # get default config
  config <- get_config()

  if(panel == 'sidebar'){
    if(tab == 'both'){
      tag <-
        tagList(
          fluidRow(
            column(6, align='left',
              tags$label(class='control-label',
                              'Analysis'
              ) # tags$label
            ), # column
            column(6, align='right',
              helpButtonUI(ns('dp_controls_help'))
            ) # column
          ), # fluidRow

          selectizeInput(ns('dp_analysis'), label=NULL,
                         choices=NULL, selected=NULL
          ) # selectizeInput
        ) # tagList
    } else if(tab == 'plot'){
      tag <-
        tagList(

          fluidRow(
            column(12, strong('Plot options')),
          ), # fluidRow

          fluidRow(
            column(6, h5('cluster column')),
            column(6,
              selectInput(ns('deg_facet'), label=NULL,
                          choices=NULL, selected=NULL
              ) # selectInput
            ) # column
          ), # fluidRow

          fluidRow(
            column(6, h5('x-axis variable')),
            column(6,
              selectInput(ns('deg_time'), label=NULL,
                          choices=NULL, selected=NULL
              ) # selectInput
            ) # column
          ), # fluidRow

          fluidRow(
            column(6, h5('group by')),
            column(6,
              selectInput(ns('deg_color'),label=NULL,
                          choices=NULL, selected=NULL
              ) # selectInput
            ) # column
          ), # fluidRow

          fluidRow(
            column(6, h5('label')),
            column(6,
              selectInput(ns('deg_label'),label=NULL,
                          choices=c('none', 'gene_scratchpad', 'upset_intersections'),
                          selected='none'
              ) # selectInput
            ) # column
          ), # fluidRow

          conditionalPanel(paste0('input["', ns('deg_label'), '"] == "upset_intersections"'),
            fluidRow(
              column(12,
                selectizeInput(ns('upset_intersect'), label=NULL,
                               choices=NULL, selected=NULL)
              ) # column
            ) # fluidRow
          ), # conditionalPanel

          bsCollapse(
            bsCollapsePanel('cluster settings',

              fluidRow(
                column(4, 'clusters to show'),
                column(8,
                  selectizeInput(ns('facet_var_levels'),
                                 label=NULL,
                                 choices=NULL,
                                 selected=NULL,
                                 multiple=TRUE
                  ) # selectizeInput
                ) # column
              ), # fluidRow

              fluidRow(
                column(6, align='center',
                       style='margin-bottom: 10px;',
                  actionButton(ns('deg_facet_all'), 'Select all')
                ), # column
                column(6, align='center',
                       style='margin-bottom: 10px;',
                  actionButton(ns('deg_facet_none'), 'Select none')
                ) # column
              ), # fluidRow

              fluidRow(
                column(6, h5('min cluster size')),
                column(6,
                  numericInput(ns('deg_minc'), label=NULL,
                    value=config$ui$pattern_analysis$min_cluster_size
                  ) # numericInput
                ) # column
              ) # fluidRow

            ), # bsCollapsePanel

            # x axis settings
            bsCollapsePanel('x-axis settings',

              fluidRow(
                column(12, align='left',
                       'variable levels'),
              ), # fluidRow

              uiOutput(ns('x_axis_bucket')),

              fluidRow(
                column(6, align='center',
                       style='margin-bottom: 10px;',
                  actionButton(ns('deg_x_all'), 'Select all')
                ), # column
                column(6, align='center',
                       style='margin-bottom: 10px;',
                  actionButton(ns('deg_x_none'), 'Select none')
                ) # column
              ), # fluidRow

              fluidRow(
                column(4, h5('rotate x labels')),
                column(8,
                  numericInput(ns('deg_rotate'), label=NULL,
                    value=config$ui$pattern_analysis$x_rotate,
                    step=15
                  ) # numericInput
                ) # column
              ) # fluidRow

            ), # bsCollapsePanel

            bsCollapsePanel('More options',

              fluidRow(
                column(6, h5('trendline')),
                column(6,
                  selectInput(ns('deg_smooth'), label=NULL,
                              choices=c('line', 'smooth', 'none')
                  ) # selectInput
                ) # column
              ), # fluidRow

              fluidRow(
                column(6, h5('show boxes?')),
                column(6,
                  selectInput(ns('deg_boxes'), label=NULL,
                              choices=c(TRUE, FALSE),
                              selected=TRUE
                  ) # selectInput
                ) # column
              ), # fluidRow

              fluidRow(
                column(6, h5('show points?')),
                column(6,
                  selectInput(ns('deg_points'), label=NULL,
                              choices=c(TRUE, FALSE),
                              selected=FALSE
                  ) # selectInput
                ) # column
              ), # fluidRow

              fluidRow(
                column(6, h5('show lines?')),
                column(6,
                  selectInput(ns('deg_lines'), label=NULL,
                              choices=c(TRUE, FALSE),
                              selected=TRUE
                  ) # selectInput
                ) # column
              ), # fluidRow

              fluidRow(
                column(6, h5('text scale')),
                column(6,
                  numericInput(ns('txt_scale'), label=NULL,
                               value=1, min=0, max=3, step=0.1
                  ) # numericInput
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
    } else if(tab == 'table'){
      tag <-
        tagList(

          fluidRow(
            column(6, 'Cluster column'),
            column(6,
              selectInput(ns('deg_cluster'), label=NULL,
                          choices=NULL,
                          selected=NULL
              ) # selectInput
            ) # column
          ), # fluidRow
          wellPanel(
            selectizeInput(ns('deg_cluster_levels'),
                           label='Clusters to show',
                           choices=NULL, selected=NULL,
                           multiple=TRUE
            ), # selectInput

            fluidRow(
              column(6, align='center',
                     style='margin-bottom: 20px;',
                actionButton(ns('deg_cluster_all'), 'Select all')
              ), # column
              column(6, align='center',
                     style='margin-bottom: 20px;',
                actionButton(ns('deg_cluster_none'), 'Select none')
              ) # column
            ) # fluidRow
          ) # wellPanel
        )  # tagList
    }
  } else if(panel == 'plot'){
    tag <-
      tagList(
        fluidRow(
          column(6, align='left',
            helpButtonUI(ns('pattern_plt_help'))
          ), # column
          column(6, align='right',
            downloadButtonUI(ns('dp_download'))
          ) # column
        ), # fluidRow

        withSpinner(
          plotlyOutput(ns('dp_plot'), height='700px')
        ), # withSpinner
        uiOutput(ns('label_tbl'))
      ) # tagList
  } else if(panel == 'table'){
    tag <-
      tagList(
        fluidRow(
           column(12, align='left',
            helpButtonUI(ns('pattern_tbl_help'))
          ) # column
        ), # fluidRow
        withSpinner(
          DTOutput(ns('dp_table'))
        ) # withSpinner
      ) # tagList
  }
  tag
}


#' @rdname degmod
#' @export
patternPlotServer <- function(id,
                              obj,
                              coldata,
                              plot_args,
                              config){
  moduleServer(
    id,

    function(input, output, session){

      ns <- NS(id)

      res_obj <- reactive({
        obj$res
      })

      pattern_obj <- eventReactive(res_obj(), {
        obj$degpatterns
      }, ignoreNULL=FALSE)

      pattern_coldata <- reactive({
        coldata$curr
      })

      # reactive to keep track of DEG plot data
      deg_plot_data <- reactiveValues(obj=NULL)

      xchoices <- reactiveValues(all=NULL, current=NULL)

      current_facet_levels <- reactiveValues(l=NULL)
      upset_data <- reactiveValues(genes=NULL, labels=NULL)
      gene_scratchpad <- reactiveValues(genes=NULL)

      # update from reactive config
      observeEvent(config(), {
        updateNumericInput(session, 'min_cluster_size',
                           value=config()$ui$pattern_analysis$min_cluster_size)
        updateNumericInput(session, 'x_rotate',
                           value=config()$ui$pattern_analysis$x_rotate)
      })

      observeEvent(pattern_obj(), {
        obj <- pattern_obj()

        if(is.null(names(obj))){
          updateSelectizeInput(session,
                               'dp_analysis',
                               choices=c(''),
                               selected=c(''))
        } else {
          updateSelectizeInput(session,
                               'dp_analysis',
                               choices=names(obj),
                               selected=names(obj)[1])
        }

        # reset reactive values & menus on initial load
        deg_plot_data$obj <- NULL
        xchoices$all <- NULL
        xchoices$current <- NULL
        current_facet_levels$l <- NULL

        updateSelectizeInput(session, 'facet_var_levels',
                             choices='',
                             selected='')

        updateSelectInput(session,
                          'deg_facet',
                          choices='',
                          selected='')

        updateSelectInput(session,
                          'deg_cluster',
                          choices='',
                          selected='')

        updateSelectizeInput(session,
                             'deg_cluster_levels',
                             choices='',
                             selected='')

        updateSelectInput(session,
                          'deg_time',
                          choices='',
                          selected='')

        updateSelectInput(session,
                          'deg_color',
                          choices='',
                          selected='')
      }, ignoreNULL=FALSE)

      # update upset intersections menu
      observeEvent(plot_args()$upset_data, {
        upset_genes <- plot_args()$upset_data$genes
        upset_labels <- plot_args()$upset_data$labels

        # only update if changed
        updt <- FALSE
        if(is.null(upset_data$labels)){
          updt <- TRUE
        } else if(length(unlist(upset_genes)) != length(unlist(upset_data$genes))){
          updt <- TRUE
        } else if(sum(unlist(upset_genes) != unlist(upset_data$genes)) > 0){
          updt <- TRUE
        }

        if(updt){
          upset_data$genes <- upset_genes
          upset_data$labels <- upset_labels

          updateSelectizeInput(session, 'upset_intersect',
                               choices=upset_data$labels,
                               server=TRUE)
        }
      })

      observeEvent(plot_args()$gene.to.plot, {
        gene_scratchpad$genes <- plot_args()$gene.to.plot
      })

      # observer to get deg plot data and update menus
      observeEvent(c(input$dp_analysis, pattern_coldata()), {
        validate(
          need(input$dp_analysis != '' & input$dp_analysis %in% names(pattern_obj()), 'Waiting for selection')
        )

        validate(
          need(!is.null(pattern_coldata()), 'Waiting for selection')
        )

        obj <- pattern_obj()[[input$dp_analysis]]
        if(!is.data.frame(obj)) obj <- obj$normalized

        # get metadata and attach if necessary
        # NOTE: we get coldata from 'all_samples', which has
        # to be present
        cdata <- pattern_coldata()[['all_samples']]
        if(!all(colnames(cdata) %in% colnames(obj))){
            obj <- add_metadata(obj,
                                cdata,
                                config()$server$cols.to.drop)
        }

        # get time variable
        cols <- colnames(obj)

        # remove common metadata columns if available
        cols <- setdiff(cols, config()$server$cols.to.drop)

        # remove more columns from color options
        cols <- setdiff(cols,
                        config()$server$pattern_analysis$cols.to.drop)

        # cluster options are 'cluster' and 'cutoff*'
        cluster.options <- c('cluster', cols[grep('cutoff', cols)])

        # remove cluster columns from time variable
        # and only keep columns found in main metadata
        col.keep <- setdiff(cols, cluster.options)
        if(is.null(cdata)) time <- NULL
        else time <- intersect(col.keep, colnames(cdata))

        # if intersection is null, keep the original columns
        if(is.null(time)) time <- setdiff(cols, cluster.options)

        # remove 'symbol' from time options
        time <- setdiff(time, 'symbol')

        # set 'group' as default time variable if available
        if(length(time) > 1 & 'group' %in% time){
          sel_time <- 'group'
        } else {
          sel_time <- time[1]
        }

        if(!is.null(input$deg_facet) & input$deg_facet %in% cluster.options){
          cluster_selected <- input$deg_facet
        } else {
          cluster_selected <- cluster.options[1]
        }

        updateSelectInput(session,
                        'deg_facet',
                        choices=cluster.options,
                        selected=cluster_selected)

        updateSelectInput(session,
                        'deg_cluster',
                        choices=cluster.options,
                        selected=cluster_selected)

        if(!is.null(input$deg_time) & input$deg_time %in% time){
            updateSelectInput(session,
                              'deg_time',
                              choices=time,
                              selected=input$deg_time)
        } else {
            updateSelectInput(session,
                              'deg_time',
                              choices=time,
                              selected=sel_time)
        }

        all.colors <- c('none', time)
        if(!is.null(input$deg_color) & input$deg_color %in% all.colors){
            updateSelectInput(session,
                              'deg_color',
                              choices=all.colors,
                              selected=input$deg_color)
        } else {
            updateSelectInput(session,
                              'deg_color',
                              choices=all.colors,
                              selected=all.colors[1])
        }

        # save degplot data
        deg_plot_data$obj <- obj

      }, ignoreNULL=FALSE) # observeEvent get degplot data

      observeEvent(c(input$deg_cluster, deg_plot_data$obj,
                     input$deg_cluster_all), {

        validate(
          need(!is.null(deg_plot_data$obj), 'Waiting for selection')
        )
        obj <- deg_plot_data$obj
        validate(
          need(input$deg_cluster %in% colnames(obj),
               'Waiting for selection')
        )

        cluster_levels <- unique(obj[[input$deg_cluster]])
        updateSelectizeInput(session, 'deg_cluster_levels',
                             choices=cluster_levels,
                             selected=cluster_levels)
      })

      observeEvent(input$deg_cluster_none, {
        updateSelectizeInput(session, 'deg_cluster_levels',
                             selected='')
      })

      ####################### facet settings ###################

      # observer to update facet var levels
      # when facet var, coldata updated or facet var reset
      observeEvent(c(input$deg_facet,
                     deg_plot_data$obj,
                     input$deg_minc,
                     input$deg_facet_all), {
        validate(
          need(!is.null(deg_plot_data$obj), 'Waiting for selection')
        )
        obj <- deg_plot_data$obj
        validate(
          need(input$deg_facet %in% colnames(obj), 'Waiting for selection')
        )

        # check min cluster size

        # get cluster sizes
        nodup.idx <- !duplicated(obj$genes)
        cluster.sizes <- table(obj[nodup.idx, input$deg_facet])

        # if input is NA, set to 1
        if(!is.na(input$deg_minc) & all(cluster.sizes < input$deg_minc)){
          facet_levels <- NULL
        } else {
          facet_levels <- names(cluster.sizes)[cluster.sizes > input$deg_minc]
        }

        updateSelectizeInput(session, 'facet_var_levels',
                             choices=facet_levels,
                             selected=facet_levels)
      }, ignoreNULL=FALSE)

      observeEvent(input$deg_facet_none, {
        updateSelectizeInput(session, 'facet_var_levels',
                             selected='')
      })

      ####################### x axis settings ###################

      # observer to save initial factor levels of x axis variable
      observeEvent(c(input$deg_time, deg_plot_data$obj), {
        validate(
          need(!is.null(input$deg_time) & input$deg_time != '',
               'Waiting for selection')
        )

        validate(
          need(!is.null(pattern_coldata()), 'Waiting for selection')
        )

        obj <- deg_plot_data$obj
        validate(
          need(input$deg_time %in% colnames(obj), 'Waiting for selection')
        )

        # get initial factor levels
        cdata <- obj[,input$deg_time]

        # if factor, use existing levels
        # if numeric, sort ascending
        # else, use unique values
        if(is.factor(cdata)){
          xchoices$all <- levels(cdata)
        } else if(is.numeric(cdata)){
          lvls <- unique(cdata)
          xchoices$all <- lvls[order(lvls)]
        } else {
          xchoices$all <- unique(cdata)
        }

      })

      # observer to update factor levels when 'Select all'
      # button is clicked
      observeEvent(c(xchoices$all, input$deg_x_all), {
        xchoices$current <- xchoices$all
      })

      observeEvent(input$deg_x_order, {
        # update xchoices if different from deg_x_order
        if(any(xchoices$current != input$deg_x_order))
          xchoices$current <- input$deg_x_order
      })

      # observer to update factor levels when 'Select none'
      # button is clicked
      observeEvent(input$deg_x_none, {
        xchoices$current <- NULL
      }, ignoreNULL=FALSE)

      # x-axis bucket list UI
      output$x_axis_bucket <- renderUI({
        tagList(
          fluidRow(
            column(12,
              bucket_list(
                header=NULL,
                group_name='xlabels_group',
                class=c('default-sortable','custom-sortable'),
                orientation='horizontal',
                add_rank_list(
                    text = 'current',
                    labels = xchoices$current,
                    input_id=ns('deg_x_order')
                ), # add_rank_list
                add_rank_list(
                    text = 'unused',
                    labels = setdiff(xchoices$all,
                                     xchoices$current),
                    input_id=ns('deg_x_other')
                ) # add_rank_list
              ) # bucket_list
            ) # column
          ) # fluidRow
        ) # tagList
      })

      observeEvent(input$dp_analysis, {
        current_facet_levels$l <- NULL
        removeNotification('clust_filt_notify')
      })
      ############################################################

      # reactive to make degPatterns plot
      degplot <- eventReactive(c(deg_plot_data$obj, input$plot_do), {

        validate(
            need(!is.null(deg_plot_data$obj),
                 'Pattern analysis not available')
        )

        validate(
            need(input$deg_time %in% colnames(deg_plot_data$obj),
            'Loading')
        )

        validate(
          need(length(input$facet_var_levels) > 0,
               paste0('Must have at least 1 cluster to plot. ',
                      'Please adjust "cluster settings"'))
        )

        obj <- deg_plot_data$obj

        # get time variable
        time <- input$deg_time

        # get color variable
        if(input$deg_color == 'none') color <- NULL
        else color <- input$deg_color

        ## check min cluster size

        ## get cluster sizes
        nodup.idx <- !duplicated(obj$genes)
        cluster.sizes <- table(obj[nodup.idx, input$deg_facet])

        # if input is NA, set to 1
        if(!is.na(input$deg_minc) & all(cluster.sizes < input$deg_minc)){
            validate(
                need(any(cluster.sizes > input$deg_minc),
                     paste0('No clusters left after filtering. ',
                            'Please choose lower min cluster size'))
            )

        }

        # number of clusters removed
        nrem <- length(cluster.sizes) - length(input$facet_var_levels)
        # check to see if any changes made to current facet levels
        if(is.null(current_facet_levels$l)){
          nchange <- length(input$facet_var_levels)
        } else {
          nchange <- length(setdiff(current_facet_levels$l, input$facet_var_levels))
        }

        if(nrem > 0 & nchange > 0){
            if(nrem == 1) suff <- ''
            else suff <- 's'

            showNotification(
              id='clust_filt_notify',
              paste0('Note: ', nrem, ' cluster', suff,
                     ' filtered out at current settings. To view more, try adjusting "min cluster size" or "facet settings"'),
              duration=10
            )
        }

        # get current x-axis order
        x_order <- xchoices$current
        validate(
          need(length(x_order) > 0, 'Must have at least one value to plot on x-axis')
        )

        current_facet_levels$l <- input$facet_var_levels

        if(input$deg_label == 'gene_scratchpad') g <- plot_args()$gene_scratchpad
        else if(input$deg_label == 'upset_intersections')
          g <- upset_data$genes[[ input$upset_intersect ]]
        else g <- NULL

        if(length(g) > 10){
          showNotification(
            'Too many labeled genes. Using a single color to label',
            type='warning'
          )
        }

        p <- get_degplot(obj, time=time, color=color,
                         cluster_column=input$deg_facet,
                         cluster_to_show=input$facet_var_levels,
                         x_order=x_order,
                         points=as.logical(input$deg_points),
                         lines=as.logical(input$deg_lines),
                         boxes=as.logical(input$deg_boxes),
                         smooth=input$deg_smooth,
                         genes_to_label=g)

        if(input$txt_scale == 0 | is.na(input$txt_scale)){
          txt_scale <- 1
          showNotification(
            'Invalid text scaling factor. Resetting to default',
            type='warning'
          )
        } else {
          txt_scale <- input$txt_scale
        }

        p + ggtitle(paste0('Analysis: ', input$dp_analysis)) +
            theme(axis.text.x=element_text(size=9*txt_scale,
                                           angle=input$deg_rotate),
                  axis.text.y=element_text(size=9*txt_scale),
                  strip.text.x=element_text(size=9*txt_scale),
                  legend.title=element_blank(),
                  legend.text=element_text(size=9*txt_scale))
      }, ignoreNULL=FALSE) # eventReactive degplot

      output$dp_plot <- renderPlotly({
        degplot()
      })

      get_label_tbl <- eventReactive(c(deg_plot_data$obj,
                                       input$plot_do), {

        validate(
          need(!is.null(deg_plot_data$obj), '')
        )

        table <- deg_plot_data$obj

        validate(
          need(!is.null(input$deg_facet) & input$deg_facet %in% colnames(table), '')
        )
        cols.to.keep <- c('genes', 'symbol', input$deg_facet)
        table <- unique(table[, cols.to.keep])

        # filter by min cluster size
        cluster.sizes <- table(table[, input$deg_facet])
        clusters.to.keep <- names(cluster.sizes)[cluster.sizes > input$deg_minc]

        # only keep clusters chosen in facet_var_levels
        clusters.to.keep <- intersect(clusters.to.keep, input$facet_var_levels)

        cidx <- table[, input$deg_facet] %in% clusters.to.keep
        table <- table[cidx,]

        if(input$deg_label == 'gene_scratchpad'){
          g <- plot_args()$gene_scratchpad
        } else if(input$deg_label == 'upset_intersections'){
          g <- plot_args()$upset_data$genes[[ input$upset_intersect ]]
        } else {
          g <- NULL
        }

        df2 <- NULL
        if(!is.null(g)){
          if('symbol' %in% colnames(table)){
              gene_col <- 'symbol'
          } else {
              gene_col <- 'genes'
          }

          gidx <- table[[ gene_col ]] %in% g
          df <- table[gidx,]
          df <- df[order(df[,input$deg_facet]),]

          all_facet <- table(df[, input$deg_facet])
          facet_genes <- unlist(lapply(names(all_facet), function(x){
                           idx <- df[, input$deg_facet] == x
                           paste(df[idx, gene_col], collapse=',')
                          }))

          df2 <- data.frame(cluster=names(all_facet),
                            count=as.numeric(all_facet),
                            genes=format_genes(facet_genes, genes.per.line=15, sep=','))
          colnames(df2)[1] <- input$deg_facet
        }

        list(tbl=df2, ngenes=nrow(df))
      })

      output$labtbl <- renderDT({
        ll <- get_label_tbl()
        ll$tbl %>%
            datatable(rownames=FALSE,
                      selection='none',
                      caption=tags$caption(style='font-weight: bold; font-size: 15px;',
                                           paste0('Labeled genes (n = ', ll$ngenes, ')')),
                      options=list(dom='tp'))
      })

      output$label_tbl <- renderUI({
        fluidRow(
          column(10, align='center',
            DTOutput(ns('labtbl'))
          )
        )
      })

      # reactive to get degpatterns tabls
      degtable <- eventReactive(c(deg_plot_data$obj,
                                  input$deg_cluster,
                                  input$deg_cluster_levels), {

        validate(
          need(!is.null(deg_plot_data$obj),
               'Pattern analysis not available')
        )

        obj <- deg_plot_data$obj
        df <- obj

        cols.to.keep <- c('genes', 'symbol', input$deg_cluster)
        df <- df[ ,cols.to.keep]
        df <- unique(df)

        # only show selected clusters
        cluster_col <- input$deg_cluster
        cluster_levels <- input$deg_cluster_levels
        validate(
          need(length(cluster_levels) > 0, 'Must have at least one cluster to show')
        )
        df <- df[as.character(df[[cluster_col]]) %in% cluster_levels,]

        df
      }, ignoreNULL=FALSE) # eventReactive degtable

      output$dp_table <- renderDT({
        degtable()
      }, rownames=FALSE)

      # pattern analysis help
      helpButtonServer('dp_controls_help', size='l')
      helpButtonServer('pattern_plt_help', size='l')
      helpButtonServer('pattern_tbl_help', size='l')

      downloadButtonServer('dp_download', degplot,
                           paste0(input$dp_analysis, '-degplot'))

    } # function
  ) # moduleServer
}
