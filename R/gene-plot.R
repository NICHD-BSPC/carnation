#' Gene plot module
#'
#' @description
#' UI & server for module to create gene plot
#'
#' @param id Module id
#' @param panel string, can be 'sidebar' or 'main'
#' @param obj reactiveValues object containing carnation object
#' @param coldata reactiveValues object containing object metadata
#' @param plot_args reactive list with 3 elements: 'gene.id' (all gene IDs) & 'gene_scratchpad'
#' (genes selected in scratchpad) & 'comp_all' (selected comparison)
#' @param config reactive list with config settings
#'
#'
#' @examples
#' library(DESeq2)
#' library(shiny)
#'
#' dds <- makeExampleDESeqDataSet()
#' rld <- varianceStabilizingTransformation(dds, blind=TRUE)
#'
#' dds <- DESeq(dds)
#' results <- results(dds, contrast=c('condition', 'A', 'B'))
#'
#' # Create reactive values to simulate app state
#' obj <- reactiveValues(
#'   dds = list(main = dds),
#'   rld = list(main = rld),
#'   res = list(comp1 = results),
#'   all_dds = dds,
#'   all_rld = rld,
#'   dds_mapping = list(comp1 = 'main')
#' )
#'
#' # Set up coldata structure that the module expects
#' coldata <- reactiveValues(
#'   curr = list(
#'     all_samples = colData(dds),
#'     main = colData(dds)
#'   )
#' )
#'
#' plot_args <- reactive({
#'   list(
#'     gene.to.plot = c("gene1", "gene2"),
#'     gene.id = rownames(dds),
#'     comp_all = "comp1"
#'   )
#' })
#'
#' config <- reactiveVal(get_config())
#'
#' if(interactive()){
#'   shinyApp(
#'     ui = fluidPage(
#'            sidebarPanel(genePlotUI('p', 'sidebar')),
#'            mainPanel(genePlotUI('p', 'main'))
#'          ),
#'     server = function(input, output, session){
#'                genePlotServer('p', obj, coldata, plot_args, config)
#'              }
#'   )
#' }
#'
#' @name geneplotmod
#' @rdname geneplotmod
NULL


#' @rdname geneplotmod
#' @export
genePlotUI <- function(id, panel){
  ns <- NS(id)

  # load default config
  config <- get_config()
  if(panel == 'sidebar'){
    tag <-
      tagList(

        fluidRow(
          column(6, strong('Plot options')),
          column(6, align='right',
            helpButtonUI(ns('geneplt_opts_help'))
          ) # column
        ), # fluidRow

        fluidRow(
          column(4, h5('sample group')),
          column(8,
            selectInput(ns('samples'), label=NULL,
                        choices=NULL, selected=NULL
            ) # selectInput
          ) # column
        ), # fluidRow

        fluidRow(
          column(4, h5('x-axis variable')),
          column(8,
            selectInput(ns('xvar'), label=NULL,
                        choices=NULL
            ) # selectInput
          ) # column
        ), # fluidRow

        fluidRow(
          column(4, h5('facet by')),
          column(8,
            selectizeInput(ns('facet'), label=NULL,
                           choices=NULL, selected=NULL,
                           multiple=TRUE
            ) # selectizeInput
          ) # column
        ), # fluidRow

        fluidRow(
          column(4, h5('color by')),
          column(8,
            selectizeInput(ns('color'), label=NULL,
                           choices=NULL,
                           selected=NULL
            ) # selectizeInput
          ) # column
        ), # fluidRow

        bsCollapse(
          # reactive ui for x axis bucket list
          bsCollapsePanel('x-axis settings',

            fluidRow(
              column(12, align='left',
                     'variable levels')
            ), # fluidRow

            uiOutput(ns('x_axis_bucket')),

            fluidRow(
              column(6, align='center',
                     style='margin-bottom: 10px;',
                actionButton(ns('x_all'), 'Select all')
              ), # column
              column(6, align='center',
                     style='margin-bottom: 10px;',
                actionButton(ns('x_none'), 'Select none')
              ) # column
            ), # fluidRow

            fluidRow(
              column(4, h5('rotate x labels')),
              column(8,
                numericInput(ns('x_rotate'), label=NULL,
                  value=config$ui$de_analysis$gene_plot$x_rotate,
                  step=15)
              ) # column
            ) # fluidRow

          ), # bsCollapsePanel

          bsCollapsePanel('y-axis settings',

            fluidRow(
              column(4, h5('y-axis max')),
              column(8,
                numericInput(ns('ymax'), label=NULL, value=NA)
              ) # column
            ), # fluidRow

            fluidRow(
              column(4, h5('y-axis min')),
              column(8,
                numericInput(ns('ymin'), label=NULL, value=NA)
              ) # column
            ), # fluidRow

            fluidRow(
              column(4, align='left', style='margin-bottom: 10px;',
                actionButton(ns('y_auto'), label='Autoscale')
              ) # column
            ), # fluidRow

            fluidRow(
              column(4, h5('log scale')),
              column(8, align='left',
                checkboxInput(ns('logy'), label=NULL,
                  value=config$ui$de_analysis$gene_plot$logy
                ) # checkboxInput
              ) # column
            ), # fluidRow

            fluidRow(
              column(4, h5('free y axes')),
              column(8, align='left',
                checkboxInput(ns('freey'), label=NULL,
                  value=config$ui$de_analysis$gene_plot$freey
                ) # checkboxInput
              ) # column
            ) # fluidRow
          ), # bsCollapsePanel

          bsCollapsePanel('facet settings',

            fluidRow(
              column(4, 'variable'),
              column(8,
                selectInput(ns('facet_vars'),
                            label=NULL,
                            choices=NULL,
                            selected=NULL
                ) # selectInput
              ) # column
            ), # fluidRow
            fluidRow(
              column(4, 'levels'),
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
              column(12, align='right',
                actionButton(ns('reset_facet_vars'),'Reset',
                             style='margin-bottom: 20px')
              ) # column
            ) # fluidRow

          ), # bsCollapsePanel

          bsCollapsePanel('More options',

            fluidRow(
              column(4, h5('normalization')),
              column(8,
                selectInput(ns('norm_method'), label=NULL,
                  choices=unlist(config$ui$de_analysis$gene_plot$norm_method)
                ) # selectInput
              ) # column
            ), # fluidRow

            fluidRow(
              column(4, h5('trendline')),
              column(8, align='left',
                selectInput(ns('trendline'), label=NULL,
                  choices=c('line','smooth')
                ) # selectInput
              ) # column
            ), # fluidRow

            fluidRow(
              column(4, h5('boxes')),
              column(8, align='left',
                checkboxInput(ns('boxes'), label=NULL,
                  value=config$ui$de_analysis$gene_plot$boxes
                ) # checkboxInput
              ) # column
            ), # fluidRow

            fluidRow(
              column(4, h5('legend')),
              column(8, align='left',
                checkboxInput(ns('legend'), label=NULL,
                  value=config$ui$de_analysis$gene_plot$legend
                ) # checkboxInput
              ) # column
            ), # fluidRow

            fluidRow(
              column(4, h5('text scale')),
              column(8,
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
  } else if(panel == 'main'){
    tag <- tagList(
             fluidRow(
                column(6, align='left',
                 helpButtonUI(ns('de_geneplt_help'))
               ), # column
               column(6, align='right',
                 downloadButtonUI(ns('normplot_download'))
               ) # column
             ), # fluidRow

             withSpinner(
               plotlyOutput(ns('normplot'))
             ) # withSpinner
           )
  }
  tag
}

#' @rdname geneplotmod
#' @export
genePlotServer <- function(id, obj,
                           coldata,
                           plot_args,
                           config){
  moduleServer(
    id,

    function(input, output, session){

      ns <- NS(id)

      xchoices <- reactiveValues(all=NULL, current=NULL)

      gene_plot_data <- reactiveValues(all=NULL, plotted=NULL, handle=NULL)

      gene_coldata <- eventReactive(c(coldata$curr, input$samples), {
        validate(
          need(!is.null(coldata$curr) & input$samples != '', 'Waiting for selection')
        )

        coldata$curr[[ input$samples ]]
      })

      app_object <- reactive({
          list(res=obj$res,
               dds=obj$dds,
               rld=obj$rld,
               dds_mapping=obj$dds_mapping,
               all_dds=obj$all_dds,
               all_rld=obj$all_rld)
      })

      # faceting columns that are always shown
      facet.extra <- c('gene')

      cols.to.drop <- reactive({ config()$server$cols.to.drop })

      # update menus from reactive config
      observeEvent(config(), {

        updateSelectInput(session, 'norm_method',
                          choices=unlist(config()$ui$de_analysis$gene_plot$norm_method))
        updateNumericInput(session, 'x_rotate',
                           value=config()$ui$de_analysis$gene_plot$x_rotate)

        # update checkbox inputs
        for(name in c('logy', 'freey', 'boxes', 'legend')){
          updateCheckboxInput(session, name,
                              value=config()$ui$de_analysis$gene_plot[[ name ]])
        }
      })

      # observer to update samples menu for gene plot
      observeEvent(c(app_object()$res, plot_args()$comp_all), {
        validate(
          need(!is.null(app_object()$res) & !is.null(app_object()$dds_mapping) & !is.null(plot_args()$comp_all), 'Waiting for selection')
        )

        validate(
          need(plot_args()$comp_all %in% names(app_object()$dds_mapping), 'Waiting for selection')
        )

        if(!is.null(app_object()$all_dds)){
            if(length(app_object()$dds) == 1){
                choices <- 'all_samples'
            } else {
                choices <- c(names(app_object()$dds),
                             'all_samples')
            }
            if(!is.null(input$samples) & input$samples %in% choices) selected <- input$samples
            else selected <- choices[1]
        } else {
            choices <- names(app_object()$dds)
            selected <- choices[1]
        }
        updateSelectInput(session,
                          'samples',
                          choices=choices,
                          selected=selected)
      }) # observeEvent update samples menu

      # observer to update gene plot data and y-axis limits
      # NOTE: keeping this here to avoid passing all these inputs
      # to the module
      observeEvent(c(app_object(),
                     plot_args(),
                     input$samples,
                     input$norm_method), {

        validate(
          need(!is.null(app_object()$res) & !is.null(app_object()$dds_mapping) & !is.null(plot_args()$comp_all), 'Waiting for selection')
        )

        validate(
          need(all(names(app_object()$dds_mapping) %in% names(app_object()$res)), 'Waiting for selection')
        )

        sample_grp <- input$samples
        if(sample_grp != 'all_samples'){
            if(input$norm_method == 'vst') rld.i <- app_object()$rld[[ sample_grp ]]
            else if(input$norm_method == 'libsize') rld.i <- app_object()$dds[[ sample_grp ]]
        } else if(sample_grp == 'all_samples'){
            if(input$norm_method == 'vst') rld.i <- app_object()$all_rld
            else if(input$norm_method == 'libsize') rld.i <- app_object()$all_dds
        }

        gene_plot_data$all <- rld.i
      })

      # update gene plot controls when col.data changes
      observeEvent(gene_coldata(), {

        validate(
          need(!is.null(gene_coldata()), 'Waiting for selection')
        )

        # get relevant columns from metadata
        int.cols <- colnames(gene_coldata())[!colnames(gene_coldata()) %in% cols.to.drop()]

        if(!input$xvar %in% int.cols){
            selected <- ifelse('group' %in% int.cols, 'group', int.cols[1])
            updateSelectizeInput(session, 'xvar',
                                 choices=int.cols,
                                 selected=selected)
        } else {
            updateSelectizeInput(session, 'xvar',
                                 choices=int.cols,
                                 selected=input$xvar)
        }

        # update color choices
        if(is.null(input$color) || !(input$color %in% c(int.cols, 'sample', 'gene'))){
            updateSelectizeInput(session, 'color',
                                 choices=c(int.cols, 'sample','gene'),
                                 selected='gene')
        } else {
            updateSelectizeInput(session, 'color',
                                 choices=c(int.cols, 'sample','gene'),
                                 selected=input$color)
        }

        # any columns not being used to group can be used to facet
        # NOTE: this section is only used once per new comparison selected
        facet.cols <- colnames(gene_coldata())[!colnames(gene_coldata()) %in% c(input$xvar, cols.to.drop())]
        facet.cols <- c(facet.extra, facet.cols)
        if(length(facet.cols) > 0){
            if(is.null(input$facet) || !(input$facet %in% facet.cols)){
                updateSelectizeInput(session, 'facet',
                                     choices=facet.cols,
                                     selected=''
                                     )
            } else {
                updateSelectizeInput(session, 'facet',
                                     choices=facet.cols,
                                     selected=input$facet
                                     )
                facet_vars <- setdiff(input$facet, c('', 'none', 'gene'))
                if(length(facet_vars) > 0){
                    if(input$facet_vars %in% facet_vars) selected <- input$facet_vars
                    else selected <- facet_vars[1]
                    updateSelectInput(session, 'facet_vars',
                                      choices=facet_vars,
                                      selected=selected)
                }
            }
        }

        # reset ymax/ymin for gene plot
        updateNumericInput(session, 'ymax', value=NA)
        updateNumericInput(session, 'ymin', value=NA)
      }) # observeEvent update coldata()$gene

      # update faceting variable when new grouping variable is selected
      observeEvent(c(input$xvar, gene_coldata(), app_object()), {
        validate(
          need(!is.null(input$xvar) & input$xvar != '',
               'Waiting for selection')
        )
        # any columns not being used to group can be used to facet
        facet.cols <- colnames(gene_coldata())[!colnames(gene_coldata()) %in% c(input$xvar, cols.to.drop())]
        facet.cols <- c(facet.extra, facet.cols)
        if(length(facet.cols) > 0){
            if(is.null(input$facet) || !(input$facet %in% facet.cols)){
                updateSelectizeInput(session, 'facet',
                                     choices=facet.cols,
                                     selected='')

            } else {
                updateSelectizeInput(session, 'facet',
                                     choices=facet.cols,
                                     selected=input$facet)
                facet_vars <- setdiff(input$facet, c('', 'none', 'gene'))
                if(length(facet_vars) > 0){
                    if(input$facet_vars %in% facet_vars) selected <- input$facet_vars
                    else selected <- facet_vars[1]
                    updateSelectInput(session, 'facet_vars',
                                      choices=facet_vars)
                }
            }
        }
      }, ignoreNULL=FALSE) # observeEvent update faceting var

      # observer to update facet var levels
      # when facet var, coldata updated or facet var reset
      observeEvent(c(input$facet, input$facet_vars, gene_coldata(), input$reset_facet_vars), {
        validate(
          need(!is.null(gene_coldata()), 'Waiting for selection')
        )
        cdata <- gene_coldata()
        validate(
          need(input$facet_vars %in% colnames(cdata), 'Waiting for selection')
        )

        levels <- unique(cdata[,input$facet_vars])
        updateSelectizeInput(session, 'facet_var_levels',
                             choices=levels,
                             selected=levels)
      }, ignoreNULL=FALSE)

      # observer to update facet vars & facet var levels
      # when facet or coldata changes
      observeEvent(c(input$facet, gene_coldata()), {
        validate(
            need(!is.null(gene_coldata()), 'Waiting for selection')
        )

        choices <- setdiff(input$facet, c('', 'none', 'gene'))
        if(!is.null(choices)){
            updateSelectInput(session, 'facet_vars',
                              choices=choices)
        } else {
            updateSelectInput(session, 'facet_vars',
                              choices='')
            updateSelectizeInput(session, 'facet_var_levels',
                                 choices='', selected='')
        }
      }, ignoreNULL=FALSE)

      # observer to save initial factor levels of x axis variable
      observeEvent(c(input$xvar, gene_coldata()), {
        validate(
          need(!is.null(input$xvar) & input$xvar != '',
               'Waiting for selection')
        )

        validate(
          need(!is.null(gene_coldata()), 'Waiting for selection')
        )

        validate(
          need(input$xvar %in% colnames(gene_coldata()), 'Waiting for selection')
        )

        # get initial x-axis labels
        cdata <- gene_coldata()[, input$xvar]

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
      observeEvent(c(xchoices$all, input$x_all), {
        xchoices$current <- xchoices$all
      })

      observeEvent(input$x_order, {
        if(length(setdiff(xchoices$current, input$x_order)) != 0)
          xchoices$current <- input$x_order
      })

      # observer to update factor levels when 'Select none'
      # button is clicked
      observeEvent(input$x_none, {
        xchoices$current <- NULL
      })

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
                    input_id=ns('x_order')
                ), # add_rank_list
                add_rank_list(
                    text = 'unused',
                    labels = setdiff(xchoices$all,
                                     xchoices$current),
                    input_id=ns('x_other')
                ) # add_rank_list
              ) # bucket_list
            ) # column
          ) # fluidRow
        ) # tagList
      })

      observeEvent(gene_plot_data$all, {
        validate(
            need(!is.null(gene_plot_data$all), 'Waiting for selection')
        )

        #gene_plot_data$all <- gene_plot_df$df

        coldata <- colData(gene_plot_data$all)

        gene.to.plot <- plot_args()$gene.to.plot

        validate(
            need(input$xvar %in% colnames(coldata), 'Waiting for selection')
        )

        validate(
          need(!is.null(gene.to.plot) & gene.to.plot != '', 'No genes in scratchpad')
        )

        validate(
          need(any(gene.to.plot %in% rownames(gene_plot_data$all)), '')
        )

        # get counts for the genes
        df <- get_gene_counts(gene_plot_data$all, gene.to.plot,
                              input$xvar,
                              norm_method=input$norm_method)

        # check to make sure that most points remain on plot
        # if norm_method is changed
        #
        # NOTE: this is only run when the data set is loaded
        #       to prevent autoscaling when scale is manually changed
        if(!is.null(gene_plot_data$plotted) & !is.na(input$ymax)){
            # add pseudocount
            pseudocount <- config()$server$de_analysis$gene_plot$pseudocount
            df$count <- df$count + pseudocount

            # set initial limits to y_delta % outside max & min
            y_delta <- config()$server$de_analysis$gene_plot$y_delta

            # get min and max init values for y-axis
            y_init <- get_y_init(df, y_delta, pseudocount)

            # get y-limits
            df.max <- max(df$count)
            df.min <- min(df$count)

            # get fraction of points within y-limits
            include.idx <- df$count <= input$ymax & df$count >= input$ymin
            frac <- sum(include.idx)/nrow(df)

            # get fraction of y-limits covered by points
            yrange <- (max(df$count[include.idx]) - min(df$count[include.idx]))/(input$ymax - input$ymin)

            # get thresholds for frac and yrange from config()
            min_frac <- config()$server$de_analysis$gene_plot$min_fraction
            min_yrange <- config()$server$de_analysis$gene_plot$min_yrange_fraction
            if(frac < min_frac | yrange < min_yrange){
                if(frac < min_frac){
                    min_frac <- round(100*(1 - min_frac))
                    showNotification(
                        paste0('More than ', min_frac, '% of points fall',
                               'outside y-limits. Autoscaling y-axis')
                    )
                } else {
                    min_yrange <- round(100*min_yrange)
                    showNotification(
                        paste0('Points cover < ', min_yrange, '% of y-axis',
                               'range. Autoscaling y-axis')
                    )
                }

                updateNumericInput(session, 'ymax', value=y_init[2])
                updateNumericInput(session, 'ymin', value=y_init[1])
            }
        }
      }) # observeEvent update gene_plot y-limits

      # reset y-axis limits
      observeEvent(input$y_auto, {
        validate(
            need(!is.null(gene_plot_data$plotted), 'Waiting for selection')
        )

        df <- gene_plot_data$plotted

        # get pseudocount
        pseudocount <- config()$server$de_analysis$gene_plot$pseudocount

        # set initial limits to y_delta % outside max & min
        y_delta <- config()$server$de_analysis$gene_plot$y_delta

        y_init <- get_y_init(df, y_delta, pseudocount)

        updateNumericInput(session, 'ymax', value=y_init[2])
        updateNumericInput(session, 'ymin', value=y_init[1])

      }) # observeEvent reset y-limits

      all_settings <- reactive({
        list(
          input$ymax,
          input$ymin,
          gene_plot_data$all,
          plot_args(),
          input$plot_do
        )
      })

      # code to get gene plot
      normplot <- eventReactive(all_settings(), {

        # first check to see if some genes were selected
        gene.to.plot <- plot_args()$gene.to.plot
        validate(
          need(!is.null(gene.to.plot) & gene.to.plot != '',
               '\u25c1\u25c1 Please select genes in the scratchpad to visualize them here!')
        )

        # get data set for gene plot
        rld.i <- gene_plot_data$all

        validate(
            need(!is.null(rld.i) & !is.null(gene_coldata()), 'Waiting for selection')
        )

        # get metadata and attach to data
        coldata <- gene_coldata()
        colData(rld.i) <- coldata

        validate(
            need(input$xvar %in% colnames(coldata), 'Waiting for selection')
        )

        # if faceting variable is specified, make sure it is not 'none'
        facet.test <- setdiff(input$facet, c('','none'))
        if(length(facet.test) == 0) facet <- NULL
        else facet <- facet.test

        # if more than 1 facet variable, use first 2
        if(length(facet) > 2){
          showNotification('More than 2 faceting variables specified. Using first two ...',
                           duration=5)
          facet <- facet[seq_len(2)]
        }

        if(!is.null(facet) & !all(facet %in% c('gene'))){
            validate(
              need(input$facet_vars != '' & !is.null(input$facet_var_levels), 'Waiting for selection')
            )

            validate(
              need(input$facet_vars %in% colnames(coldata), 'Waiting for selection')
            )
            validate(
              need(all(input$facet_var_levels %in% coldata[,input$facet_vars]), 'Waiting for selection')
            )
        }

        # get genes to plot
        g <- gene.to.plot

        all_gids <- plot_args()$gene.id
        if(!all(g %in% rownames(rld.i)) | !all(g %in% all_gids)){
            idx <- !g %in% rownames(rld.i) | !g %in% all_gids
            g.missing <- g[idx]
            showNotification(
                paste0('Skipping genes not found in data set: ',
                       paste(g.missing, collapse=',')))
        }

        validate(
          need(any(g %in% rownames(rld.i)) && any(g %in% all_gids), '')
        )

        # get counts for the genes
        df <- get_gene_counts(rld.i, g, input$xvar,
                              norm_method=input$norm_method)

        # add pseudocount
        pseudocount <- config()$server$de_analysis$gene_plot$pseudocount
        df$count <- df$count + pseudocount

        # set initial limits to y_delta % outside max & min
        y_delta <- config()$server$de_analysis$gene_plot$y_delta

        y_init <- get_y_init(df, y_delta, pseudocount)

        # for the first run, ymax is NA
        # This forces an update and reruns this reactive
        if(is.na(input$ymax) | is.null(input$ymax) | input$ymax == 0){
            updateNumericInput(session, 'ymax', value=y_init[2])
            updateNumericInput(session, 'ymin', value=y_init[1])

            validate(
              need(is.na(input$ymax), 'Loading ...')
            )
        }

        # set gene column levels to have same order as input
        df$gene <- factor(df$gene, levels=g)

        # add metadata to count data frame
        df <- add_metadata(df, coldata, cols.to.drop())

        # reset facet levels based on input
        if(!is.null(facet) & !all(facet %in% c('gene'))){
            idx <- df[,input$facet_vars] %in% input$facet_var_levels
            df <- df[idx,]
            df[,input$facet_vars] <- factor(df[, input$facet_vars], levels=input$facet_var_levels)
        }

        # gather params for gene plot
        logy <- input$logy
        freey <- input$freey
        rotate_x_labels <- input$x_rotate
        color <- input$color
        trendline <- input$trendline
        xvar <- input$xvar
        gene_nrow <- config()$ui$de_analysis$gene_plot$nrow
        legend <- input$legend
        ymax <- input$ymax
        ymin <- input$ymin
        ylab <- config()$server$de_analysis$gene_plot$y_labels
        boxes <- input$boxes
        ht <- config()$ui$de_analysis$gene_plot$height

        if(logy){
          validate(
            need(ymin > 0, 'Error: if using log-scale, y-axis min must be > 0')
          )
        }

        # if more than 1 facet variable specified,
        # paste facet columns together
        if(length(facet) > 1){
            facet_pasted <- apply(df[, facet], 1,
                                  function(x) paste(x, collapse=''))
        } else {
            facet_pasted <- df[, facet]
        }

        # find sizes of groups and show warnings if
        # all groups are n=1
        # NOTE: inner table() counts the pasted values
        #       outer table() makes a histogram of counts
        grp_sizes <- table(table(paste(df[,xvar], df[,color],
                                       facet_pasted)))
        if(!all(as.integer(names(grp_sizes)) >= 2) & boxes){
          showNotification('All groups smaller than 2, consider switching off boxes', type='warning')
        }

        if(!is.null(facet)){
          # number of subplots
          plot_num <- length(unique(facet_pasted))

          # fix number of plots per row
          # - free y-axis: 6 (since y-axis labels are drawn for every plot)
          # - otherwise: 8
          if(freey) gene_nrow <- ceiling(plot_num/6)
          else gene_nrow <- ceiling(plot_num/8)

          nrow_def <- config()$ui$de_analysis$gene_plot$nrow
          if(gene_nrow == 0 | is.na(gene_nrow)){
            gene_nrow <- nrow_def
            showNotification(
              'Invalid number of rows. Resetting to default',
              type='warning'
            )
          } else if(gene_nrow > nrow_def*2){
            # adjust height for many rows
            ht <- ht*(gene_nrow/(nrow_def*2))
          }
        }

        # set default x axis order; later take from bucket list
        if(is.null(input$x_order) | !(all(input$x_order %in% xchoices$all))) x_order <- xchoices$all
        else x_order <- input$x_order

        validate(
          need(length(x_order) > 0, 'Must have at least one value to plot on x-axis')
        )

        # save plotted data
        gene_plot_data$plotted <- df

        # reset y-axis limits if most points fall outside
        pts_inside <- sum(df$count >= ymin & df$count <= ymax)
        if(pts_inside > 0 & pts_inside/nrow(df) < 0.75){
          showNotification(
            'More than 25% of points fall outside y-limits, consider auto-scaling'
          )
        }

        validate(
          need(pts_inside > 0, 'No points within y-axis limits. Please adjust limits in "y-axis settings" or click "Autoscale"')
        )

        # if more than 2 faceting variables are specified
        p <- getcountplot(df, intgroup=xvar, ylab=ylab,
                     log=logy, freey=freey,
                     color=color, ymax=ymax, ymin=ymin,
                     factor.levels=x_order, rotate_x_labels=rotate_x_labels,
                     nrow=gene_nrow, trendline=trendline,
                     facet=facet, legend=legend, boxes=boxes)

        if(input$txt_scale == 0 | is.na(input$txt_scale)){
          txt_scale <- 1
          showNotification(
            'Invalid text scaling factor. Resetting to default',
            type='warning'
          )
        } else {
          txt_scale <- input$txt_scale
        }

        p <- p + theme(axis.text.x=element_text(size=9*txt_scale),
                       axis.text.y=element_text(size=8*txt_scale),
                       strip.text.x=element_text(size=9*txt_scale),
                       legend.title=element_text(size=10*txt_scale),
                       legend.text=element_text(size=9*txt_scale))

        gene_plot_data$handle <- p

        ggplotly(p, height=ht)
      }) # eventReactive normplot

      output$normplot <- renderPlotly({
        normplot()
      })

      ################## buttons #####################

      helpButtonServer('de_geneplt_help', size='l')
      helpButtonServer('geneplt_opts_help', size='l')

      downloadButtonServer('normplot_download', reactive({ gene_plot_data$handle }),
                           'gene_plot')

    } # function
  ) # moduleServer
}
