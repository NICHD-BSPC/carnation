#' Carnation
#'
#' Interactive shiny dashboard for exploring RNA-Seq analysis.
#'
#' @param credentials path to encrypted sqlite db with user credentials.
#' @param passphrase passphrase for credentials db.
#' @param enable_admin if TRUE, admin view is shown. Note, this is only available
#'        if credentials have sqlite backend.
#' @param ... parameters passed to shinyApp() call
#'
#' @export
run_carnation <- function(credentials=NULL, passphrase=NULL, enable_admin=TRUE, ...){

  # read config yaml
  config <- get_config()

  # set some options
  oopt <- options(spinner.type = 4)
  options(shiny.maxRequestSize = config$max_upload_size*1024^2)

  # reset to previous options on exit
  on.exit(options(oopt))

  if(!is.null(credentials)){

    # check to see if shinymanager is available
    if(!requireNamespace('shinymanager', quietly=TRUE)){
      stop(
        paste('Login functionality using SQL/sqlite credentials requires "shinymanager".',
              'Please install using "install.packages(\'shinymanager\')"'),
        .call=FALSE
      )
    } else if(!file.exists(credentials)){
      stop(
        paste0('Credentials specified, but file not found: "',
               credentials, '"')
      )
    }
  }

  ui <- fluidPage(
    theme=shinytheme('united'),

    introjsUI(),

    # custom CSS styles
    tags$head(
      tags$style(
        HTML(config$style$global)
      )  # tags$style
    ), # tags$head

    titlePanel(
      fluidRow(
        # add spacer to center heading
        column(4, span()),
        column(4,
          tags$div(
            HTML(
              paste0(
                'ca',
                span('rna', style='color: #e95420;'), # primary color from united theme
                'tion'
              ) # paste
            ) # HTML
          ),
          align='center',
          style='font-family: Helvetica; font-size: 40px;'
        ), # column
        column(2,
          actionButton('intro', label='Take a tour!',
                       icon=icon('info'))
        ), # column

        column(2,
          introBox(
            saveUI('save_object'),
            data.step=7,
            data.intro='Use this button at any point to save the modified object'
          ), # introBox
          offset=-1
        ) # column

      ), # fluidRow
      windowTitle='Carnation'
    ), # titlePanel

    # side bar
    sidebarPanel(width=1, id='sidebar',

      introBox(
        fluidRow(align='center',
          introBox(
            dropdownButton(
              tagList(
                  fluidRow(
                    column(6, align='left',
                      tags$label(class='control-label',
                                      'DE Filters'
                      ) # tags$label
                    ), # column
                    column(6, align='right',
                      helpButtonUI('de_filters_help')
                    ) # column
                  ), # fluidRow

                  fluidRow(
                    column(4, h5('FDR threshold')),
                    column(8,
                      numericInput("fdr.thres", label=NULL,
                        value=config$ui$de_analysis$filters$fdr_threshold,
                        min=0, max=1
                      ) # numericInput
                    ) # column
                  ), # fluidRow
                  fluidRow(
                    column(4, h5('log2FC threshold')),
                    column(8,
                      numericInput("fc.thres", label=NULL,
                        value=config$ui$de_analysis$filters$log2fc_threshold,
                        min=0
                      ) # numericInput
                    ) # column
                  ) # fluidRow

              ), # tagList

              icon = icon("sliders-h"), width = "400px",
              size='sm',

              tooltip = tooltipOptions(title = "Global settings")

            ), # dropdownButton
            data.step=3,
            data.intro='Click this button to access global settings'
          ), # introBox

          br(),

          introBox(
            dropdownButton(
              tagList(

                conditionalPanel("input.mode == 'Load data'",
                  span('Load your data here. You can choose existing or new datasets',
                    style='font-size: 15px;'
                  )
                ), # conditionalPanel

                conditionalPanel("input.mode == 'DE analysis' & input.de_mode == 'Summary'",
                  span('Here we show a summary of the DE comparisons performed',
                    style='font-size: 15px;'
                  )
                ), # conditionalPanel

            conditionalPanel("input.mode == 'DE analysis' & input.de_mode == 'Table'",

              fluidRow(
                column(6, align='left',
                  tags$label(class='control-label',
                             'Comparison')
                ), # column
                column(6, align='right',
                  helpButtonUI('de_cmp_help')
                ) # column
              ), # fluidRow

              selectizeInput('comp_all',
                             label=NULL,
                             choices=NULL,
                             selected=NULL,
                             multiple=TRUE),

              # For 'Select all' and 'Select none' buttons
              fluidRow(style='margin-bottom: 20px',
                column(6, align='center',
                  actionButton('select_all_comp', 'Select all', class = "btn-secondary")
                ), # column for 'Select all'
                column(6, align='center',
                  actionButton('select_none_comp', 'Select none', class = "btn-secondary")
                ) # column for 'Select none'
              ), # fluidRow

              fluidRow(
                column(4, h5('Only DE genes')),
                column(8,
                  checkboxInput("toggle_filters", label=NULL,
                    value=config$ui$de_analysis$filters$only_de_toggle
                  ) # checkboxInput
                ) # column
              ) # fluidRow

            ), # conditionalPanel

            ########################## Metadata ######################################

            conditionalPanel("input.mode == 'DE analysis' & input.de_mode == 'Metadata'",

              fluidRow(
                column(6, align='left',
                  tags$label(class='control-label',
                             'Sample group')
                ), # column
                column(6, align='right',
                  helpButtonUI('metadata_sidebar_help')
                ) # column
              ), # fluidRow

              metadataUI('metadata', panel='sidebar')

            ), # conditionalPanel

            ########################## PCA plot ######################################

            conditionalPanel("input.mode == 'DE analysis' & input.de_mode == 'PCA plot'",
              pcaPlotUI('pcaplot', panel='sidebar')
            ), # conditionalPanel

            ######################### MA plot ########################################

            conditionalPanel("input.mode == 'DE analysis' & input.de_mode == 'MA plot'",

              maPlotUI('maplot', panel='sidebar')

            ), # conditionalPanel

            ######################### Scatter plot ########################################

            conditionalPanel("input.mode == 'DE analysis' & input.de_mode == 'Scatter plot'",

              scatterPlotUI('scatterplot', panel='sidebar')

            ), # conditionalPanel

            ######################### Upset plot ####################################

            conditionalPanel("input.mode == 'DE analysis' & input.de_mode == 'Upset plot'",

              upsetPlotUI('upset_plot', panel='sidebar')

            ), # conditionalPanel

            ######################### Heatmap ########################################

            conditionalPanel("input.mode == 'DE analysis' & input.de_mode == 'Heatmap'",

              heatmapUI('heatmap', panel='sidebar')

            ), # conditionalPanel

            ########################## Functional enrichment #########################

            conditionalPanel("input.mode == 'Functional enrichment' & input.func == 'Table'",
              enrichUI('func_enrich',
                       panel='sidebar',
                       tab='table')
            ), # conditionalPanel

            conditionalPanel("input.mode == 'Functional enrichment' & input.func == 'Plots'",
              enrichUI('func_enrich',
                       panel='sidebar',
                       tab='plots')
            ), # conditionalPanel

            # compare results
            conditionalPanel("input.mode == 'Functional enrichment' & input.func == 'Compare results'",
              enrichUI('func_enrich',
                       panel='sidebar',
                       tab='compare_results')
            ), # conditionalPanel

            ########################## Pattern analysis ##############################

            conditionalPanel(condition = "input.mode == 'Pattern analysis'",

              patternPlotUI('deg_plot', 'sidebar', 'both')

            ), # conditionalPanel

            conditionalPanel(condition = "input.mode == 'Pattern analysis' & input.pattern_mode == 'Plot'",

              patternPlotUI('deg_plot', 'sidebar', 'plot')

            ), # conditionalPanel

            conditionalPanel(condition = "input.mode == 'Pattern analysis' & input.pattern_mode == 'Cluster membership'",

              patternPlotUI('deg_plot', 'sidebar', 'table')

            ), # conditionalPanel

            ########################## Gene plot #####################################

            conditionalPanel("input.mode == 'DE analysis' & input.de_mode == 'Gene plot'",

              genePlotUI('gene_plot', panel='sidebar')

            ), # conditionalPanel

            ########################## Settings ######################################

            conditionalPanel(condition = "input.mode == 'Settings'",
              uiOutput('settings_sidebar')
            ) # conditionalPanel

              ), # tagList

              icon = icon("gear"), width = "400px",
              size='sm',

              tooltip = tooltipOptions(title = "Settings")

            ), # dropdownButton

            data.step=4,
            data.intro='Here you will find specific controls to adjust plots or tables'
          ), # introBox

          br(),

          ########################## Gene Scratchpad ###############################

          introBox(
            dropdownButton(
              tagList(
                conditionalPanel("input.mode == 'DE analysis' | input.mode == 'Functional enrichment' | input.mode == 'Pattern analysis'",

                  fluidRow(
                    column(6, strong('Gene scratchpad')),
                    column(6, align='right',
                       helpButtonUI('gene_scratchpad_help')
                    ) # column
                  ), # fluidRow

                  selectizeInput('gene.to.plot',
                                 label=NULL,
                                 choices=NULL, selected=NULL,
                                 options=list(create=TRUE, delimiter=','),
                                 multiple=TRUE
                  ), # selectizeInput

                  fluidRow(
                    column(12, align='center',
                      style='margin-bottom: 20px;',
                      splitLayout(cellWidths=c('50%','50%'),
                        actionButton('reset.genes', label='Reset',
                                   class='btn-primary'),
                        actionButton('quick_add', label='Add top genes')
                      ) # splitLayout
                    ) # column
                  ), # fluidRow

                  bsCollapse(open='settings',
                    bsCollapsePanel(span(icon('gear'), 'Top genes settings'),
                      value='settings',

                      fluidRow(
                        column(4, 'Comparison'),
                        column(8,
                          selectInput('scratchpad_comp',
                                      label=NULL,
                                      choices=NULL, selected=NULL)
                        ) # column
                      ), # fluidRow

                      fluidRow(
                        column(4, '# of genes'),
                        column(8,
                          numericInput('scratchpad_ngenes',
                                      label=NULL,
                                      value=config$server$de_analysis$gene_scratchpad$ngenes,
                                      step=1)
                        ) # column
                      ), # fluidRow

                      fluidRow(
                        column(4, 'ranking metric'),
                        column(8,
                          selectInput('top_genes_by',
                                      label=NULL,
                                      choices=c('padj', 'log2FoldChange'))
                        ) # column
                      ) # fluidRow

                    ) # bsCollapsePanel
                  ) # bsCollapse

                ) # conditionalPanel

              ), # tagList

              icon = icon("clipboard"), width = "400px",
              size='sm',

              tooltip = tooltipOptions(title = "Gene scratchpad")

            ), # dropdownButton
            data.step=6,
            data.intro='Keep track of your favorite genes or quickly select top genes with the "Gene scratchpad" here'
          ), # introBox

          data.step=2,
          data.intro='You can use controls shown in this area to filter data or adjust plots and tables.'
        ) # fluidRow
      ) # introBox
    ), # sidebarPanel

    mainPanel(width=11,

      introBox(
        tabsetPanel(type='tabs',
                    id='mode',

          tabPanel('Load data',
            selectInput('data_type', label='Type of data',
                        choices=c('Existing', 'New', 'Edit')),

            conditionalPanel('input.data_type == "Existing"',
                selectizeInput('dds',
                               label=h5('Available projects'),
                               choices=NULL,
                               selected=NULL
                ), # selectizeInput
                selectizeInput('assay',
                               label=h5('Available analyses'),
                               choices=NULL,
                               selected=NULL
                ), # selectizeInput
                actionButton('assay_do', label='Go!',
                             class='btn-primary')
            ), # conditionalPanel

            withSpinner(
              uiOutput('load_ui')
            ) # withSpinner
          ), # tabPanel

          tabPanel('DE analysis',
            tabsetPanel(id='de_mode',

              ##### DE Analysis subtabs #####

              tabPanel('Summary',
                fluidRow(
                  column(1, align='left',
                    helpButtonUI('de_summary_help')
                  ) # column
                ), # fluidRow

                fluidRow(
                  column(12, align='left',
                    withSpinner(
                      DTOutput('summary_tbl')
                    ) # withSpinner
                  ) # column
                ) # fluidRow
              ), # tabPanel summary

              tabPanel('Metadata',
                # help button
                fluidRow(
                  column(12, align='left',
                    helpButtonUI('de_meta_help')
                  ) # column
                ), # fluidRow

                metadataUI('metadata', panel='main')
              ), # tabPanel metadata

              tabPanel('PCA plot',
                pcaPlotUI('pcaplot', panel='main')
              ), # tabPanel PCA

              tabPanel('Table',
                fluidRow(
                  column(12, align='left',
                    helpButtonUI('de_tbl_help')
                  ) # column
                ), # fluidRow

                fluidRow(
                  column(12,
                    withSpinner(
                      DTOutput('detable')
                    ) # withSpinner
                  ) # column
                ), # fluidRow

                fluidRow(align='center', style='margin-top: 25px;',
                  column(12,
                    splitLayout(
                      cellWidths=c('15%', '25%', '20%'),

                      strong('Selection options'),
                      actionButton('add_selected', 'Add to scratchpad'),
                      actionButton('reset_detable',
                                   'Reset selection',
                                   class='btn-primary')
                    ) # splitLayout
                  ) # column
                ) # fluidRow
              ), # tabPanel table

              tabPanel('MA plot',
                maPlotUI('maplot', panel='main')
              ), # tabPanel ma_plot

              tabPanel('Scatter plot',
                scatterPlotUI('scatterplot', panel='main')
              ), # tabPanel scatterplot

              tabPanel('Gene plot',
                genePlotUI('gene_plot', panel='main')
              ), # tabPanel gene_plot

              tabPanel('Upset plot',
                upsetPlotUI('upset_plot', panel='main')
              ), # tabPanel upset

              tabPanel('Heatmap',
                heatmapUI('heatmap', panel='main')
              ) # tabPanel

            ) # tabsetPanel de_mode
          ), # tabPanel de_analysis


          tabPanel('Functional enrichment',
            tabsetPanel(id='func',

              ###### Functional enrichment subtabs #####

              tabPanel('Table',
                enrichUI('func_enrich',
                         panel='main',
                         tab='table')
              ), # tabPanel table

              tabPanel('Plots',
                enrichUI('func_enrich',
                         panel='main',
                         tab='plots')
              ), # tabPanel plots

              tabPanel('Compare results',
                enrichUI('func_enrich',
                         panel='main',
                         tab='compare_results')
              ) # tabPanel compare_results

            ) # tabsetPanel func
          ), # tabPanel functional_enrichment

          tabPanel('Pattern analysis',
            tabsetPanel(id='pattern_mode',

              tabPanel('Plot',
                patternPlotUI('deg_plot', 'plot')
              ), # tabPanel plot

              tabPanel('Cluster membership',
                patternPlotUI('deg_plot', 'table')
              ) # tabPanel cluster_membership

            ) # tabsetPanel pattern_mode
          ), # tabPanel pattern_analysis

          tabPanel('Settings',
            uiOutput('settings_main')
          ) # tabPanel settings

        ), # tabsetPanel mode
        data.step=1,
        data.intro='Choose project and analysis and click "Go".\n\nOnce the data is loaded, use tabs to explore and navigate'
      ) # introBox
    ) # mainPanel

  ) # ui

  # add custom login page, shown only if credentials are specified
  if(!is.null(credentials)){
    ui <- secure_app(
            ui,

            # customize login page
            tags_top =
              tags$div(
                HTML(
                  paste0(
                  'ca',
                  span('rna', style='color: #e95420;'), # primary color from united theme
                  'tion'
                  )
                ),
                align='center',
                style='font-family: Helvetica; font-size: 40px;'


            ),
            # add information on bottom ?
            tags_bottom = tags$div(
              tags$p(
                "For any question, please  contact carnation authors",
              )
            ),

            enable_admin=enable_admin,
            theme=shinytheme('united')
          )
  }

  server <- function(input, output, session){

    oopt <- options(spinner.type = 4)

    # get admin group from config
    admin_group <- config$server$admin_group

    # metadata columns to ignore
    cols.to.drop <- config$server$cols.to.drop

    # list to hold dds, rld and id mapping
    app_object <- reactiveValues(dds=NULL, rld=NULL, res=NULL,
                                 dds_mapping=NULL, labels=NULL,
                                 enrich=NULL, genetonic=NULL,
                                 degpatterns=NULL,
                                 all_dds=NULL, all_rld=NULL)

    # list to hold original object and file path
    original <- reactiveValues(obj=NULL, path=NULL)

    # list to hold user details
    user_details <- reactiveValues(username=NULL, admin=FALSE)

    #################### authentication ####################

    ########## shinymanager login #############

    res_auth <- secure_server(
      check_credentials = check_credentials(
          db=credentials,
          passphrase=passphrase
      )
    )

    observeEvent(res_auth, {

      req(res_auth$user)

      user_details$username <- res_auth$user
      user_details$admin <- res_auth$user_info
    })

    ########## proxy login ###########

    # get username from http request header
    observeEvent(session$request, {
      user_details$username <- session$request[[ config$http_request_header ]]
    })

    #################### Intro ####################

    observeEvent(input$intro, {

      introjs(session,
              options = list("nextLabel"="Next",
                             "prevLabel"="Back")
              )

    }) # observeEvent

    ########### Load data #####################

    output$load_ui <- renderUI({
      if(input$data_type == 'Edit'){
        if(is.null(original$obj)){
          showNotification(
            'Must load analysis to edit!',
            type='error'
          )
        }

        validate(
          need(!is.null(original$obj), 'Must load analysis to edit!')
        )
        loadDataUI('edit_obj')
      } else if(input$data_type == 'New'){
        loadDataUI('load_new_data')
      }
    })

    reload_new <- loadDataServer('load_new_data',
                                 username=reactive({ user_details$username }))

    observeEvent(reload_new(), {
      flag <- reload_new()

      if(flag) session$reload()
    })

    reload_edit <- loadDataServer('edit_obj',
                                  username=reactive({ user_details$username }),
                                  rds=original)

    observeEvent(reload_edit(), {
      flag <- reload_edit()

      if(flag) session$reload()
    })

    ########### Settings module ###############

    pattern <- config$server$pattern

    # reactive values to keep assay list
    assay.list <- reactiveValues(l=NULL)

    settings <- settingsServer('settings',
                               details=reactive({ list(username=user_details$username, where=input$shinymanager_where) }),
                               depth=2,
                               end_offset=0,
                               assay_fun=function(x)
                                 sub(paste0(pattern, '\\.rds$'), '',
                                     basename(x),
                                     ignore.case=TRUE)
                               )

    # update assay list
    observeEvent(settings(), {
      l <- settings()

      validate(
         need(!is.null(l$assay_list), 'No projects found')
      )
      # get project names and alphabetically sort
      project_names <- names(l$assay_list)

      # if more than 1 project, add 'Choose one' to options
      dds_choices <- project_names[order(project_names)]
      if(length(dds_choices) > 1)
        dds_choices <- c('Choose one', dds_choices)
      updateSelectizeInput(session, 'dds',
                           choices=dds_choices)
      assay.list$l <- l$assay_list

      if(l$reload_parent) session$reload()

      # hide admin tab if user not in admin group
      if(!l$is_admin){
        hideTab(inputId = 'mode', target='Settings')
      }
    })

    output$settings_main <- renderUI({
      settingsUI('settings', panel='main', username=reactive({ user_details$username }))
    })

    output$settings_sidebar <- renderUI({
      settingsUI('settings', panel='sidebar', username=reactive({ user_details$username }))
    })

    ############### Metadata #################

    # reactive values to hold metadata
    coldata.all <- reactiveValues(init=NULL, curr=NULL, staging=NULL,
                                  gene=NULL)

    # get metadata from module
    metadata <- metadataServer('metadata', app_object,
                               cols.to.drop)

    observeEvent(c(app_object$dds, metadata()), {
      clist <- metadata()

      for(name in names(clist)){
          coldata.all[[ name ]] <- clist[[ name ]]
      }
    }) # observeEvent

    ############# Save object button ##################

    save_event <- saveServer('save_object',
                             original=original,
                             current=app_object,
                             coldata=coldata.all$curr,
                             pattern=pattern,
                             username=reactive({ user_details$username }))

    observeEvent(save_event(), {
      validate(
          need(!is.null(original$obj), 'Waiting for selection')
      )

      l <- save_event()

      if(l$reload_parent) session$reload()

    }) # observeEvent

    # function to reset app data
    reset_data <- function(){
      # reset data
      original$obj <- NULL
      original$path <- NULL
      app_object$dds <- NULL
      app_object$rld <- NULL
      app_object$res <- NULL
      app_object$enrich <- NULL
      app_object$genetonic <- NULL
      app_object$degpatterns <- NULL
      app_object$labels <- NULL
      app_object$all_dds <- NULL
      app_object$all_rld <- NULL
      app_object$dds_mapping <- NULL

      # reset coldata
      coldata.all$init <- NULL
      coldata.all$curr <- NULL
      coldata.all$staging <- NULL
      coldata.all$gene <- NULL

      # reset gene ids
      gene.id$gene <- NULL

      # reset de table
      res_data$tbl <- NULL
    }

    ############### Initial load #################

    # update assay list & reset data
    observeEvent(c(input$dds, assay.list$l), {
      l <- assay.list$l

      validate(
        need(!is.null(input$dds), '')
      )

      assay.choices <- l[[input$dds]]

      # if more than one assay found & autoload_first_analysis not set
      if(length(assay.choices) > 1 & !config$server$autoload_first_analysis){
        assay.choices <- c('Choose one', assay.choices)
      }

      updateSelectizeInput(session, 'assay',
                           choices=assay.choices,
                           selected=assay.choices[1])

      reset_data()
    }) # observeEvent

    # observer to load data
    observeEvent(input$assay_do, {
      validate(
        need(!is.null(input$assay) & input$assay != '' & input$assay != 'Choose one',
             'No assay selected!')
      )

      reset_data()

      # get file size for message
      fs <- file.size(input$assay)
      if(fs < (1024**3)){
        fs_human <- paste0('(', round(fs/(1024**2), digits=1), ' MB)')
      } else {
        fs_human <- paste0('(', round(fs/(1024**3), digits=1), ' GB)')
      }
      showModal(
        modalDialog(
          span(paste('Loading data', fs_human, 'please wait')),
          footer=NULL
        )
      )

      # load data from RDS object
      obj <- readRDS(input$assay)
      original$obj <- obj
      original$path <- input$assay

      # if 'carnation-ready' object is uploaded these elements will always be present
      all_obj_names <- c('res', 'dds', 'rld', 'labels', 'dds_mapping')
      if(all(all_obj_names %in% names(obj))){
        app_object$res <- obj$res
        app_object$dds <- obj$dds
        app_object$rld <- obj$rld
        app_object$labels <- obj$labels
        app_object$dds_mapping <- obj$dds_mapping
        if('all_dds' %in% names(obj)){
          app_object$all_dds <- obj$all_dds
        } else {
          app_object$all_dds <- NULL
        }

        if('all_rld' %in% names(obj)){
          app_object$all_rld <- obj$all_rld
        } else {
          app_object$all_rld <- NULL
        }

        # check for the optional 'enrich' & 'degpatterns' slots
        if('enrich' %in% names(obj)){
          if(!'genetonic' %in% names(obj)){
            stop('Object is missing "genetonic" slot')
          }
          app_object$enrich <- obj$enrich
          app_object$genetonic <- obj$genetonic
        }
        if('degpatterns' %in% names(obj)) app_object$degpatterns <- obj$degpatterns
      } else {

      # perform vst if necessary
      orig.names <- names(obj)

      # TODO: add checks for presence of expected elements
      # res, dds, enrich, degpatterns
      res.name <- orig.names[grep('res', orig.names)]
      dds.name <- orig.names[setdiff(grep('dds', orig.names),
                               c(grep('all_dds', orig.names),
                                 grep('dds_mapping', orig.names)))]
      rld.name <- orig.names[setdiff(grep('rld', orig.names),
                                     grep('all_rld', orig.names))]
      enrich.name <- orig.names[grep('enrich', orig.names)]
      degpatterns.name <- orig.names[grep('degpatterns', orig.names)]

      if(length(res.name) == 0 | length(dds.name) == 0){
        if(length(res.name) == 0 & length(dds.name) == 0){
          missing <- 'both'
        } else if(length(res.name) == 0){
          missing <- 'DE result'
        } else {
          missing <- 'counts table'
        }
        showNotification(
          paste0('Loaded object must have both DE results & counts table. Missing:',
                 missing),
          type='error'
        )

        validate(
          need(length(res.name) != 0 & length(dds.name) != 0,
               'Missing DE result and/or counts table')
        )
      }

      validate(
          need(length(dds.name) == 1,
               'Uploaded object must have exactly 1 counts list')
      )

      # TODO: sanitize res, dds
      # TODO: sanitize degpatterns
      # - check for 'normalized' in names if list
      # - check for 'genes' column

      if(length(rld.name) == 0){

        showModal(
          modalDialog(
            span('Normalizing data'),
            footer=NULL
          )
        )

        rld_list <- lapply(obj[[dds.name]],
                           function(x) varianceStabilizingTransformation(x, blind=TRUE))

        obj$rld_list <- rld_list
        rld.name <- 'rld_list'

        removeModal()
      } else if(length(rld.name) > 1){
        validate(
            need(length(rld.name) == 1,
                 'Uploaded object must have exactly 1 normalized counts list')
        )
      }

      showModal(
        modalDialog(
          span('Sanitizing object'),
          footer=NULL
        )
      )

      # add 'sample' to colData of rld_list & dds_list
      obj[[ rld.name ]] <- lapply(obj[[ rld.name ]],
                             function(x){
                               colData(x)$sample <- rownames(colData(x))
                               x
                           })

      obj[[ dds.name ]] <- lapply(obj[[ dds.name ]],
                             function(x){
                               colData(x)$sample <- rownames(colData(x))
                               x
                           })

      # add obj slots to reactive values
      obj <- make_final_object(obj)

      # get final names of obj
      n <- names(obj)
      res.name <- n[grep('res', n)]
      enrich.name <- n[grep('enrich', n)]
      degpatterns.name <- n[grep('degpatterns', n)]

      # get dds and rld elements
      # - need to distinguish from all_dds, dds_mapping & all_rld
      dds.idx <- setdiff(grep('dds', n),
                         c(grep('all_dds', n),
                           grep('dds_mapping', n)))
      rld.idx <- setdiff(grep('rld', n), grep('all_rld', n))
      dds.name <- n[dds.idx]
      rld.name <- n[rld.idx]

      # map to elements
      app_object$dds <- obj[[dds.name]]
      app_object$rld <- obj[[rld.name]]
      app_object$res <- obj[[res.name]]
      if(length(enrich.name) != 0)
        app_object$enrich <- obj[[enrich.name]]
      if(length(degpatterns.name) != 0)
        app_object$degpatterns <- obj[[degpatterns.name]]
      app_object$labels <- obj$labels
      app_object$dds_mapping <- obj$dds_mapping

      # add element with genetonic objects
      # NOTE: using loop to keep names
      if(!'genetonic' %in% names(obj) | is.null(obj$genetonic)){

        showModal(
          modalDialog(
            span('Converting FE results to GeneTonic format'),
            footer=NULL
          )
        )

        start_time <- Sys.time()

        # flatten list
        elem_names <- NULL
        sep <- '*'
        res_keys <- list()
        for(x in names(app_object$enrich)){
          for(y in names(app_object$enrich[[x]])){
            # NOTE: if key is 'res' save & skip
            if(y == 'res'){
              res_keys[[ x ]] <- app_object$enrich[[ x ]][[ 'res' ]]
              next
            }
            for(z in names(app_object$enrich[[x]][[y]])){
              elem_names <- c(elem_names, paste(x, y, z, sep=sep))
              if(!is.data.frame(app_object$enrich[[ x ]][[ y ]][[ z ]])){
                app_object$enrich[[ x ]][[ y ]][[ z ]] <- app_object$enrich[[ x ]][[ y ]][[ z ]]@result
              }
            }
          }
        }
        names(elem_names) <- elem_names

        # run conversion
        #
        # NOTE: making local copies here to avoid 'reactive' errors
        #       in parallel lapply
        #
        res_names <- names(app_object$res)
        res_list <- app_object$res
        enrich_list <- app_object$enrich

        # TODO: add check for cores
        flat_obj <- BiocParallel::bplapply(elem_names, function(x){
                      toks <- strsplit(x, split=sep, fixed=TRUE)[[1]]
                      # look for toks[1] in res_names & res_keys, else NULL
                      if(toks[1] %in% res_names){
                        res <- res_list[[ toks[1] ]]
                      } else if(toks[1] %in% names(res_keys)){
                        res <- res_list[[ res_keys[[ toks[1] ]] ]]
                      } else{
                        return(NULL)
                      }

                      eres <- enrich_list[[ toks[1] ]][[ toks[2] ]][[ toks[3] ]]

                      # NOTE: if obj is dataframe then convert, else continue
                      if(is.data.frame(eres)){
                        # NOTE: if df contains 'core_enrichment' then it's a 'gseaResult'
                        #       else it's an enrichResult
                        if(!'core_enrichment' %in% colnames(eres)){
                          eres <- makeEnrichResult(eres, type='enrichResult')
                        } else {
                          eres <- makeEnrichResult(eres, type='gseaResult')
                        }
                      }

                      df <- enrich_to_genetonic(eres, res)
                      df
                    }, BPPARAM=BiocParallel::MulticoreParam(config$server$cores))

        # reconstitute & clean up
        app_object$genetonic <- list()
        rm(res_list)
        rm(enrich_list)

        for(x in names(app_object$enrich)){
          app_object$genetonic[[x]] <- list()

          for(y in names(app_object$enrich[[x]])){
            # NOTE: if key is 'res', skip
            if(y == 'res'){
              next
            } else {
              app_object$genetonic[[x]][[y]] <- list()
            }

            for(z in names(app_object$enrich[[x]][[y]])){
              key <- paste(x, y, z, sep=sep)
              app_object$genetonic[[x]][[y]][[z]] <- flat_obj[[key]]
            }
          }
        }

        end_time <- Sys.time()

        delta <- end_time - start_time

      } else {
        app_object$genetonic <- obj[['genetonic']]
      }

      if('all_dds' %in% names(obj)){
          app_object$all_dds <- obj$all_dds
          app_object$all_rld <- obj$all_rld
      } else {
          app_object$all_dds <- NULL
          app_object$all_rld <- NULL
      }

        showNotification(
          'Remember to save object to avoid preprocessing again next time!',
          type='warning'
      )
      }

      if(is.null(app_object$all_dds)){
          gene.id$gene <- rownames(app_object$dds[[1]])
      } else {
          gene.id$gene <- rownames(app_object$all_dds)
      }

      removeModal()

      updateTabsetPanel(session, inputId='mode',
                        selected='DE analysis')
      updateTabsetPanel(session, inputId='de_mode',
                        selected='Summary')
    }) # observeEvent load data

    # update comparison menus after load
    observeEvent(app_object$res, {
      validate(
        need(!is.null(app_object$res) & !is.null(app_object$dds_mapping), 'Waiting for selection')
      )
      # update names of comparisons for all tabs
      updateSelectizeInput(session,
                           'comp_all',
                           choices=names(app_object$res),
                           selected=names(app_object$res)[1])

      updateSelectizeInput(session,
                           'scratchpad_comp',
                           choices=names(app_object$res))

    }) # observeEvent update comparisons

    ############### Gene scratchpad #################

    # gene to plot
    gene.id <- reactiveValues(gene=NULL)

    # observer for adding top genes
    observeEvent(input$quick_add, {
      validate(
        need(!is.null(input$scratchpad_comp) & input$scratchpad_comp != '' & !is.null(app_object$res) & input$scratchpad_comp %in% names(app_object$res),
             'Waiting for selection')
      )

      validate(
        need(!is.null(gene.id$gene), 'Waiting for selection')
      )

      fc.thres <- ifelse(input$fc.thres == '' | is.na(input$fc.thres),
                         config$ui$de_analysis$filters$log2fc_threshold, input$fc.thres)
      fdr.thres <- ifelse(input$fdr.thres == '' | is.na(input$fdr.thres),
                          config$ui$de_analysis$filters$fdr_threshold, input$fdr.thres)

      res <- app_object$res[[input$scratchpad_comp]]

      validate(
        need(input$scratchpad_ngenes > 0, '')
      )
      # get top genes
      g.top <- top.genes(res,
                         fdr.thres=fdr.thres,
                         fc.thres=fc.thres,
                         n=input$scratchpad_ngenes,
                         by=input$top_genes_by)

      showNotification(paste0('Adding ', length(g.top),
                              ' genes to scratchpad'))

      updateSelectizeInput(session, 'gene.to.plot',
                           choices=gene.id$gene,
                           selected=unique(c(input$gene.to.plot, g.top)),
                           server=TRUE)

    }) # observeEvent add_top_genes by lfc

    # observer to reset gene scratchpad
    observeEvent(input$reset.genes, {
      validate(
        need(!is.null(gene.id$gene), 'Waiting for selection')
      )

      updateSelectizeInput(session, 'gene.to.plot',
                           choices=gene.id$gene,
                           selected='',
                           server=TRUE)

    }) # observeEvent

    # observer to update gene scratchpad choices
    observeEvent(gene.id$gene, {
      validate(
        need(!is.null(gene.id$gene), 'Waiting for selection')
      )

      # update gene selector with genes in object
      updateSelectizeInput(session, 'gene.to.plot',
                           choices=gene.id$gene,
                           selected='',
                           server=TRUE)
    }) # observeEvent

    ############### Summary table #################

    # reactive function to return df with summary of results
    get_summary <- eventReactive(c(app_object$res, input$fdr.thres, input$fc.thres),{
      validate(
        need(!is.null(app_object$res) & !is.null(app_object$dds_mapping), 'Waiting for selection')
      )

      validate(
        need(all(names(app_object$dds_mapping) %in% names(app_object$res)), 'Waiting for selection')
      )

      fc.thres <- ifelse(input$fc.thres == '' | is.na(input$fc.thres),
                         config$ui$de_analysis$filters$log2fc_threshold, input$fc.thres)
      fdr.thres <- ifelse(input$fdr.thres == '' | is.na(input$fdr.thres),
                          config$ui$de_analysis$filters$fdr_threshold, input$fdr.thres)

      df <- summarize.res.list(app_object$res, app_object$dds, app_object$dds_mapping,
                               alpha=fdr.thres,
                               lfc.thresh=fc.thres,
                               app_object$labels)
      df
    }) # eventReactive get_summary

    # summary table output
    output$summary_tbl <- renderDT({

      df <- get_summary()
      if(!is.null(app_object$labels)){
          df %>%
              relocate(.data$description, .after=.data$down) %>%
              datatable(rownames=FALSE,
                        selection='none',
                        options=list(autoWidth=TRUE)) %>%
              formatStyle(c('description', 'design', 'contrast'), 'white-space'='nowrap')
      } else {
          df %>%
              datatable(rownames=FALSE,
                        selection='none',
                        options=list(autoWidth=TRUE)) %>%
              formatStyle(c('design', 'contrast'), 'white-space'='nowrap')
      }
    }) # renderDT

    ############### DE table #################

    # reactiveValues to keep track of results tbl data
    res_data <- reactiveValues(tbl=NULL)

    # Reactive expression to update res_data
    processed_table <- eventReactive(
      c(input$comp_all,
        input$fdr.thres,
        input$fc.thres,
        input$toggle_filters
      ), {

      validate(
        need(!is.null(app_object$res),
             'Waiting for selection')
      )

      if (is.null(input$comp_all) || length(input$comp_all) == 0) {
        res_data$tbl <- NULL

      } else {

        if(length(input$comp_all) > 1){
          # get shared & all columns for multi-comparison
          shared_columns <- NULL
          all_columns <- NULL

          for(comp in input$comp_all){
            cnames <- colnames(app_object$res[[ comp ]])
            if(is.null(shared_columns)){
              shared_columns <- cnames
              all_columns <- cnames
            } else {
              shared_columns <- intersect(shared_columns, cnames)
              all_columns <- unique(c(all_columns, cnames))
            }
          }

          validate(
            need(length(shared_columns) > 0,
                 'Error: No columns shared between tables!')
          )

          # output notif about missing cols
          missing_cols <- setdiff(all_columns, shared_columns)
          if(length(missing_cols) > 0){
            msg <- paste(missing_cols, collapse=', ')
            showNotification(
              paste('Warning: Dropping columns not present in all tables -',
                    msg),
              type='warning'
            )
          }

          # Loop through each selected comparison and bind them into a single data frame
          df <- do.call(rbind, lapply(input$comp_all, function(comp) {
            # Get results object df
            res_df <- app_object$res[[comp]]
            # Add a 'comparison' column if it doesn't exist
            res_df$comparison <- comp
            # Combine results dfs
            res_df[, c(shared_columns, 'comparison')]
          }))

        } else {
          df <- app_object$res[[ input$comp_all ]]
          df$comparison <- input$comp_all
        }


        fc.thres <- ifelse(input$fc.thres == '' | is.na(input$fc.thres),
                           config$ui$de_analysis$filters$log2fc_threshold, input$fc.thres)
        fdr.thres <- ifelse(input$fdr.thres == '' | is.na(input$fdr.thres),
                            config$ui$de_analysis$filters$fdr_threshold, input$fdr.thres)

        # filter df
        if(input$toggle_filters){
            idx <- df$padj < fdr.thres & !is.na(df$padj) & abs(df$log2FoldChange) >= fc.thres
            if(sum(idx, na.rm=TRUE) == 0){
              res_data$tbl <- NULL
              validate(
                  need(sum(idx, na.rm=TRUE) > 0, 'No DE genes found!')
              )
            }
            df <- df[idx,]
        }

        sidx <- grep('symbol', colnames(df), ignore.case=TRUE)
        colnames(df)[sidx] <- tolower(colnames(df)[sidx])

        if(!'gene' %in% colnames(df))
          df$gene <- sub('\\.\\d+$','',rownames(df))

        gene.col <- 'gene'
        symbol.col <- ifelse('symbol' %in% colnames(df), 'symbol', 'SYMBOL')

        cols.to.drop <- config$server$de_analysis$de_table$cols.to.drop
        df <- df %>% as.data.frame %>%
          relocate(!!symbol.col, .before=.data$baseMean) %>%
          relocate(!!gene.col, .before=!!symbol.col) %>%
          select(-any_of(c(cols.to.drop, toupper(cols.to.drop))))

        res_data$tbl <- df
      }
    }) # observeEvent update res_data

    # Update res_data$tbl based on processed data
    observe({
      res_data$tbl <- processed_table()
    })

    # Observer for Select all checkbox
    observeEvent(input$select_all_comp, {
      # Get all possible comparison options
      all_comps <- names(app_object$res)
      # Update the select_none checkbox
      updateSelectizeInput(session, 'comp_all', selected=all_comps)
    })

    # Observer for Select none checkbox
    observeEvent(input$select_none_comp, {
      # Update comp_all with no selected comparisons
      updateSelectizeInput(session, 'comp_all', selected=character(0))
    })

    # render DE table
    output$detable <- renderDT({
      validate(
        need(input$comp_all != '', 'Waiting for selection')
      )

      validate(
        need(!is.null(res_data$tbl), 'No DE genes to display! Loosen DE criteria or uncheck "Only DE genes"')
      )

      tbl <- res_data$tbl

      tbl %>%
        datatable(rownames=FALSE,
                  selection=list(mode='multiple')) %>%
        formatSignif(columns=config$server$de_analysis$de_table$format_significant$columns,
                     digits=config$server$de_analysis$de_table$format_significant$digits)
    }) # renderDT

    detable_proxy <- dataTableProxy('detable')

    observeEvent(input$reset_detable, {
      detable_proxy %>% selectRows(NULL)
    })

    observeEvent(input$add_selected, {
      tbl <- res_data$tbl
      sel <- input$detable_rows_selected

      # handle NAs in symbol
      s <- tbl$symbol
      s[is.na(s)] <- tbl$gene[is.na(s)]

      if(is.null(sel)){
        showNotification(
          'Cannot add genes, no rows selected', type='warning'
        )
        validate(
          need(!is.null(sel), '')
        )
      } else if(all(s[sel] %in% input$gene.to.plot)){
        showNotification(
          'Selected genes already present in scratchpad, skipping', type='warning'
        )
      } else {
        if(is.null(input$gene.to.plot)) selected <- s[sel]
        else selected <- unique(c(input$gene.to.plot, s[sel]))

        new_genes <- setdiff(selected, input$gene.to.plot)
        showNotification(
          paste('Adding', length(new_genes),
                'new genes to scratchpad')
        )

        # update gene selector with clicked genes
        updateSelectizeInput(session, 'gene.to.plot',
                             choices=gene.id$gene,
                             selected=selected,
                             server=TRUE)
      }
    }) # observeEvent

    ####################### Gene plot ##########################

    gene_plot_args <- reactive({
      list(gene.to.plot=gene_scratchpad(),
           gene.id=gene.id$gene,
           comp_all=input$comp_all
      )
    })

    gene_scratchpad <- reactive({
      if(!is.null(input$gene.to.plot)){
        input$gene.to.plot
      } else {
        ''
      }
    })

    genePlotServer('gene_plot', app_object,
                   coldata=coldata.all,
                   gene_plot_args)

    ####################### PCA Plot ###########################

    pcaPlotServer('pcaplot', app_object, coldata.all)

    ####################### MA plot #############################

    maplot_args <- reactive({
      list(fdr.thres=input$fdr.thres,
           fc.thres=input$fc.thres,
           gene.to.plot=gene_scratchpad())
    })

    maPlotServer('maplot', app_object, maplot_args)

    ####################### Scatter plot #############################

    scatterplot_args <- reactive({
      list(fdr.thres=input$fdr.thres,
           fc.thres=input$fc.thres,
           gene.to.plot=gene_scratchpad())
    })

    scatterPlotServer('scatterplot', app_object, scatterplot_args)

    ##################### UpSet plot #########################

    upset_plot_args <- reactive({
      list(fdr.thres=input$fdr.thres,
           fc.thres=input$fc.thres)
    })

    upset_data <- upsetPlotServer('upset_plot',
                                  app_object,
                                  upset_plot_args,
                                  gene_scratchpad,
                                  reactive({ input$reset.genes }))

    upset_table <- reactiveValues(tbl=NULL, intersections=NULL, set_labels=NULL)

    observeEvent(upset_data(), {
      g <- upset_data()$genes

      # only update scratchpad if different genes returned
      if(length(setdiff(g, input$gene.to.plot)) != 0){
        # update gene selector with clicked genes
        updateSelectizeInput(session, 'gene.to.plot',
                             choices=gene.id$gene,
                             selected=g,
                             server=TRUE)
      }

      # save the set labels
      upset_table$tbl <- upset_data()[[ 'tbl' ]][[ 'all' ]]
      labels <- upset_data()[[ 'tbl' ]][[ 'set_labels' ]]
      upset_table$set_labels <- labels

      # get gene symbols in upset intersections
      intersections <- lapply(unname(labels), function(x){
                         idx <- upset_table$tbl[, 'set'] == x
                         upset_table$tbl[idx, 'symbol']
                       })
      names(intersections) <- unname(labels)
      upset_table$intersections <- intersections

    })

    ######################## Heatmap #############################

    hmap_plot_args <- reactive({
      list(fdr.thres=input$fdr.thres,
           fc.thres=input$fc.thres,
           upset_data=list(genes=upset_table$intersections,
                           labels=upset_table$set_labels))
    })

    heatmapServer('heatmap', app_object, coldata.all, hmap_plot_args, gene_scratchpad)

    #################### Functional enrichment #######################

    enrich_data <- enrichServer('func_enrich', app_object,
                                upset_table,
                                gene_scratchpad,
                                reactive({ input$reset.genes }))

    observeEvent(enrich_data(), {
      g <- enrich_data()$genes

      # only update scratchpad if different genes returned
      if(length(setdiff(g, input$gene.to.plot)) != 0){
        # update gene selector with clicked genes
        updateSelectizeInput(session, 'gene.to.plot',
                             choices=gene.id$gene,
                             selected=g,
                             server=TRUE)
      }
    })


    ######################## DEG patterns ########################

    pattern_plot_args <- reactive({
      list(
        gene_scratchpad=gene_scratchpad(),
        upset_data=list(genes=upset_table$intersections,
                        labels=upset_table$set_labels)
      )
    })

    patternPlotServer('deg_plot', app_object, coldata.all,
                      pattern_plot_args)

    ######################### Help buttons #######################

    # general help
    helpButtonServer(id='gene_scratchpad_help', size='l')
    helpButtonServer('de_filters_help')
    helpButtonServer('metadata_sidebar_help')
    helpButtonServer('de_cmp_help')

    # de analysis help
    helpButtonServer('de_summary_help', size='l')
    helpButtonServer('de_meta_help', size='l')
    helpButtonServer('de_tbl_help', size='l')

  }

  shinyApp(ui, server, ...)

}
