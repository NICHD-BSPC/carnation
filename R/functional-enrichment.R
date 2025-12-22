#' Functional enrichment module
#'
#' @description
#' UI & module to show functional enrichment tables & plots.
#'
#' @param id ID string used to match the ID used to call the module server function
#' @param panel string, can be 'sidebar' or 'main'
#' @param tab string, if 'table' show table settings, if 'plots' show plot settings;
#' if 'compare_results', show comparison settings.
#' @param id ID string used to match the ID used to call the module UI function
#' @param obj reactiveValues object containing carnation object
#' @param upset_table reactive, data from upset plot module
#' @param gene_scratchpad reactive, genes selected in gene scratchpad
#' @param reset_genes reactive to reset genes in scratchpad
#' @param config reactive list with config settings
#'
#' @returns
#' UI returns tagList with plot UI
#' server returns reactive with gene selected from functional enrichment tables.
#'
#' @rdname funenrichmod
#' @name funenrichmod
NULL

#' @rdname funenrichmod
#' @export
enrichUI <- function(id, panel, tab='none'){
  ns <- NS(id)

  # get default config
  config <- get_config()

  if(panel == 'sidebar'){
    if(tab == 'table'){
      tag <-
        tagList(
          fluidRow(
            column(6, align='left',
              tags$label(class='control-label',
                         'Comparison')
            ),
            column(6, align='right',
              helpButtonUI(ns('tbl_ctrls_help'))
            ) # column
          ), # fluidRow

          selectizeInput(ns('comp_tbl'),
                         label=NULL,
                         choices=NULL,
                         selected=NULL
          ), # selectizeInput

          fluidRow(
            column(5, h5('Direction')
           ), # fluidRow
            column(7,
              selectInput(ns('geneset_tbl'), label=NULL,
                       choices=NULL,
                       selected=NULL
              ) # selectInput
            ) # column
          ), # fluidRow

          fluidRow(
            column(5, h5('Database')), # column
            column(7,
              selectInput(ns('pathway_tbl'), label=NULL,
                       choices=NULL,
                       selected=NULL
              ) # selectInput
            ) # column
          ), # fluidRow

          bsCollapse(
            bsCollapsePanel(span(icon('search'), 'Search table'),
              value='search',
              fluidRow(
                column(5, h5('Search in')),
                column(7,
                  selectInput(ns('search_opts'), label=NULL,
                    choices=c('genes + description',
                              'description',
                              'genes')
                  ) # selectInput
                ) # column
              ), # fluidRow

              fluidRow(
                column(12,
                  selectizeInput(ns('search_txt'), label=NULL,
                    choices='', selected='', multiple=TRUE,
                    options=list(create=TRUE, placeholder='Enter text to search')
                  ) # selectizeInput
                ) # column
              ), # fluidRow

              fluidRow(align='center',
                actionButton(ns('search_func_do'),
                             label='Refresh',
                             icon=icon('arrows-rotate'),
                             class='btn-primary',
                             style='margin-bottom: 5px;')
              ) # fluidRow
            ) # bsCollapsePanel
          ), # bsCollapse

          bsCollapse(
            bsCollapsePanel(span(icon('filter'), 'Subset table'),
              value='subset',
              fluidRow(
                column(5, h5('Subset by')),
                column(7,
                  selectInput(ns('subset_opts'), label=NULL,
                    choices=c('none', 'gene_scratchpad', 'upset_intersections')
                  ) # selectInput
                ) # column
              ), # fluidRow

              conditionalPanel(paste0('input["', ns('subset_opts'), '"] == "upset_intersections"'),
                fluidRow(
                  column(12,
                    selectizeInput(ns("upset_intersect"), label=NULL,
                                   choices=NULL, selected=NULL
                    ) # selectInput
                  ) # column
                ) # fluidRow

              ), # conditionalPanel

              fluidRow(align='center',
                actionButton(ns('subset_func_do'),
                             label='Refresh',
                             icon=icon('arrows-rotate'),
                             class='btn-primary',
                             style='margin-bottom: 5px;')
              ) # fluidRow
            ) # bsCollapsePanel
          ), # bsCollapse

          fluidRow(
            column(12, strong('More options'))
          ),

          fluidRow(
            column(5, h5('Genes per line')),
            column(7,
              numericInput(ns('genes.per.line'), label=NULL,
                value=config$ui$functional_enrichment$table$genes.per.line,
                min=1, step=1
              ) # numericInput
            ) # column
          ), # fluidRow


          conditionalPanel(paste0('input["', ns('functable'), '"] == "Distill enrichment"'),

            fluidRow(
              column(5, h5('# of terms')),
              column(7,
                numericInput(ns('numcat.distill.tbl'), label=NULL,
                  value=config$ui$functional_enrichment$plots$emap_distill$numcat
                ) # numericInput
              ) # column
            ) # fluidRow

          ), # conditionalPanel

          conditionalPanel(paste0('input["', ns('functable'), '"] == "Fuzzy clustering"'),
            fluidRow(
              column(5, h5('# of terms')),
              column(7,
                numericInput(ns('numcat.fuzzy.tbl'), label=NULL,
                  value=config$ui$functional_enrichment$plots$emap_fuzzy$numcat
                ) # numericInput
              ) # column
            ) # fluidRow

          )

        ) # tagList
      } else if(tab == 'plots'){
        tag <-
          tagList(
            fluidRow(
              column(6, align='left',
                tags$label(class='control-label',
                           'Comparison')
              ),
              column(6, align='right',
                helpButtonUI(ns('func_plt_help'))
              ) # column
            ), # fluidRow

            selectizeInput(ns('comp_fun'),
                           label=NULL,
                           choices=NULL,
                           selected=NULL
            ), # selectizeInput

            fluidRow(
              column(5, h5('Direction')
             ), # fluidRow
              column(7,
                selectInput(ns('geneset'), label=NULL,
                         choices=NULL,
                         selected=NULL
                ) # selectInput
              ) # column
            ), # fluidRow

            fluidRow(
              column(5, h5('Database')), # column
              column(7,
                selectInput(ns('pathway'), label=NULL,
                         choices=NULL,
                         selected=NULL
                ) # selectInput
              ) # column
            ), # fluidRow

            wellPanel(
              strong('Plot options'),
              uiOutput(ns('plot_opts'))
            ) # wellPanel

          ) # tagList
    } else if(tab == 'compare_results'){
      tag <-
        tagList(

          fluidRow(
            column(6, align='left',
              tags$label(class='control-label',
                         'Comparison 1')
            ),
            column(6, align='right',
              helpButtonUI(ns('func_cmp_help'))
            ) # column
          ), # fluidRow

          selectizeInput(ns('comp_fun1'),
                         label=NULL,
                         choices=NULL,
                         selected=NULL
          ), # selectizeInput

          fluidRow(
            column(5, h5('Direction')),
            column(7,
              selectInput(ns('geneset1'), label=NULL,
                          choices=NULL,
                          selected=NULL
              ) # selectInput
            ) # column
          ), # fluidRow

          selectizeInput(ns('comp_fun2'),
                         label='Comparison 2',
                         choices=NULL,
                         selected=NULL
          ), # selectizeInput

          fluidRow(
            column(5, h5('Direction')),
            column(7,
              selectInput(ns('geneset2'), label=NULL,
                          choices=NULL,
                          selected=NULL
              ) # selectInput
            ) # column
          ), # fluidRow

          fluidRow(
            column(5,
              tags$label(class='control-label', 'Database')
            ), # column
            column(7,
              selectInput(ns('pathway1'), label=NULL,
                          choices=NULL,
                          selected=NULL
              ) # selectInput
            ) # column
          ), # fluidRow

          fluidRow(
            column(6, align='right', style='margin-bottom: 10px;',
              actionButton(ns('swap_comp'), label='Swap comparisons',
                           icon=icon('arrows-rotate'),
                           status='info')
            ), # column

            column(6, align='left',
              actionButton(ns('swap_geneset'), label='Swap directions',
                           icon=icon('arrows-rotate'),
                           status='info')
            )
          ),

          wellPanel(
            strong('Plot options'),
            uiOutput(ns('cmp_plot_opts'))
          ) # wellPanel

        ) # tagList
    } # tab
  } else if(panel == 'main'){
    if(tab == 'table'){
      tag <-
      tagList(
        tabsetPanel(id=ns('functable'),

          tabPanel('Enrichment',

            fluidRow(
              column(1, align='left',
               helpButtonUI(ns('func_tbl_help'))
              ), # column
              column(11,
                uiOutput(ns('func_table_title'))
              ) # column
            ), # fluidRow

            fluidRow(
              column(12,
                withSpinner(
                  DTOutput(ns('func_table'))
                ) # withSpinner
              ) # column
            ), # fluidRow

            fluidRow(align='center',
              br(),
              column(12,
                splitLayout(
                  cellWidths=c('15%', '25%', '20%'),

                  strong('Selection options'),
                  actionButton(ns('add_func_selected'), 'Add to scratchpad'),
                  actionButton(ns('reset_functable'),
                               'Reset selection',
                               class='btn-primary')
                ) # splitLayout
              ) # column
            ) # fluidRow

          ), # tabPanel

          tabPanel('Distill enrichment',
            fluidRow(
              column(1, align='left',
               helpButtonUI(ns('distill_tbl_help'))
              ), # column
              column(11,
                uiOutput(ns('distill_table_title'))
              ) # column
            ), # fluidRow

            fluidRow(
              column(12,
                withSpinner(
                  DTOutput(ns('distill_table'))
                ) # withSpinner
              ) # column
            ), # fluidRow

            fluidRow(align='center',
              br(),
              column(12,
                splitLayout(
                  cellWidths=c('15%', '25%', '20%'),

                  strong('Selection options'),
                  actionButton(ns('add_distill_selected'), 'Add to scratchpad'),
                  actionButton(ns('reset_distilltable'),
                               'Reset selection',
                               class='btn-primary')
                ) # splitLayout
              ) # column
            ) # fluidRow

          ),

          tabPanel('Fuzzy clustering',
            fluidRow(
              column(1, align='left',
                helpButtonUI(ns('fuzzy_tbl_help'))
              ), # column
              column(11,
                uiOutput(ns('fuzzy_table_title'))
              ) # column
            ), # fluidRow

            fluidRow(
              column(12,
                withSpinner(
                  DTOutput(ns('fuzzy_table'))
                ) # withSpinner
              ) # column
            ), # fluidRow

            fluidRow(align='center',
              br(),
              column(12,
                splitLayout(
                  cellWidths=c('15%', '25%', '20%'),

                  strong('Selection options'),
                  actionButton(ns('add_fuzzy_selected'), 'Add to scratchpad'),
                  actionButton(ns('reset_fuzzytable'),
                               'Reset selection',
                               class='btn-primary')
                ) # splitLayout
              ) # column
            ) # fluidRow

          ) # tabPanel
        ) # tabsetPanel

      ) # tagList
    } else if(tab == 'plots'){
      tag <-
        tagList(
          fluidRow(style='margin-top: 10px;',
            column(2, align='right', style='margin-top: 5px;',
                   strong('Type of plot')),
            column(2, align='left',
              selectInput(ns('plottype'),
                label=NULL,
                choices=config$ui$functional_enrichment$plottype$choices,
                selected=config$ui$functional_enrichment$plottype$default
              ) # selectInput
            ), # column
            column(6, align='left',
              uiOutput(ns('enrichplot_title'))
            ),
            column(2,
              # enrichplot dload & help btns
              uiOutput(ns('enrichplot_btns'))
            )
          ),
          uiOutput(ns('enrichplot'))
        ) # tagList
    } else if(tab == 'compare_results'){
      tag <-
        tagList(
          fluidRow(style='margin-top: 10px;',
            column(2, align='right', style='margin-top: 5px;',
                   strong('Type of plot')),
            column(2, align='left',
              selectInput(ns('comptype'), label=NULL,
                          choices=c('summary_overview','radar','horizon')
              ) # selectInput
            ), # column
            column(6, align='left',
              uiOutput(ns('compenrich_title'))
            ),
            column(2,
              # enrichplot dload & help btns
              uiOutput(ns('compplot_btns'))
            )
          ),
          uiOutput(ns('compare_enrich'))
        ) # tagList
    } # tab
  } # panel
  tag
}


#' @rdname funenrichmod
#' @export
enrichServer <- function(id, obj, upset_table,
                         gene_scratchpad, reset_genes, config){

  moduleServer(
    id,

    function(input, output, session){

      ns <- NS(id)

      app_object <- reactive({
        list(res=obj$res,
             enrich=obj$enrich,
             genetonic=obj$genetonic)
      })

      enrich_data <- reactiveValues(enrich_list=NULL,
                                    genetonic_list=NULL,
                                    distill=NULL,
                                    fuzzy=NULL)

      flags <- reactiveValues(data_loaded=0)

      upset_data <- reactive({
        list(intersections=upset_table$intersections,
             labels=upset_table$set_labels)
      })

      # reactive values to track genes clicked/labeled
      genes_clicked <- reactiveValues(g=NULL)

      observeEvent(config(), {
        updateNumericInput(session, 'genes.per.line',
                           value=config()$ui$functional_enrichment$table$genes.per.line)

        updateNumericInput(session, 'numcat.distill.tbl',
                           value=config()$ui$functional_enrichment$plots$emap_distill$numcat)

        updateNumericInput(session, 'numcat.fuzzy.tbl',
                           value=config()$ui$functional_enrichment$plots$emap_fuzzy$numcat)

        updateSelectInput(session, 'plottype',
                          choices=config()$ui$functional_enrichment$plottype$choices,
                          selected=config()$ui$functional_enrichment$plottype$default)

      })

      observeEvent(app_object(), {
        enrich_data$enrich_list <- NULL
        enrich_data$genetonic_list <- NULL
        enrich_data$distill <- NULL
        enrich_data$fuzzy <- NULL


        if(is.null(names(app_object()$enrich))){
          updateSelectizeInput(session,
                               'comp_tbl',
                               choices=c(''),
                               selected=c(''))
          updateSelectizeInput(session,
                               'comp_fun',
                               choices=c(''),
                               selected=c(''))
          updateSelectizeInput(session,
                               'comp_fun1',
                               choices=c(''),
                               selected=c(''))
          updateSelectizeInput(session,
                               'comp_fun2',
                               choices=c(''),
                               selected=c(''))
        } else {
          updateSelectizeInput(session,
                               'comp_tbl',
                               choices=names(app_object()$enrich),
                               selected=names(app_object()$enrich)[1])
          updateSelectizeInput(session,
                               'comp_fun',
                               choices=names(app_object()$enrich),
                               selected=names(app_object()$enrich)[1])
          updateSelectizeInput(session,
                               'comp_fun1',
                               choices=names(app_object()$enrich),
                               selected=names(app_object()$enrich)[1])
          updateSelectizeInput(session,
                               'comp_fun2',
                               choices=names(app_object()$enrich),
                               selected=names(app_object()$enrich)[2])
        }

        # reset genes when new data loaded
        genes_clicked$g <- NULL

        flags$data_loaded <- flags$data_loaded + 1
      })

      # observer to update genes_clicked
      observeEvent(gene_scratchpad(), {

        if(any(gene_scratchpad() != ''))
          genes_clicked$g <- gene_scratchpad()

      })

      observeEvent(reset_genes(), {
        genes_clicked$g <- NULL
      })

      # observer to update upset intersections
      observeEvent(upset_data(), {
        updateSelectizeInput(session, 'upset_intersect',
                             choices=upset_data()$labels,
                             server=TRUE)
      })

      # observer to update comparison menu
      observeEvent(input$comp_fun, {
        validate(
          need(input$comp_fun != '', 'No comparisons selected')
        )

        validate(
          need(!is.null(app_object()$res), 'Waiting for selection')
        )

        validate(
          need(!is.null(app_object()$enrich), 'Functional enrichment results not available')
        )

        validate(
          need(input$comp_fun %in% names(app_object()$enrich),
               'Functional enrichment results not available for this comparison')
        )

        elements <- names(app_object()$enrich[[input$comp_fun]])

        # skip 'res' element
        elements <- setdiff(elements, 'res')

        if(!input$geneset %in% elements){
            # choose pathway results to show
            updateSelectInput(session,
                              'geneset',
                              choices=elements,
                              selected=elements[1])
        } else {
            # choose pathway results to show
            updateSelectInput(session,
                              'geneset',
                              choices=elements,
                              selected=input$geneset)
        }
      }) # observeEvent

      # observer to update comparison menu
      observeEvent(input$comp_tbl, {
        validate(
          need(input$comp_tbl != '', 'No comparisons selected')
        )

        validate(
          need(!is.null(app_object()$res), 'Waiting for selection')
        )

        validate(
          need(!is.null(app_object()$enrich), 'Functional enrichment results not available')
        )

        validate(
          need(input$comp_tbl %in% names(app_object()$enrich),
               'Functional enrichment results not available for this comparison')
        )

        elements <- names(app_object()$enrich[[input$comp_tbl]])

        # skip 'res' element
        elements <- setdiff(elements, 'res')

        if(!input$geneset_tbl %in% elements){
            # choose pathway results to show
            updateSelectInput(session,
                              'geneset_tbl',
                              choices=elements,
                              selected=elements[1])
        } else {
            # choose pathway results to show
            updateSelectInput(session,
                              'geneset_tbl',
                              choices=elements,
                              selected=input$geneset_tbl)
        }
      }) # observeEvent

      # observer to update pathway menu when comparison/geneset
      # is selected
      observeEvent(c(input$comp_fun, input$geneset), {
        # get names of available pathway analyses
        # named vector of pathway names
        p <- unlist(config()$server$functional_enrichment$pathways)

        all.p.names <- names(app_object()$enrich[[input$comp_fun]][[input$geneset]])
        p.choices <- p[p %in% all.p.names]
        # if database names are not in known options
        if(length(p.choices) < length(all.p.names)){
            p.choices <- c(p.choices, all.p.names[!all.p.names %in% p])
        }

        if(!input$pathway %in% p.choices){
            updateSelectInput(session,
                              'pathway',
                              choices=p.choices,
                              selected=p.choices[1])
        } else {
            updateSelectInput(session,
                              'pathway',
                              choices=p.choices,
                              selected=input$pathway)
        }

      }) # observeEvent

      # observer to update pathway menu when comparison/geneset
      # is selected
      observeEvent(c(input$comp_tbl, input$geneset_tbl), {
        # get names of available pathway analyses
        # named vector of pathway names
        p <- unlist(config()$server$functional_enrichment$pathways)

        all.p.names <- names(app_object()$enrich[[input$comp_tbl]][[input$geneset_tbl]])
        p.choices <- p[p %in% all.p.names]
        # if database names are not in known options
        if(length(p.choices) < length(all.p.names)){
            p.choices <- c(p.choices, all.p.names[!all.p.names %in% p])
        }

        if(!input$pathway_tbl %in% p.choices){
            updateSelectInput(session,
                              'pathway_tbl',
                              choices=p.choices,
                              selected=p.choices[1])
        } else {
            updateSelectInput(session,
                              'pathway_tbl',
                              choices=p.choices,
                              selected=input$pathway_tbl)
        }

      }) # observeEvent

      observeEvent(c(input$comp_fun, input$geneset, input$pathway, flags$data_loaded), {
        validate(
          need(!is.null(app_object()$enrich), '')
        )
        df <- app_object()$enrich[[input$comp_fun]][[input$geneset]][[input$pathway]]

        # get genetonic df
        if('res' %in% names(app_object()$enrich[[input$comp_fun]])){
            res <- app_object()$res[[ app_object()$enrich[[ input$comp_fun ]][['res']] ]]
        } else if(input$comp_fun %in% names(app_object()$enrich)){
            res <- app_object()$res[[input$comp_fun]]
        } else {
            res <- NULL
        }

        genetonic.list <- app_object()$genetonic[[input$comp_fun]][[input$geneset]][[input$pathway]]
        title <- paste0(input$comp_fun, ' | ',
                        input$geneset, ' | ',
                        input$pathway)
        enrich_data$enrich_list <- list(df=df, res=res, title=title)
        enrich_data$genetonic_list <- genetonic.list
      })

      observeEvent(c(input$comp_tbl, input$geneset_tbl, input$pathway_tbl), {
        df <- app_object()$enrich[[input$comp_tbl]][[input$geneset_tbl]][[input$pathway_tbl]]

        # get genetonic df
        if('res' %in% names(app_object()$enrich[[input$comp_tbl]])){
            res <- app_object()$res[[ app_object()$enrich[[ input$comp_tbl ]][['res']] ]]
        } else if(input$comp_tbl %in% names(app_object()$enrich)){
            res <- app_object()$res[[input$comp_tbl]]
        } else {
            res <- NULL
        }

        genetonic.list <- app_object()$genetonic[[input$comp_tbl]][[input$geneset_tbl]][[input$pathway_tbl]]
        title <- paste0(input$comp_tbl, ' | ',
                        input$geneset_tbl, ' | ',
                        input$pathway_tbl)
        enrich_data$enrich_list <- list(df=df, res=res, title=title)
        enrich_data$genetonic_list <- genetonic.list
      })

      ################### Functional enrichment table #########################

      # function to find & return case-insensitive matches in a column
      # with the matches in bold (optional)
      # - can do exact or partial matching
      match_text_tbl <- function(df, search_txt, gene_col, gene_sep, highlight, exact_match){
        ll <- lapply(seq_len(nrow(df)), function(x){
                tmp <- strsplit(df[x, gene_col], gene_sep)[[1]]

                # case-insensitive search
                if(exact_match){
                  gg <- which(tolower(tmp) %in% tolower(search_txt))
                } else {
                  gg <- unique(unlist(lapply(tolower(search_txt), function(x) grep(x, tolower(tmp)))))
                }

                if(length(gg) > 0){
                  if(highlight) tmp[ gg ] <- paste0('<b>', tmp[ gg ], '</b>')
                  return(list(idx=x, txt=paste(tmp, collapse=gene_sep)))
                }
              })

        return(
          list(
            idx = unlist(lapply(ll, function(x) x$idx)),
            txt = unlist(lapply(ll, function(x) x$txt))
          )
        )
      }

      # function to subset FE table
      #
      # - case-insensitive search with exact or partial matching
      #   for genes or description
      # - can subset by gene scratchpad genes or upset intersections
      # - matches can be optionally highlighted in the table view
      #
      subset_func_tbl <- function(df, gene_col, desc_col, gene_sep,
                                  search_opts, search_txt,
                                  subset_opts, upset_intersect,
                                  highlight=TRUE, quiet=TRUE){
        if(!is.null(search_txt)){
          if(search_opts == 'genes + description'){
            # search in description
            # turn off highlighting for description
            ll1 <- match_text_tbl(df, search_txt, desc_col, gene_sep, highlight, exact_match=FALSE)

            # search in genes
            ll2 <- match_text_tbl(df, search_txt, gene_col, gene_sep, highlight, exact_match=FALSE)

            # keep union of search hits
            keep_rows <- unique(c(ll1$idx, ll2$idx))
            df[ll1$idx, desc_col] <- ll1$txt
            df[ll2$idx, gene_col] <- ll2$txt
          } else if(search_opts == 'description'){

            # turn off highlighting for description
            ll <- match_text_tbl(df, search_txt, desc_col, gene_sep, highlight, exact_match=FALSE)

            keep_rows <- ll$idx
            df[ll$idx, desc_col] <- ll$txt
          } else if(search_opts == 'genes'){
            ll <- match_text_tbl(df, search_txt, gene_col, gene_sep, highlight, exact_match=FALSE)

            keep_rows <- ll$idx
            df[ll$idx, gene_col] <- ll$txt
          }

          if(length(keep_rows) == 0){
            if(!quiet){
              showNotification(
                'Search text not found in table! Not subsetting', type='warning'
              )
            }
          } else {
            if(!quiet){
              showNotification(
                paste('Subsetting table to keep', length(keep_rows), 'rows')
              )
            }

            keep_rows <- keep_rows[order(keep_rows)]
            df <- df[keep_rows,]
          }
        }

        # subset table
        if(subset_opts != 'none'){
          if(subset_opts == 'gene_scratchpad'){
            g <- gene_scratchpad()
            if(length(g) == 0 || all(g == '')){
              if(!quiet){
                showNotification(
                  'No genes selected in scratchpad! Not subsetting',
                  type='warning'
                )
              }
              g <- NULL
            }
          } else if(subset_opts == 'upset_intersections'){
            if(upset_intersect == ''){
              if(!quiet){
                showNotification(
                  'Select an Upset intersection to subset with',
                  type='warning'
                )
              }
              g <- NULL
            } else {
              g <- upset_data()$intersections[[ upset_intersect ]]
            }
          }

          if(!is.null(g)){
            ll <- match_text_tbl(df, g, gene_col, gene_sep, highlight, exact_match=TRUE)

            if(length(ll$idx) > 0){
              df <- df[ll$idx, ]
              df[, gene_col] <- ll$txt
            } else {
              if(!quiet){
                showNotification(
                  'Selected gene not found in table! Not subsetting',
                  type='warning'
                )
              }
            }

          }
        }
        df
      } # function

      # observer to get functional enrichment tbl
      get_func_table <- eventReactive(c(enrich_data$enrich_list,
                                        input$search_func_do,
                                        input$subset_func_do), {
        validate(
          need(!is.null(app_object()$enrich), 'Functional enrichment results not available')
        )

        validate(
          need(input$comp_tbl != '', 'Waiting for selection')
        )

        validate(
          need(input$comp_tbl %in% names(app_object()$enrich),
               'Functional enrichment results not available for this comparison')
        )

        df <- as.data.frame(enrich_data$enrich_list$df)

        validate(
          need(nrow(df) > 0, 'No DE genes found!')
        )

        # support gseaResult and enrichResult
        gene.id.col <- if('geneID' %in% colnames(df)) 'geneID'
          else if('core_enrichment' %in% colnames(df)) 'core_enrichment'

        # replace '/' with ',' in gene id column
        # - this allows us to add <b> </b> tags to highlight genes
        df[[gene.id.col]] <- gsub('/', ',', df[[ gene.id.col ]])

        df <- subset_func_tbl(df, gene_col=gene.id.col, desc_col='Description', gene_sep=',',
                              input$search_opts, input$search_txt,
                              input$subset_opts, input$upset_intersect)

        df
      }) # eventReactive get_func_table

      output$func_table_title <- renderUI({
        fluidRow(align='center',
          h4(enrich_data$enrich_list$title)
        )
      })


      output$func_table <- renderDT({
        df <- get_func_table()

        float_idx <- vapply(df, function(x) typeof(x) %in% c('double', 'float'), logical(1))
        format_cols <- colnames(df)[float_idx]
        digits <- config()$server$functional_enrichment$table$enrichment$format_significant$digits

        # support gseaResult and enrichResult
        gene.id.col <- if('geneID' %in% colnames(df)) 'geneID'
          else if('core_enrichment' %in% colnames(df)) 'core_enrichment'

        df[, gene.id.col] <- format_genes(df[, gene.id.col],
                                          genes.per.line=input$genes.per.line, sep=',')

        if('Count' %in% colnames(df)){
          rem.cols <- setdiff(colnames(df), c('Count',
                                              gene.id.col))
          desc.idx <- which(rem.cols == 'Description')
          col.order <- c(rem.cols[seq_len(desc.idx)], gene.id.col,
                         'Count', rem.cols[(desc.idx+1):length(rem.cols)])
        } else {
          rem.cols <- setdiff(colnames(df), gene.id.col)
          desc.idx <- which(rem.cols == 'Description')
          col.order <- c(rem.cols[seq_len(desc.idx)], gene.id.col,
                         rem.cols[(desc.idx+1): length(rem.cols)])
        }
        cols_to_format <- intersect(colnames(df), format_cols)

        df %>%
          datatable(rownames=FALSE,
                    escape=FALSE,
                    options=list(
                      autoWidth=TRUE,
                      columnDefs=list(list(width='40%', targets=c(3)))
                    )) %>%
          formatSignif(columns=cols_to_format,
                       digits=digits) %>%
          formatStyle('Description', 'white-space'='nowrap')
      })

      # selection from FE table
      functable_proxy <- dataTableProxy('func_table')

      observeEvent(input$reset_functable, {
        functable_proxy %>% selectRows(NULL)
      })

      observeEvent(input$add_func_selected, {
        tbl <- get_func_table()
        sel <- input$func_table_rows_selected

        # handle NAs in symbol
        if('geneID' %in% colnames(tbl)){
          s <- tbl$geneID
        } else if('core_enrichment' %in% colnames(tbl)){
          s <- tbl$core_enrichment
        }

        if(is.null(sel)){
          showNotification(
            'Cannot add genes, no rows selected', type='warning'
          )
          validate(
            need(!is.null(sel), '')
          )
        } else {
          sel_genes <- unique(unlist(lapply(s[sel], function(x) strsplit(x, ',')[[1]])))
          if(all(sel_genes %in% genes_clicked$g)){
            showNotification(
              'Selected genes already present in scratchpad, skipping', type='warning'
            )
          } else {
            if(is.null(genes_clicked$g)) selected <- sel_genes
            else selected <- unique(c(genes_clicked$g, sel_genes))

            new_genes <- setdiff(selected, genes_clicked$g)
            showNotification(
              paste('Adding', length(new_genes),
                    'new genes to scratchpad')
            )

            genes_clicked$g <- selected
          }
        }
      })

      # this observer looks only at the controls shown on the Plots tab
      observeEvent(c(input$comp_fun,
                     input$geneset,
                     input$pathway,
                     flags$data_loaded,
                     input$search_func_do,
                     input$subset_func_do), {

        validate(
          need(!is.null(app_object()$res), 'Waiting for selection')
        )

        validate(
          need(!is.null(app_object()$enrich), 'Functional enrichment results not available')
        )

        validate(
          need(input$comp_fun != '' & input$geneset != '',
               'Waiting for selection')
        )

        isolate({
          updateSelectInput(session, 'comp_tbl',
                            selected=input$comp_fun)
          updateSelectInput(session, 'geneset_tbl',
                            selected=input$geneset)
          updateSelectInput(session, 'pathway_tbl',
                            selected=input$pathway)
        })

        if(is.null(input$numcat.distill.tbl))
          numcat <- config()$ui$functional_enrichment$plots$enrichment_map$numcat
        else numcat <- input$numcat.distill.tbl

        validate(
          need(numcat > 0, 'Number of terms must be > 0')
        )

        # get obj to visualize
        l <- enrich_data$enrich_list$df

        validate(
          need(nrow(as.data.frame(l)) > 0, 'No DE genes found!')
        )

        res <- enrich_data$enrich_list$res

        validate(
            need(!is.null(res), 'DE results not found for this geneset')
        )

        genetonic.list <- enrich_data$genetonic_list
        l_gs <- genetonic.list$l_gs
        anno_df <- genetonic.list$anno_df

        # subset table
        l_gs <- subset_func_tbl(l_gs, gene_col='gs_genes', desc_col='gs_description', gene_sep=',',
                                input$search_opts, input$search_txt,
                                input$subset_opts, input$upset_intersect)

        n_gs <- min(nrow(l_gs), numcat)

        df <- tryCatch(
                distill_enrichment(l_gs,
                                 res,
                                 anno_df,
                                 n_gs = n_gs,
                                 cluster_fun = "cluster_markov"),
                warning = function(w){ w },
                error = function(e){ e }
              )

        if(inherits(df, 'error')){
           showNotification(
             paste('Distill enrichment error:', as.character(df)), type='error'
           )
        }

        enrich_data$distill$tbl <- df
        enrich_data$distill$title <- paste0(input$comp_fun, ' | ',
                                            input$geneset, ' | ',
                                            input$pathway)

      }) # eventReactive distill


      # observer to get distill enrichment table
      observeEvent(c(input$comp_tbl,
                     input$geneset_tbl,
                     input$pathway_tbl,
                     flags$data_loaded,
                     input$numcat.distill.tbl,
                     input$search_func_do,
                     input$subset_func_do), {

        validate(
          need(!is.null(app_object()$res), 'Waiting for selection')
        )

        validate(
          need(!is.null(app_object()$enrich), 'Functional enrichment results not available')
        )

        validate(
          need(input$comp_tbl != '' & input$geneset_tbl != '',
               'Waiting for selection')
        )

        isolate({
          updateSelectInput(session, 'comp_fun',
                            selected=input$comp_tbl)
          updateSelectInput(session, 'geneset',
                            selected=input$geneset_tbl)
          updateSelectInput(session, 'pathway',
                            selected=input$pathway_tbl)
        })

        if(is.null(input$numcat.distill.tbl))
          numcat <- config()$ui$functional_enrichment$plots$enrichment_map$numcat
        else numcat <- input$numcat.distill.tbl

        validate(
          need(numcat > 0, 'Number of terms must be > 0')
        )

        # get obj to visualize
        l <- enrich_data$enrich_list$df

        validate(
          need(nrow(as.data.frame(l)) > 0, 'No DE genes found!')
        )

        res <- enrich_data$enrich_list$res

        validate(
            need(!is.null(res), 'DE results not found for this geneset')
        )

        genetonic.list <- enrich_data$genetonic_list
        l_gs <- genetonic.list$l_gs
        anno_df <- genetonic.list$anno_df

        # subset table
        l_gs <- subset_func_tbl(l_gs, gene_col='gs_genes', desc_col='gs_description', gene_sep=',',
                                input$search_opts, input$search_txt,
                                input$subset_opts, input$upset_intersect, quiet=FALSE)

        n_gs <- min(nrow(l_gs), numcat)

        df <- tryCatch(
                distill_enrichment(l_gs,
                                 res,
                                 anno_df,
                                 n_gs = n_gs,
                                 cluster_fun = "cluster_markov"),
                warning = function(w){ w },
                error = function(e){ e }
              )

        enrich_data$distill$tbl <- df
        enrich_data$distill$title <- paste0(input$comp_tbl, ' | ',
                                            input$geneset_tbl, ' | ',
                                            input$pathway_tbl)

      }) # eventReactive distill

      # output distill table
      output$distill_table <- renderDT({
        df <- enrich_data$distill$tbl

        validate(
            need(!is.null(nrow(df$distilled_table)), 'Choose lower number of terms')
        )
        tbl <- df$distilled_table

        # rename columns
        # TODO: move to config
        colnames(tbl) <- c('cluster', '#_terms', 'genes', '#_genes',
                           'term_id_list', 'term_description_list', 'most_significant_term',
                           'strongest_term')

        cols.to.drop <- config()$server$functional_enrichment$table$distill_tbl$cols.to.drop

        tbl %>%
            mutate(genes=format_genes(.data$genes, genes.per.line=input$genes.per.line, sep=',')) %>%
            mutate(term_description_list=format_genes(.data$term_description_list, genes.per.line=5, sep=',')) %>%
            relocate('term_description_list', .after='strongest_term') %>%
            relocate('genes', .after='strongest_term') %>%
            select(-any_of(cols.to.drop)) %>%
            datatable(rownames=FALSE,
                      escape=FALSE,
                      options=list(
                      autoWidth=TRUE
                    ))
      }) # renderDT distill_table

      output$distill_table_title <- renderUI({
        fluidRow(align='center',
          h4(enrich_data$distill$title)
        )
      })

      # selection from distill table
      distilltable_proxy <- dataTableProxy('distill_table')

      observeEvent(input$reset_distilltable, {
        distilltable_proxy %>% selectRows(NULL)
      })

      observeEvent(input$add_distill_selected, {
        tbl <- enrich_data$distill$tbl$distilled_table
        sel <- input$distill_table_rows_selected

        s <- tbl$metags_genes

        if(is.null(sel)){
          showNotification(
            'Cannot add genes, no rows selected', type='warning'
          )
          validate(
            need(!is.null(sel), '')
          )
        } else {
          sel_genes <- unique(unlist(lapply(s[sel], function(x) strsplit(x, ',')[[1]])))
          if(all(sel_genes %in% genes_clicked$g)){
            showNotification(
              'Selected genes already present in scratchpad, skipping', type='warning'
            )
          } else {
            if(is.null(genes_clicked$g)) selected <- sel_genes
            else selected <- unique(c(genes_clicked$g, sel_genes))

            new_genes <- setdiff(selected, genes_clicked$g)
            showNotification(
              paste('Adding', length(new_genes),
                    'new genes to scratchpad')
            )

            genes_clicked$g <- selected
          }
        }
      })

      # fuzzy clustering
      observeEvent(c(input$comp_fun,
                     input$geneset,
                     input$pathway,
                     flags$data_loaded,
                     input$search_func_do,
                     input$subset_func_do), {
        validate(
          need(!is.null(app_object()$res), 'Waiting for selection')
        )

        validate(
          need(!is.null(app_object()$enrich), 'Functional enrichment results not available')
        )

        validate(
          need(input$comp_fun != '' & input$geneset != '',
               'Waiting for selection')
        )

        # get obj to visualize
        l <- enrich_data$enrich_list$df

        validate(
          need(nrow(as.data.frame(l)) > 0, 'No DE genes found!')
        )

        if(is.null(input$numcat.fuzzy.tbl))
          numcat <- config()$ui$functional_enrichment$plots$enrichment_map$numcat
        else numcat <- input$numcat.fuzzy.tbl

        validate(
          need(numcat > 0, 'Number of terms must be > 0')
          )

        res <- enrich_data$enrich_list$res

        genetonic.list <- enrich_data$genetonic_list
        l_gs <- genetonic.list$l_gs
        anno_df <- genetonic.list$anno_df

        # subset table
        l_gs <- subset_func_tbl(l_gs, gene_col='gs_genes', desc_col='gs_description', gene_sep=',',
                                input$search_opts, input$search_txt,
                                input$subset_opts, input$upset_intersect,
                                highlight=FALSE)


        subset_rows <- min(nrow(l_gs), numcat)

        fuzzy_clusters <- tryCatch(
                            gs_fuzzyclustering(l_gs[seq_len(subset_rows),],
                              # n_gs = nrow(res_enrich_subset),
                              # gs_ids = NULL,
                              # similarity_matrix = NULL,
                              similarity_threshold = 0.35,
                              fuzzy_seeding_initial_neighbors = 3,
                              fuzzy_multilinkage_rule = 0.5),
                            warning = function(w){ w },
                            error = function(e){ e })

        if(inherits(fuzzy_clusters, 'error')){
           showNotification(
             paste('Fuzzy clustering error:', as.character(fuzzy_clusters)), type='error'
           )
        } else {
          # add highlights
          fuzzy_clusters <- subset_func_tbl(fuzzy_clusters, gene_col='gs_genes', desc_col='gs_description', gene_sep=',',
                                  input$search_opts, input$search_txt,
                                  input$subset_opts, input$upset_intersect,
                                  quiet=TRUE)
        }

        enrich_data$fuzzy$tbl <- list(clusters=fuzzy_clusters,
                                      res=res, anno_df=anno_df)
        enrich_data$fuzzy$title <- paste0(input$comp_fun, ' | ',
                                          input$geneset, ' | ',
                                          input$pathway)

        isolate({
          updateSelectInput(session, 'comp_tbl',
                            selected=input$comp_fun)
          updateSelectInput(session, 'geneset_tbl',
                            selected=input$geneset)
          updateSelectInput(session, 'pathway_tbl',
                            selected=input$pathway)
        })

      }) # eventReactive fuzzy


      observeEvent(c(input$comp_tbl,
                     input$geneset_tbl,
                     input$pathway_tbl,
                     flags$data_loaded,
                     input$numcat.fuzzy.tbl,
                     input$search_func_do,
                     input$subset_func_do), {
        validate(
          need(!is.null(app_object()$res), 'Waiting for selection')
        )

        validate(
          need(!is.null(app_object()$enrich), 'Functional enrichment results not available')
        )

        validate(
          need(input$comp_tbl != '' & input$geneset_tbl != '',
               'Waiting for selection')
        )

        # get obj to visualize
        l <- enrich_data$enrich_list$df

        validate(
          need(nrow(as.data.frame(l)) > 0, 'No DE genes found!')
        )

        if(is.null(input$numcat.fuzzy.tbl))
          numcat <- config()$ui$functional_enrichment$plots$enrichment_map$numcat
        else numcat <- input$numcat.fuzzy.tbl

        validate(
          need(numcat > 0, 'Number of terms must be > 0')
          )

        res <- enrich_data$enrich_list$res

        genetonic.list <- enrich_data$genetonic_list
        l_gs <- genetonic.list$l_gs
        anno_df <- genetonic.list$anno_df

        # subset table
        l_gs <- subset_func_tbl(l_gs, gene_col='gs_genes', desc_col='gs_description', gene_sep=',',
                                input$search_opts, input$search_txt,
                                input$subset_opts, input$upset_intersect,
                                highlight=FALSE)

        subset_rows <- min(nrow(l_gs), numcat)

        fuzzy_clusters <- tryCatch(
                            gs_fuzzyclustering(l_gs[seq_len(subset_rows),],
                              # n_gs = nrow(res_enrich_subset),
                              # gs_ids = NULL,
                              # similarity_matrix = NULL,
                              similarity_threshold = 0.35,
                              fuzzy_seeding_initial_neighbors = 3,
                              fuzzy_multilinkage_rule = 0.5),
                            warning = function(w){ w },
                            error = function(e){ e })

        if(!inherits(fuzzy_clusters, 'error')){
          # add highlights
          fuzzy_clusters <- subset_func_tbl(fuzzy_clusters, gene_col='gs_genes', desc_col='gs_description', gene_sep=',',
                                  input$search_opts, input$search_txt,
                                  input$subset_opts, input$upset_intersect,
                                  quiet=TRUE)
        }

        enrich_data$fuzzy$tbl <- list(clusters=fuzzy_clusters,
                                      res=res, anno_df=anno_df)
        enrich_data$fuzzy$title <- paste0(input$comp_tbl, ' | ',
                                          input$geneset_tbl, ' | ',
                                          input$pathway_tbl)


        isolate({
          updateSelectInput(session, 'comp_fun',
                            selected=input$comp_tbl)
          updateSelectInput(session, 'geneset',
                            selected=input$geneset_tbl)
          updateSelectInput(session, 'pathway',
                            selected=input$pathway_tbl)
        })

      }) # eventReactive fuzzy

      # render fuzzy table
      output$fuzzy_table <- renderDT({
        l <- enrich_data$fuzzy$tbl

        tbl <- l$clusters
        validate(
          need(nrow(tbl) > 0, '')
        )
        colnames(tbl) <- sub('^gs_', '', colnames(tbl))

        float_idx <- vapply(tbl, function(x) typeof(x) %in% c('double', 'float'), logical(1))
        format_cols <- colnames(tbl)[float_idx]
        digits <- config()$server$functional_enrichment$table$fuzzy_tbl$format_significant$digits

        cols_to_format <- intersect(colnames(tbl), format_cols)

        tbl %>%
            mutate(genes=format_genes(.data$genes, genes.per.line=input$genes.per.line, sep=',')) %>%
            relocate('fuzzycluster', .before='pvalue') %>%
            relocate('cluster_status', .before='pvalue') %>%
            relocate('genes', .before='pvalue') %>%
            datatable(rownames=FALSE, escape=FALSE) %>%
          formatSignif(columns=cols_to_format,
                       digits=digits) %>%
          formatStyle('description', 'white-space'='nowrap')
      }) # renderDT fuzzy_table

      output$fuzzy_table_title <- renderUI({
        fluidRow(align='center',
          h4(enrich_data$fuzzy$title)
        )
      })

      # selection from distill table
      fuzzytable_proxy <- dataTableProxy('fuzzy_table')

      observeEvent(input$reset_fuzzytable, {
        fuzzytable_proxy %>% selectRows(NULL)
      })

      observeEvent(input$add_fuzzy_selected, {
        tbl <- enrich_data$fuzzy$tbl$clusters
        sel <- input$fuzzy_table_rows_selected

        s <- tbl$gs_genes

        if(is.null(sel)){
          showNotification(
            'Cannot add genes, no rows selected', type='warning'
          )
          validate(
            need(!is.null(sel), '')
          )
        } else {
          sel_genes <- unique(unlist(lapply(s[sel], function(x) strsplit(x, ',')[[1]])))
          if(all(sel_genes %in% genes_clicked$g)){
            showNotification(
              'Selected genes already present in scratchpad, skipping', type='warning'
            )
          } else {
            if(is.null(genes_clicked$g)) selected <- sel_genes
            else selected <- unique(c(genes_clicked$g, sel_genes))

            new_genes <- setdiff(selected, genes_clicked$g)
            showNotification(
              paste('Adding', length(new_genes),
                    'new genes to scratchpad')
            )

            genes_clicked$g <- selected
          }
        }
      })

     ######################### Functional enrichment plots ###########################
      output$enrichplot_title <- renderUI({
        txt <- enrich_data$enrich_list$title

        fluidRow(align='center',
          h4(txt)
        )
      })

      enrich_obj <- reactive({
        obj <- enrich_data$enrich_list$df

        # support gseaResult and enrichResult
        gene.id.col <- if('geneID' %in% colnames(obj)) 'geneID'
          else if('core_enrichment' %in% colnames(obj)) 'core_enrichment'

        obj <- subset_func_tbl(obj, gene_col=gene.id.col, desc_col='Description', gene_sep='/',
                               input$search_opts, input$search_txt,
                               input$subset_opts, input$upset_intersect, highlight=FALSE)

        validate(
          need(nrow(as.data.frame(obj)) > 0, 'No DE genes found!')
        )

        obj
      })

      genetonic_obj <- reactive({
        genetonic.list <- enrich_data$genetonic_list

        # subset before passing to plots
        genetonic.list$l_gs <- subset_func_tbl(genetonic.list$l_gs, gene_col='gs_genes',
                                               desc_col='gs_description', gene_sep=',',
                                               input$search_opts, input$search_txt,
                                               input$subset_opts, input$upset_intersect,
                                               highlight=FALSE)

        list(l_gs = genetonic.list$l_gs,
             anno_df = genetonic.list$anno_df,
             label = input$comp_fun)
      })

      # plot modules
      sumovPlotServer('sumov', genetonic_obj, config)
      enrichmapServer('enrichmap', genetonic_obj, enrich_obj, config)
      cnetPlotServer('cnetplot', enrich_obj, config)
      radarServer('radar', genetonic_obj, config)
      alluvialServer('alluvial', genetonic_obj,
                     enrich_obj, config)
      dendrogramServer('dendrogram', genetonic_obj, config)

      ##################### Distill plot ####################

      distill_args <- reactive({
        list(numcat=input$numcat.distill.tbl)
      })

      distill_data <- distillPlotServer('emap_distill',
                                        reactive({ enrich_data$distill$tbl }),
                                        distill_args,
                                        config)

      observeEvent(distill_data(), {
        numcat <- distill_data()

        if(numcat != input$numcat.distill.tbl){
          updateNumericInput(session, 'numcat.distill.tbl',
                             value=numcat)
        }
      })

      ##################### Fuzzy plot ####################

      fuzzy_args <- reactive({
        list(numcat=input$numcat.fuzzy.tbl)
      })

      fuzzy_data <- fuzzyPlotServer('emap_fuzzy',
                                    reactive({ enrich_data$fuzzy$tbl }),
                                    fuzzy_args,
                                    config)

      observeEvent(fuzzy_data(), {
        numcat <- fuzzy_data()

        if(numcat != input$numcat.fuzzy.tbl){
          updateNumericInput(session, 'numcat.fuzzy.tbl',
                             value=numcat)
        }
      })

      # output slot with enrichplot
      output$enrichplot <- renderUI({

          if(input$plottype == 'cnetplot'){

            cnetPlotUI(ns('cnetplot'), panel='main')

          } else if(input$plottype %in% c('radar')){

            radarUI(ns('radar'), panel='main')

          } else if(input$plottype %in% c('alluvial')){

            alluvialUI(ns('alluvial'), panel='main')

          } else if(input$plottype %in% c('dendrogram')){

            dendrogramUI(ns('dendrogram'), panel='main')

          } else if(input$plottype %in% c('enrichment_map')){

            enrichmapUI(ns('enrichmap'), panel='main')

          } else if(input$plottype %in% c('emap_distill')){

            distillPlotUI(ns('emap_distill'), panel='main')

          } else if(input$plottype %in% c('emap_fuzzy')){

            fuzzyPlotUI(ns('emap_fuzzy'), panel='main')

          } else if(input$plottype == 'summary_overview'){

            sumovPlotUI(ns('sumov'), panel='main')

          }
      }) # renderUI enrichplot

      # output slot with enrichplot
      output$enrichplot_btns <- renderUI({

          if(input$plottype == 'cnetplot'){

            cnetPlotUI(ns('cnetplot'), panel='main_btns')

          } else if(input$plottype %in% c('radar')){

            radarUI(ns('radar'), panel='main_btns')

          } else if(input$plottype %in% c('alluvial')){

            alluvialUI(ns('alluvial'), panel='main_btns')

          } else if(input$plottype %in% c('dendrogram')){

            dendrogramUI(ns('dendrogram'), panel='main_btns')

          } else if(input$plottype %in% c('enrichment_map')){

            enrichmapUI(ns('enrichmap'), panel='main_btns')

          } else if(input$plottype %in% c('emap_distill')){

            distillPlotUI(ns('emap_distill'), panel='main_btns')

          } else if(input$plottype %in% c('emap_fuzzy')){

            fuzzyPlotUI(ns('emap_fuzzy'), panel='main_btns')

          } else if(input$plottype == 'summary_overview'){

            sumovPlotUI(ns('sumov'), panel='main_btns')

          }
      }) # renderUI enrichplot


      #################### Plot options ############################

      output$plot_opts <- renderUI({
        if(input$plottype == 'cnetplot'){

          cnetPlotUI(ns('cnetplot'), panel='sidebar')

        } else if(input$plottype == 'summary_overview'){

          sumovPlotUI(ns('sumov'), panel='sidebar')

        } else if(input$plottype == 'enrichment_map'){

          enrichmapUI(ns('enrichmap'), panel='sidebar')

        } else if(input$plottype == 'emap_distill'){

          distillPlotUI(ns('emap_distill'), panel='sidebar')

        } else if(input$plottype == 'emap_fuzzy'){

          fuzzyPlotUI(ns('emap_fuzzy'), panel='sidebar')

        } else if(input$plottype == 'radar'){

          radarUI(ns('radar'), panel='sidebar')

        } else if(input$plottype == 'alluvial'){

          alluvialUI(ns('alluvial'), panel='sidebar')

        } else if(input$plottype == 'dendrogram'){

          dendrogramUI(ns('dendrogram'), panel='sidebar')

        }
      }) # renderUI

      ###################### Compare functional enrichment ############################

      # plot modules
      sumovPlotServer('sumov_comp', genetonic_comp_obj, config, type='comp')
      radarServer('radar_comp', genetonic_comp_obj, config, type='comp')
      horizonServer('horizon', genetonic_comp_obj, config)

      # enrichment objects
      genetonic_comp_obj <- reactive({
        # get this value to make sure DE gene checks are made
        obj <- enrich_comp_obj()

        genetonic.list1 <- app_object()$genetonic[[input$comp_fun1]][[input$geneset1]][[input$pathway1]]
        genetonic.list2 <- app_object()$genetonic[[input$comp_fun2]][[input$geneset2]][[input$pathway1]]

        # subset before passing to plots
        genetonic.list1$l_gs <- subset_func_tbl(genetonic.list1$l_gs, gene_col='gs_genes',
                                                desc_col='gs_description', gene_sep=',',
                                                input$search_opts, input$search_txt,
                                                input$subset_opts, input$upset_intersect,
                                                highlight=FALSE)
        genetonic.list2$l_gs <- subset_func_tbl(genetonic.list2$l_gs, gene_col='gs_genes',
                                                desc_col='gs_description', gene_sep=',',
                                                input$search_opts, input$search_txt,
                                                input$subset_opts, input$upset_intersect,
                                                highlight=FALSE)

        list(
          obj1=list(
             l_gs = genetonic.list1$l_gs,
             anno_df = genetonic.list1$anno_df,
             label = input$comp_fun1,
             geneset = input$geneset1
          ),
          obj2=list(
             l_gs = genetonic.list2$l_gs,
             anno_df = genetonic.list2$anno_df,
             label = input$comp_fun2,
             geneset = input$geneset2
          )
        )
      })

      enrich_comp_obj <- reactive({
        validate(
          need(input$comp_fun2 != 'Choose one', 'Please select both comparisons to continue ...')
        )
        obj1 <- app_object()$enrich[[input$comp_fun1]][[input$geneset1]][[input$pathway1]]
        obj2 <- app_object()$enrich[[input$comp_fun2]][[input$geneset2]][[input$pathway1]]

        validate(
          need(nrow(as.data.frame(obj1)) > 0,
               'No DE genes found in Comparison 1!')
        )

        validate(
          need(nrow(as.data.frame(obj2)) > 0,
               'No DE genes found in Comparison 2!')
        )

        list(
          obj1=obj1, obj2=obj2
        )
      })

      observeEvent(input$swap_comp, {
        if(input$comp_fun1 != input$comp_fun2){
          updateSelectInput(session, 'comp_fun1',
                            selected=input$comp_fun2)
          updateSelectInput(session, 'comp_fun2',
                            selected=input$comp_fun1)
        } else {
          showNotification(
            'Comparison 1 is the same as Comparison 2, no changes made'
          )
        }
      })

      observeEvent(input$swap_geneset, {
        if(input$geneset1 != input$geneset2){
          updateSelectInput(session, 'geneset1',
                            selected=input$geneset2)
          updateSelectInput(session, 'geneset2',
                            selected=input$geneset1)
        } else {
          showNotification(
            'Gene set 1 is the same as gene set 2, no changes made'
          )
        }
      })
      # observer to update geneset & pathway menus
      observeEvent(input$comp_fun1, {
        validate(
          need(input$comp_fun1 != '', 'No comparisons selected')
        )

        validate(
          need(!is.null(app_object()$enrich), 'Functional enrichment results not available')
        )

        validate(
          need(input$comp_fun1 %in% names(app_object()$enrich),
               'Functional enrichment results not available for this comparison')
        )

        elements <- names(app_object()$enrich[[input$comp_fun1]])

        # skip 'res' element
        elements <- setdiff(elements, 'res')

        # if selected geneset1 is not in elements, reset to 1st element
        # else, keep selection
        if(!input$geneset1 %in% elements){
            # choose pathway results to show
            updateSelectInput(session,
                              'geneset1',
                              choices=elements,
                              selected=elements[1])
        } else {
            # choose pathway results to show
            updateSelectInput(session,
                              'geneset1',
                              choices=elements,
                              selected=input$geneset1)
        }

        comp.choices2 <- c('Choose one', names(app_object()$enrich))

        # if selected comp_fun2 is not in choices, reset to 1st element
        # else, keep selection
        if(!input$comp_fun2 %in% comp.choices2){
            selected <- comp.choices2[1]
        } else {
            selected <- input$comp_fun2
        }

        updateSelectInput(session,
                          'comp_fun2',
                          choices=comp.choices2,
                          selected=selected)
      }) # observeEvent

      # observer to update geneset2 when comparison 2 is selected
      observeEvent(input$comp_fun2, {

        validate(
          need(input$comp_fun2 %in% names(app_object()$enrich),
               'Functional enrichment results not available for this comparison')
        )

        elements2 <- names(app_object()$enrich[[input$comp_fun2]])
        elements2 <- setdiff(elements2, 'res')

        # if selected geneset2 is not in elements2, reset to 1st element
        # else, keep selection
        if(!input$geneset2 %in% elements2){
          # if comparing within contrast, skip geneset1
          if(input$comp_fun1 == input$comp_fun2){
            selected <- setdiff(elements2, input$geneset1)[1]
          } else {
            selected <- elements2[1]
          }
        } else {
            selected <- input$geneset2
        }

        updateSelectInput(session,
                          'geneset2',
                          choices=elements2,
                          selected=selected)
      }) # observeEvent
      # observer to update pathway menu on comparison/geneset selection
      observeEvent(c(input$comp_fun1, input$geneset1,
                     input$comp_fun2, input$geneset2), {
        # get names of available pathway analyses
        # named vector of pathway names
        p <- unlist(config()$server$functional_enrichment$pathways)

        # NOTE: only pathway results in both sets can be compared
        all.p.names <- intersect(names(app_object()$enrich[[input$comp_fun1]][[input$geneset1]]),
                                 names(app_object()$enrich[[input$comp_fun2]][[input$geneset2]]))

        validate(
            need(length(all.p.names) > 0, 'Comparisons must share at least 1 pathway analysis')
        )

        p.choices <- p[p %in% all.p.names]
        # if database names are not in known options
        if(length(p.choices) < length(all.p.names)){
            p.choices <- c(p.choices, all.p.names[!all.p.names %in% p])
        }

        if(!input$pathway1 %in% p.choices){
            updateSelectInput(session,
                              'pathway1',
                              choices=p.choices,
                              selected=p.choices[1])
        } else {
            updateSelectInput(session,
                              'pathway1',
                              choices=p.choices,
                              selected=input$pathway1)
        }

      }) # observeEvent


      # reactive to render compare enrichment plot
      output$compare_enrich <- renderUI({
        #p <- compenrich()

        if(input$comptype %in% c('summary_overview')){

          sumovPlotUI(ns('sumov_comp'), panel='main')

        } else if(input$comptype %in% c('radar')){

          radarUI(ns('radar_comp'), panel='main')

        } else if(input$comptype %in% c('horizon')){

          horizonUI(ns('horizon'), panel='main')

        }
      }) # renderUI compare_enrich

      # reactive to render compare enrichment plot
      output$compplot_btns <- renderUI({
        #p <- compenrich()

        if(input$comptype %in% c('summary_overview')){

          sumovPlotUI(ns('sumov_comp'), panel='main_btns', type='comp')

        } else if(input$comptype %in% c('radar')){

          radarUI(ns('radar_comp'), panel='main_btns', type='comp')

        } else if(input$comptype %in% c('horizon')){

          horizonUI(ns('horizon'), panel='main_btns')

        }
      }) # renderUI compare_enrich

      output$compenrich_title <- renderUI({
        validate(
          need(input$comp_fun2 != 'Choose one', '')
        )
        txt1 <- paste0(input$comp_fun1, ' | ',
                      input$geneset1, ' | ',
                      input$pathway1, '   vs')
        txt2 <- paste0(input$comp_fun2, ' | ',
                      input$geneset2, ' | ',
                      input$pathway1)

        fluidRow(align='center',
          h4(txt1), h4(txt2)
        )

      })

      #################### Plot options ############################

      # compare results plot opts
      output$cmp_plot_opts <- renderUI({
        if(input$comptype == 'summary_overview'){

          sumovPlotUI(ns('sumov_comp'), panel='sidebar')

        } else if(input$comptype == 'radar'){

          radarUI(ns('radar_comp'), panel='sidebar')

        } else if(input$comptype == 'horizon'){

          horizonUI(ns('horizon'), panel='sidebar')

        }
      }) # renderUI

      ######################## buttons #######################

      # functional enrichment controls help
      helpButtonServer('tbl_ctrls_help', size='l')
      helpButtonServer('func_plt_help', size='l')
      helpButtonServer('func_cmp_help', size='l')

      # functional enrichment table help
      helpButtonServer('func_tbl_help', size='l')
      helpButtonServer('distill_tbl_help', size='l')
      helpButtonServer('fuzzy_tbl_help', size='l')

      # build reactive to return
      # - this is only updated when selected genes are changed
      return_data <- eventReactive(c(genes_clicked$g), {
        list(
          genes=genes_clicked$g
        )
      })

      return(return_data)

    } # function
  ) # moduleServer
} # enrichServer


