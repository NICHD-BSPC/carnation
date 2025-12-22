#' Upset plot module
#'
#' @description
#' Module UI & server to generate upset plots.
#'
#' @param id Module id
#' @param panel string, can be 'sidebar' or 'main'
#' @param obj reactiveValues object containing carnation object
#' @param plot_args reactive containing 'fdr.thres' (padj threshold) & 'fc.thres' (log2FC)
#' @param gene_scratchpad reactiveValues object containing genes selected in scratchpad
#' @param reset_genes reactive to reset gene scratchpad selection
#' @param config reactive list with config settings
#'
#' @returns
#' UI returns tagList with upset plot UI.
#' Server returns reactive with list containing upset table, intersections
#' & selected genes.
#'
#' @rdname upsetmod
#' @name upsetmod
NULL

#' @rdname upsetmod
#' @export
upsetPlotUI <- function(id, panel){
  ns <- NS(id)

  # get default config
  config <- get_config()

  if(panel == 'sidebar'){
    tag <- tagList(
             # reactive ui for selected comparisons bucket list
             fluidRow(
               column(6, strong('Comparisons')),
               column(6, align='right',
                 helpButtonUI(ns('upset_controls_help'))
               ) # column
             ), # fluidRow

             uiOutput(ns('upset_bucket')),

             fluidRow(style='margin-bottom: 20px',
               column(6, align='right',
                 actionButton(ns('upset_all'), 'Select all')
               ), # column
               column(6, align='left',
                 actionButton(ns('upset_none'), 'Select none')
               ) # column
             ), # fluidRow

             fluidRow(align='center',
               actionButton(ns('updt_do'),
                            label='Refresh',
                            icon=icon('arrows-rotate'),
                            class='btn-primary',
                            style='margin-bottom: 10px;')
             ), # fluidRow

            bsCollapse(
               bsCollapsePanel('UpSet options',
                 fluidRow(
                   column(4, h5('Direction of change')),
                   column(8,
                     selectInput(ns('upset_type'), label=NULL,
                                 choices=c('up & down (changed)'='changed',
                                           'up'='up',
                                           'down'='down',
                                           'custom'='custom')
                     ) # selectInput
                   ) # column
                 ), # fluidRow

                 conditionalPanel(
                   paste0('input["', ns('upset_type'), '"] == "custom"'),
                   fluidRow(align='right',
                     column(12,
                       actionButton(ns('edit_custom'),
                                    label='Customize',
                                    icon=icon('edit'),
                                    style='margin-bottom: 10px;')
                     ) # column
                   ) # fluidRow
                 ), # conditionalPanel

                 conditionalPanel(
                   paste0('input["', ns('tab_selected'), '"] == "Table"'),
                   fluidRow(
                     column(12,
                       selectizeInput(ns('upset_intersections'),
                                      label=h4('Intersections'),
                                      choices=NULL, selected=NULL,
                                      multiple=TRUE
                       ) # selectInput
                     ) # column
                   ), # fluidRow

                   fluidRow(
                     column(6, align='center',
                       actionButton(ns('intersect_all'), 'Select all')
                     ), # column
                     column(6, align='center',
                       actionButton(ns('intersect_none'), 'Select none')
                     ) # column
                   ) # fluidRow
                 ), # conditionalPanel

                 conditionalPanel(
                   paste0('input["', ns('tab_selected'), '"] == "Plot"'),
                   fluidRow(
                     column(4, h5('# of intersections')),
                     column(8,
                       numericInput(ns('n_intersections'), label=NULL,
                         value=config$ui$de_analysis$upset_plot$n_intersections,
                         min=0, step=10
                       ) # numericInput
                     ) # column
                   ), # fluidRow

                   fluidRow(
                     column(4, h5('Min intersection size')),
                     column(8,
                       numericInput(ns('min_size'), label=NULL,
                         value=config$ui$de_analysis$upset_plot$min_size,
                         min=0, step=1
                       ) # numericInput
                     ) # column
                   ), # fluidRow

                   fluidRow(
                     column(4, h5('Text scale')),
                     column(8,
                       numericInput(ns('text_scale'), label=NULL,
                         value=config$ui$de_analysis$upset_plot$text_scale,
                         min=0.1, step=0.1
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
                 ) # conditionalPanel

               ) # bsCollapsePanel

             ) # bsCollapse

           ) # tagList
  } else if(panel == 'main'){
    tag <- tagList(
             tabsetPanel(id=ns('tab_selected'),
               tabPanel('Plot',
                 fluidRow(
                    column(6, align='left',
                     helpButtonUI(ns('de_upset_help'))
                   ), # column
                   column(6, align='right',
                     downloadButtonUI(ns('upset_download'))
                   ) # column
                 ), # fluidRow

                 withSpinner(
                   plotOutput(ns('upsetplot'), height='800px')
                 ), # withSpinner
                 DTOutput(ns('upset_labeled'))

               ), # tabPanel

               tabPanel('Table',
                 fluidRow(
                    column(6, align='left',
                     helpButtonUI(ns('de_upset_tbl_help'))
                   ) # column
                 ), # fluidRow

                 fluidRow(
                   column(12,
                     DTOutput(ns('intertable'))
                   ) # column
                 ), # fluidRow

                 fluidRow(align='center',
                   br(),
                   column(12,
                     splitLayout(
                       cellWidths=c('15%', '25%', '20%'),

                       strong('Selection options'),
                       actionButton(ns('add_selected'), 'Add to scratchpad'),
                       actionButton(ns('reset_intertable'),
                                    'Reset selection',
                                    class='btn-primary')
                     ) # splitLayout
                   ) # column
                 ) # fluidRow

               ) # tabPanel
             ) # tabsetPanel

           ) # tagList
  }

}

#' @rdname upsetmod
#' @export
upsetPlotServer <- function(id, obj, plot_args, gene_scratchpad, reset_genes, config){

  moduleServer(
    id,

    function(input, output, session){

      ns <- NS(id)

      # TODO: move to config
      comp_split_pattern <- ';'

      app_object <- reactive({
        list(res=obj$res)
      })

      # reactive values to hold upset choices
      upset_choices <- reactiveValues(current=NULL, all=NULL, staged=NULL, tbl=NULL, names=NULL, dir_list=list())

      custom_refresh <- reactiveValues(flag=0)

      # reactive values to hold current upset table
      upset_table <- reactiveValues(tbl=NULL,
                                    set_mapping=NULL,
                                    set_labels=NULL)

      # reactive values to track genes clicked/labeled
      genes_clicked <- reactiveValues(g=NULL)

      reset_data <- function(){
        upset_choices$current <- NULL
        upset_choices$all <- NULL
        upset_choices$staged <- NULL
        upset_choices$tbl <- NULL
        upset_choices$names <- NULL
        upset_choices$dir_list <- list()
        upset_table$tbl <- NULL
        upset_table$set_mapping <- NULL
        upset_table$set_labels <- NULL

        genes_clicked$g <- NULL
      }

      # update from reactive config
      observeEvent(config(), {
        updateNumericInput(session, 'n_intersections',
                           value=config()$ui$de_analysis$upset_plot$n_intersections)
        updateNumericInput(session, 'min_size',
                           value=config()$ui$de_analysis$upset_plot$min_size)
        updateNumericInput(session, 'text_scale',
                           value=config()$ui$de_analysis$upset_plot$text_scale)
      })

      # observer to initialize/reset upset_plot choices
      # NOTE: this is only run once per loaded assay
      observeEvent(app_object()$res, {
        reset_data()
        validate(
          need(!is.null(app_object()$res), 'Waiting for selection')
        )
        choices <- names(app_object()$res)
        upset_choices$all <- choices

        comp_num <- config()$server$de_analysis$upset_plot$comp_num
        if(length(choices) > comp_num){
          choices <- choices[seq_len(comp_num)]
        }

        upset_choices$current <- choices
        upset_choices$staged <- choices
        upset_choices$tbl <- data.frame(row.names=choices,
                                        changed=rep('', length(choices)),
                                        up=rep('', length(choices)),
                                        down=rep('', length(choices)))
        upset_choices$names <- choices

        custom_tbl_proxy %>% selectCells(NULL)
        updateSelectInput(session, 'upset_type',
                          selected='changed')

      })

      output$custom_upset <- renderDT({
        df <- upset_choices$tbl
        df <- df[rownames(df) %in% upset_choices$current,]

        df %>%
          datatable(selection=list('multiple', target='cell'),
                    caption=tags$caption(style='font-weight: bold; font-size: 15px;',
                                         'Customize directions'),
                    options=list(dom='t', pageLength=length(app_object()$res)))
      })

      custom_tbl_proxy <- dataTableProxy('custom_upset')

      observeEvent(input$upset_all, {
        validate(
          need(!is.null(app_object()$res), 'Waiting for selection')
        )

        upset_choices$staged <- upset_choices$all
      })


      # observer to select none
      observeEvent(input$upset_none, {
        validate(
          need(!is.null(app_object()$res), 'Waiting for selection')
        )
        upset_choices$staged <- NULL
      })

      # observer to reset data for empty selection & chng current choices
      observeEvent(input$upset_comps, {
        if(length(input$upset_comps) < 2){
          updateSelectizeInput(session, 'upset_intersections',
                               choices='',
                               server=TRUE)
          upset_table$tbl <- NULL
        }

        upset_choices$staged <- input$upset_comps
      })

      observeEvent(input$updt_do, {
        upset_choices$current <- upset_choices$staged
      })

      observeEvent(upset_choices$current, {
        choices <- upset_choices$current

        prev_tbl <- upset_choices$tbl

        new_tbl <- data.frame(row.names=choices,
                              changed=rep('', length(choices)),
                              up=rep('', length(choices)),
                              down=rep('', length(choices)))
        for(rn in rownames(prev_tbl)){
          new_tbl[rn, ] <- prev_tbl[rn, ]
        }
        upset_choices$tbl <- new_tbl
        upset_choices$names <- choices

      })

      # upset plot bucket list observer
      output$upset_bucket <- renderUI({
        tagList(
          fluidRow(
            column(12,
              bucket_list(
                header=NULL,
                group_name='upset_group',
                class=c('default-sortable','custom-sortable'),
                orientation='horizontal',
                add_rank_list(
                  text = 'current',
                  labels = upset_choices$staged,
                  input_id=ns('upset_comps')
                ), # add_rank_list
                add_rank_list(
                  text = 'unused',
                  labels = setdiff(names(app_object()$res),
                                   upset_choices$staged),
                  input_id=ns('upset_other')
                ) # add_rank_list
              ) # bucket_list
            ) # column
          ) # fluidRow
        ) # tagList
      }) # renderUI

      observeEvent(input$edit_custom, {
        showModal(
          modalDialog(
            DTOutput(ns('custom_upset')),
            footer=tagList(
              actionButton(ns('custom_ok'), label='Apply'),
              actionButton(ns('custom_reset'), label='Reset'),
              modalButton('Cancel'),
              actionButton(ns('custom_do'),
                           label='Plot',
                           icon=icon('arrows-rotate'),
                           class='btn-primary')
            )
          )
        )
      })

      observeEvent(input$custom_do, {
        removeModal()
      })

      observeEvent(input$custom_reset, {
        for(i in seq_len(nrow(upset_choices$tbl))){
          upset_choices$tbl[i, ] <- ''
        }
        custom_tbl_proxy %>% selectCells(NULL)
      })

      observeEvent(input$custom_ok, {

        # convert selected cells to comparison - direction
        # direction order: changed, up, down
        sel_mat <- input$custom_upset_cells_selected

        dir_list <- list()
        # get current selections
        for(cmp in rownames(upset_choices$tbl)){
          cmp.i <- upset_choices$tbl[cmp,] != ''
          if(sum(cmp.i) > 0){
            dir_list[[ cmp ]] <- colnames(upset_choices$tbl)[ which(cmp.i)[1] ]
          }
        }

        # add new selections
        if(!is.null(sel_mat)){
          if(nrow(sel_mat) > 1 | length(dir_list) > 1){

            comp_names <- rownames(upset_choices$tbl)
            dir_names <- colnames(upset_choices$tbl)

            # decode selections
            comp_sel <- comp_names[ sel_mat[,1] ]
            dir_sel <- dir_names[ sel_mat[,2] ]

            # update current selections
            for(cmp in unique(comp_sel)){
              idx <- comp_sel %in% cmp
              dir.i <- dir_sel[idx]

              # - if 'changed' is selected, that takes precedence
              # - if 'up' and 'down' are both selected, it's converted to 'changed'
              if('changed' %in% dir.i) dir.i <- 'changed'
              else if(all(c('up','down') %in% dir.i)) dir.i <- 'changed'
              dir_list[[ cmp ]] <- dir.i
            }

            if(length(dir_list) == 1){
              showNotification(
                'Must select at least two different comparisons!',
                type='warning'
              )

              validate(
                need(length(dir_list) > 1, '')
              )
            }

            # empty tbl
            new_tbl <- upset_choices$tbl
            for(i in seq_len(nrow(upset_choices$tbl))){
              new_tbl[i, ] <- ''
            }

            # update selection table
            upset_comp_names <- NULL
            for(cmp in names(dir_list)){
              new_tbl[cmp, dir_list[[ cmp ]] ] <- '\u2713'
              upset_comp_names <- c(upset_comp_names, paste0(cmp, '_', dir_list[[ cmp ]]))
            }

            upset_choices$tbl <- new_tbl
            upset_choices$names <- upset_comp_names
            upset_choices$dir_list <- dir_list

            custom_refresh$flag <- custom_refresh$flag + 1
          } else if(nrow(sel_mat) == 1){
            showNotification(
              'Must select at least two different comparisons!',
              type='warning'
            )
          }
        }
      })

      # reactive function to return list of gene sets for
      # upset plot
      upset_gene_lists <- eventReactive(c(plot_args()$fdr.thres,
                                        plot_args()$fc.thres,
                                        custom_refresh$flag,
                                        upset_choices$tbl,
                                        input$plot_do,
                                        input$updt_do), {

        validate(
          need(!is.null(app_object()$res), 'Waiting for selection')
        )

        fc.thres <- ifelse(plot_args()$fc.thres == '' | is.na(plot_args()$fc.thres),
                           config()$ui$de_analysis$filters$log2fc_threshold, plot_args()$fc.thres)
        fdr.thres <- ifelse(plot_args()$fdr.thres == '' | is.na(plot_args()$fdr.thres),
                            config()$ui$de_analysis$filters$fdr_threshold, plot_args()$fdr.thres)

        # gene lists to compare
        # TODO: move to function?
        upset_comps <- upset_choices$current
        validate(
          need(length(upset_comps) > 1 & upset_comps %in% names(app_object()$res),
               'Must select at least two gene lists to compare!')
        )

        # store upset_type
        upset_type <- input$upset_type

        dir_list <- list()
        if(upset_type == 'custom'){

          dir_list <- upset_choices$dir_list
          upset_comp_names <- upset_choices$names

        } else {
          upset_comp_names <- paste0(upset_comps, '_', upset_type)
        }

        gene.lists <- lapply(upset_comps, function(u){
          x <- app_object()$res[[ u ]]
          idx <- x$padj < fdr.thres & !is.na(x$padj)
          if(upset_type == 'up') idx <- idx & x$log2FoldChange > fc.thres
          else if(upset_type == 'down') idx <- idx & x$log2FoldChange < -fc.thres
          else if(upset_type == 'changed') idx <- idx & abs(x$log2FoldChange) > fc.thres
          else if(upset_type == 'custom'){
            if(length(dir_list) == 0) return(character(0))
            else if(!u %in% names(dir_list)) return(character(0))

            direc <- dir_list[[ u ]]
            if(direc == 'up') idx <- idx & x$log2FoldChange > fc.thres
            else if(direc == 'down') idx <- idx & x$log2FoldChange < -fc.thres
            else if(direc == 'changed') idx <- idx & abs(x$log2FoldChange) > fc.thres
          }

          if(sum(idx, na.rm=TRUE) == 0) return(character(0))

          # TODO: remove this section
          if('symbol' %in% tolower(colnames(x))){
            sidx <- tolower(colnames(x)) %in% c('symbol')
            s <- as.character(x[,sidx])
            s[is.na(s)] <- rownames(x)[is.na(s)]

            setNames(s[idx], rownames(x)[idx])
          } else {
            setNames(rownames(x)[idx], rownames(x)[idx])
          }
        })

        names(gene.lists) <- upset_comps

        # this is when you click 'Cancel' from 'Customize' dialog, without any prior selections
        if(upset_type == 'custom' & length(dir_list) == 0){
          showNotification(
            'Direction of change is "custom", but no selections made! Please make at least two custom selections and refresh',
            type='warning'
          )
        } else if(length(upset_comp_names) < length(gene.lists)){
          gene.lists <- gene.lists[ names(dir_list) ]
          names(gene.lists) <- upset_comp_names
        } else {
          names(gene.lists) <- upset_comp_names
        }

        gene.lists
      })

      # reactive to generate upset table & plot from gene lists
      observeEvent(c(upset_gene_lists(),
                     upset_choices$current), {
        validate(
          need(!is.null(app_object()$res), 'Waiting for selection')
        )

        gene.lists <- upset_gene_lists()

        # check length of gene lists and skip empty gene lists
        glen <- unlist(lapply(gene.lists, function(x) length(x)))
        if(any(glen == 0)){
            zero.idx <- which(glen == 0)
            if(input$upset_type != 'custom'){
              showNotification(
                  paste('Upset plot warning: Skipping', length(zero.idx), 'gene lists with 0 genes:',
                        paste(names(gene.lists)[zero.idx], collapse=', '))
              )
            }
            gene.lists <- gene.lists[-zero.idx]
            if(length(gene.lists) == 0){
              upset_table$tbl <- NULL
            }
            validate(
              need(length(gene.lists) > 1,
                   'Have less than two gene lists to compare!')
            )
        }

        # get matrix of intersections, add set column & save
        df <- fromList.with.names(gene.lists)
        df <- add.set.column(df)

        # get mapping of set names to intersects
        df_unique <- unique(df[, setdiff(colnames(df), 'symbol')])
        rownames(df_unique) <- df_unique$set
        df_unique <- df_unique[, setdiff(colnames(df_unique), 'set')]
        if(nrow(df_unique) == 1){
          set_mapping <- setNames(list(colnames(df_unique)), rownames(df_unique))
        } else {
          set_mapping <- apply(df_unique, 1, function(x){
                           n <- colnames(df_unique)[x == 1]
                           n
              })
        }
        upset_table$set_mapping <- set_mapping

        # get columns with degree & comparisons & add to df
        comps <- unlist(
                   lapply(set_mapping[df$set],
                     function(x) paste(x, collapse=comp_split_pattern)
                   )
                 )
        degree <- unlist(lapply(set_mapping[df$set], length))

        upset_table$tbl <- cbind(df, comparisons=comps, degree=degree)

        # get intersection sizes
        inter_counts <- table(df$set)

        # create named vector of choices & add groups + size to names
        # e.g. if 'set1' has 150 genes,
        # then the name is 'set1 (n = 150)'
        inter_choices <- names(inter_counts)
        names(inter_choices) <- paste0(inter_choices, ' (',
                                unlist(
                                  lapply(set_mapping,
                                    function(x)
                                      paste(x, collapse=', ')
                                  )
                                ), '; n = ', inter_counts, ')')

        upset_table$set_labels <- inter_choices
        if(!is.null(input$upset_intersections) & all(input$upset_intersections %in% inter_choices)) selected <- input$upset_intersections
        else selected <- inter_choices

        intersect_num <- config()$server$de_analysis$upset_plot$intersect_num
        if(length(selected) > intersect_num){
          selected <- selected[seq_len(intersect_num)]
        }
        updateSelectizeInput(session, 'upset_intersections',
                             choices=inter_choices,
                             selected=selected,
                             server=TRUE
        )

      }) # observeEvent

      # observer to update genes_clicked
      observeEvent(gene_scratchpad(), {

        if(any(gene_scratchpad() != ''))
          genes_clicked$g <- gene_scratchpad()

      })

      observeEvent(reset_genes(), {
        genes_clicked$g <- NULL
      })

      observeEvent(input$intersect_none, {
        updateSelectizeInput(session, 'upset_intersections',
                             selected='')
      })

      observeEvent(input$intersect_all, {
        inter_choices <- upset_table$set_labels
        if(length(inter_choices) > 100){
          showNotification(
            'Warning: large number of intersections selected. This can take a while ...',
            type='warning'
          )
        }
        updateSelectizeInput(session, 'upset_intersections',
                             selected=inter_choices)
      })

      intersect_tbl <- eventReactive(c(plot_args()$fdr.thres,
                                       plot_args()$fc.thres,
                                       app_object()$res,
                                       input$custom_do,
                                       input$updt_do,
                                       input$plot_do,
                                       input$upset_intersections), {
        validate(
          need(!is.null(app_object()$res), 'Waiting for selection')
        )

        validate(
          need(!is.null(upset_table$set_labels), 'Loading ...')
        )

        validate(
          need(!is.null(upset_table$tbl), '')
        )

        validate(
          need(length(upset_choices$current) > 1,
               'Must select at least two gene lists to compare!')
        )

        validate(
          need(!is.null(input$upset_intersections),
               'Please select intersection(s) to view here')
        )

        gdf <- upset_table$tbl
        gdf <- gdf[order(rownames(gdf)),]

        # attach results columns
        res_list <- app_object()$res
        tmp <- lapply(res_list[names(res_list) %in% upset_choices$current],
                 function(x){
                   idx <- rownames(x) %in% rownames(gdf)

                   # save symbol column
                   df <- x[idx, c('log2FoldChange', 'padj')]

                   # if some rows are missing, impute NA's
                   if(nrow(df) < nrow(gdf)){
                     add.rows <- setdiff(rownames(gdf), rownames(df))
                     df[add.rows, ] <- NA
                   }

                   # order by rownames - this is needed for
                   # combining data frames
                   df[order(rownames(df)),]
                 }
               )
        tmp <- do.call('cbind', tmp)
        tmp <- tmp[order(rownames(tmp)),]

        # get final df by combining
        final_df <- cbind(gdf, tmp)

        # filter df by intersections & upset_choices
        ridx <- final_df$set %in% input$upset_intersections
        cidx <- which(!colnames(final_df) %in% upset_choices$current)
        final_df <- final_df[ridx, cidx]
        final_df <- cbind(gene=rownames(final_df),
                          final_df)

        # sort by set and then first LFC column
        #lfc_idx <- grep('log2FoldChange', colnames(final_df))
        #final_df <- final_df[order(final_df$set,
        #                           final_df[, lfc_idx[1]]),]
        #final_df
        df2 <- do.call('rbind',
                 lapply(unique(final_df$set), function(x){
                   tmp <- final_df[final_df$set == x,]
                   comps <- unlist(strsplit(
                              tmp$comparisons[1], comp_split_pattern))
                   lfc_idx <- which(colnames(tmp) %in% paste0(comps, '.log2FoldChange'))
                   if(length(lfc_idx) > 0){
                     tmp <- tmp[order(tmp[, lfc_idx[1]],
                                      decreasing=TRUE),]
                   }
                   tmp
               }))
        df2
      })

      output$intertable <- renderDT({

        final_df <- intersect_tbl()

        isolate({
          curr_choices <- upset_choices$current
          nchoices <- length(curr_choices)
          cnames <- setdiff(colnames(upset_table$tbl), curr_choices)
          curr_names <- upset_choices$names
        })

        # create custom container
        sketch <- withTags(table(
          class = 'display',
          tags$thead(
            tags$tr(
              tags$th(rowspan=2, 'gene'),
              lapply(cnames,
                     function(x) tags$th(rowspan=2, x)),
              lapply(curr_names,
                     function(x) tags$th(class='dt-center', colspan=2, x))
            ),
            tags$tr(
              lapply(rep(c('log2FoldChange','padj'), nchoices), tags$th)
            )
          )
        ))

        format_idx <- unique(
                        unlist(lapply(c('log2FoldChange', 'padj'),
                          function(x) grep(x, colnames(final_df)))
                        )
                      )

        border_cols <- grep('log2FoldChange', colnames(final_df))

        final_df %>%
            datatable(rownames=FALSE,
              container=sketch,
              options=list(
                columnDefs=list(list(className='dt-center',
                                     targets=format_idx-1)))
              ) %>%
            formatSignif(columns=format_idx,
                         digits=4) %>%
            formatStyle(columns=border_cols,
                      'border-left'='solid 1px')
      })

      intertable_proxy <- dataTableProxy('intertable')

      observeEvent(input$reset_intertable, {
        intertable_proxy %>% selectRows(NULL)
      })

      observeEvent(input$add_selected, {
        tbl <- intersect_tbl()
        sel <- input$intertable_rows_selected

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

          genes_clicked$g <- selected
        }
      })

      # reactive function to generate upset plot
      upsetplot <- eventReactive(c(plot_args()$fdr.thres,
                                   plot_args()$fc.thres,
                                   app_object()$res,
                                   input$custom_do,
                                   input$updt_do,
                                   input$plot_do), {

        validate(
          need(!is.null(app_object()$res), 'Waiting for selection')
        )

        validate(
          need(!is.null(upset_table$tbl) & length(upset_choices$current) > 1,
               'Must select at least two gene lists to compare!\n\nIf you\'ve already done so, try refreshing the plot from the settings menu')
        )

        gdf <- upset_table$tbl

        intersect_sets <- unique(
                            unlist(
                              strsplit(unique(gdf$comparisons),
                                       split=';', fixed=TRUE)
                            )
                          )

        # only keep intersect_sets present in gdf
        intersect_sets <- intersect(intersect_sets,
                                    colnames(gdf))

        ts <- input$text_scale

        if(input$min_size < 0) min_size <- 0
        else min_size <- input$min_size

        if(input$n_intersections < 0) n_intersections <- 0
        else n_intersections <- input$n_intersections

        inter_sizes <- table(gdf$set)
        if(max(inter_sizes) < min_size){
          validate(
            need(max(inter_sizes) >= min_size,
                 paste0('No intersections left after filtering!\n\n',
                        'Please adjust "Min intersection size" in settings menu to at least ', max(inter_sizes))
            )
          )
        }

        p <- upset(gdf, intersect=intersect_sets,
               width_ratio=0.2,
               height_ratio=0.2,
               base_annotations=list(
                 'Intersection size'=intersection_size(
                                       text=list(size=5*ts,
                                                 angle=90,
                                                 hjust=-0.2,
                                                 vjust=0.5))
               ),
               min_size=min_size,
               n_intersections=n_intersections,
               themes=upset_modify_themes(
                 list(
                   'Intersection size'=
                       theme(axis.title.y=element_text(size=20*ts),
                             axis.text.y=element_text(size=15*ts)),
                   'intersections_matrix'=
                       theme(axis.text.y=element_text(size=15*ts),
                             axis.title.x=element_blank(),
                             panel.grid.major=element_blank()),
                   'overall_sizes'=
                       theme(axis.text.x=element_text(size=15*ts,
                                                      hjust=1,
                                                      angle=90),
                             axis.title.x=element_text(size=15*ts,
                                                       face='bold'),
                             panel.grid.major=element_blank())
                 )
               )
             )
        p
      })

      output$upsetplot <- renderPlot({
        upsetplot()
      })

      ###################### buttons ###################
      helpButtonServer('upset_controls_help', size='l')
      helpButtonServer('de_upset_help', size='l')
      helpButtonServer('de_upset_tbl_help', size='l')
      downloadButtonServer('upset_download', upsetplot, 'upsetplot')

      ################# get genes from table clicks ############

      # update genes to plot based on table selection
      #observeEvent(c(input$intertable_rows_selected,
      #               input$intertable_row_last_clicked), {
      #  validate(
      #      need(!is.null(intersect_tbl()), 'Waiting for selection')
      #  )
      #  tbl <- intersect_tbl()
      #  sel <- input$intertable_rows_selected
      #  click <- input$intertable_row_last_clicked

      #  # handle NAs in symbol
      #  s <- tbl$symbol
      #  s[is.na(s)] <- tbl$gene[is.na(s)]

      #  gene.to.plot <- plot_args()$gene.to.plot
      #  if(!is.null(click)){
      #      # current genes
      #      g.now <- unique(c( as.character(s)[sel],
      #                        gene.to.plot))

      #      # should a gene be removed?
      #      if(!click %in% sel) g.now <- setdiff(g.now, as.character(s)[click])
      #  } else if(!is.null(gene.to.plot)){
      #      g.now <- gene.to.plot
      #  } else {
      #      g.now <- ''
      #  }

      #  genes_clicked$g <- g.now

      #}) # observeEvent

      upset_label_tbl <- eventReactive(c(plot_args(),
                                         upset_choices$current,
                                         genes_clicked$g,
                                         upset_table$tbl), {
        validate(
          need(!is.null(upset_table$tbl), '')
        )
        g <- genes_clicked$g

        # - DE genes will be found in upset_table
        # - non-DE genes will have 'NA' for set/comparisons
        gdf <- upset_table$tbl

        # check in 'gene' & 'symbol' columns
        idx <- (gdf$symbol %in% g) | (rownames(gdf) %in% g)
        label_tbl <- gdf[idx,]

        # get the non DE genes (if any)
        non_idx <- (!g %in% gdf$symbol) & (!g %in% rownames(gdf))
        if(sum(non_idx) > 0){
          showNotification(
            paste0('Note: Some non-DE genes are not shown in the table - ', paste(g[non_idx], collapse=', ')),
            duration=10
          )
        }

        label_tbl <- cbind(gene=rownames(label_tbl), label_tbl)

        label_tbl %>% select(-any_of(upset_choices$current)) %>%
          relocate('set', .after='comparisons') %>%
          datatable(rownames=FALSE,
                    selection='none',
                    class='stripe')
      })

      output$upset_labeled <- renderDT({
        validate(
          need(length(genes_clicked$g) > 0, '')
        )
        upset_label_tbl()
      })

      # build reactive to return
      # - this is only updated when upset table & genes are changed
      #   *NOT* fc/fdr thres
      upset_data <- eventReactive(c(upset_table$tbl,
                                    upset_table$set_labels,
                                    genes_clicked$g), {
        list(
          tbl=list(all=upset_table$tbl,
                   set_labels=upset_table$set_labels),
          genes=genes_clicked$g
        )
      })

      return(upset_data)

    }
  )
}
