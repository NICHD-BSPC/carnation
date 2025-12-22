#' Heatmap module
#'
#' @description
#' Module UI & server to generate heatmap.
#'
#' @param id Module id
#' @param panel string, can be 'sidebar' or 'main'
#' @param obj reactiveValues object containing carnation object
#' @param coldata reactiveValues object containing object metadata
#' @param plot_args reactive containing 'fdr.thres' (padj threshold), 'fc.thres' (log2FC) &
#' 'upset_data' (list containing data from upset plot module)
#' @param gene_scratchpad reactiveValues object containing genes selected in scratchpad which will
#' be labeled
#' @param config reactive list with config settings
#'
#' @returns
#' UI returns tagList with heatmap UI.
#' Server invisibly returns NULL (used for side effects).
#'
#' @rdname heatmapmod
#' @name heatmapmod
NULL

#' @rdname heatmapmod
#' @export
heatmapUI <- function(id, panel){
  ns <- NS(id)

  # get default config
  config <- get_config()

  if(panel == 'sidebar'){
    tag <-
      tagList(

        fluidRow(
          column(6, align='left',
            tags$label(class='control-label',
                            'Geneset'
            ) # tags$label
          ), # column
          column(6, align='right',
            helpButtonUI(ns('hmap_geneset_help'))
          ) # column
        ),  # fluidRow

        selectInput(ns('geneset_type'), label=NULL,
                    choices=c('de_genes', 'upset_intersections', 'gene_scratchpad')),

        # NOTE: to make conditionalPanel work in modules, we need
        #       to use this hacky workaround
        conditionalPanel(paste0('input["', ns('geneset_type'), '"] == "de_genes"'),
          fluidRow(
            column(5, h5('Comparison')),
            column(7,
              selectInput(ns("de_comp"), label=NULL,
                choices=NULL, selected=NULL
              ) # selectInput
            ) # column
          ) # fluidRow

        ), # conditionalPanel

        conditionalPanel(paste0('input["', ns('geneset_type'), '"] == "upset_intersections"'),
          fluidRow(
            column(12,
              selectizeInput(ns("upset_intersect"), label=NULL,
                             choices=NULL, selected=NULL
              ) # selectInput
            ) # column
          ) # fluidRow

        ), # conditionalPanel


        conditionalPanel(paste0('input["', ns('geneset_type'), '"] == "upset_intersections" | input["', ns('geneset_type'), '"] == "de_genes"'),

          fluidRow(
            column(5, 'Direction of change'),
            column(7,
              selectInput(ns('hmap_type'), label=NULL,
                          choices=c('up & down'='changed',
                                    'up'='upregulated',
                                    'down'='downregulated')
              ) # selectInput
            ) # column
          ) # fluidRow

        ), # conditionalPanel

        fluidRow(
          column(12, strong('Plot options'))
        ), # fluidRow

        fluidRow(
          column(5, h5('sample group')),
          column(7,
            selectInput(ns('hmap_rld'), label=NULL,
                        choices=NULL
            ) # selectInput
          ) # column
        ), # fluidRow

        fluidRow(
          column(5, h5('cluster by')),
          column(7,
            selectInput(ns('hmap_clust'), label=NULL,
                        choices=c('row', 'column', 'both', 'none')
            ) # selectInput
          ) # column
        ), # fluidRow

        fluidRow(
          column(5, h5('scale')),
          column(7,
            selectInput(ns('hmap_scale'), label=NULL,
              choices=c('row', 'none', 'column')
            ) # selectInput
          ) # column
        ), # fluidRow

        fluidRow(
          column(5, '# genes to plot'),
          column(7,
            numericInput(ns('max_gene_num'), label=NULL,
                         value=config$server$de_analysis$heatmap$max_ngenes,
                         min=0, max=500, step=25
            ) # numericInput
          ) # column
        ), # fluidRow

        conditionalPanel(paste0('input["', ns('geneset_type'), '"] == "de_genes"'),
          fluidRow(
            column(5, h5('ranking metric')),
            column(7,
              selectInput(ns('hmap_rank'), label=NULL,
                choices=c('padj', 'log2FoldChange')
              ) # selectInput
            ) # column
          ) # fluidRow
        ), # conditionalPanel


        bsCollapse(
          bsCollapsePanel('column settings',

            fluidRow(
              column(5, h5('labels')),
              column(7,
                selectizeInput(ns('hmap_colnames'), label=NULL,
                               choices=NULL,
                               multiple=TRUE
                ) # selectInput
              ) # column
            ), # fluidRow

            fluidRow(
              column(5, 'group by'),
              column(7,
                selectInput(ns('hmap_cols'),
                            label=NULL,
                            choices=NULL,
                            selected=NULL)
              ) # column
            ), # fluidRow

            selectizeInput(ns('hmap_col_levels'),
                           label=NULL,
                           choices=NULL,
                           selected=NULL,
                           multiple=TRUE),
            actionButton(ns('reset_hmap_cols'),'Reset',
                         class='btn-primary')
          ), # bsCollapsePanel

          bsCollapsePanel('More options',

            fluidRow(
              column(5, h5('row font')),
              column(7,
                numericInput(ns('hmap_fontsize_row'),
                             label=NULL,
                  value=config$ui$de_analysis$heatmap$fontsize_row,
                  step=1
                ) # numericInput
              ) # column
            ), # fluidRow

            fluidRow(
              column(5, h5('column font')),
              column(7,
                numericInput(ns('hmap_fontsize_col'), label=NULL,
                  value=config$ui$de_analysis$heatmap$fontsize_col,
                  step=1
                ) # numericInput
              ) # column
            ), # fluidRow

            fluidRow(
              column(5, h5('colormap')),
              column(7,
                selectInput(ns('hmap_colormap'), label=NULL,
                  choices=c('viridis', 'blue-white-red', 'magma', 'inferno', 'plasma', 'RdBu', 'YlOrRd', 'YlGnBu')
                ) # selectInput
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
    tag <-
      tagList(
        fluidRow(
          column(6, align='left',
            helpButtonUI(ns('de_hmap_help'))
          ), # column
          column(6, align='right',
            downloadButtonUI(ns('hmap_download'))
          ) # column
        ), # fluidRow
        withSpinner(
          plotlyOutput(ns('heatmap'), height='800px')
        ) # withSpinner
      ) # tagList
  }

  tag
}

#' @rdname heatmapmod
#' @export
heatmapServer <- function(id, obj,
                          coldata,
                          plot_args,
                          gene_scratchpad,
                          config){

  moduleServer(
    id,

    function(input, output, session){

      ns <- NS(id)

      coldata.all <- reactive({
        list(all=coldata$all,
             curr=coldata$curr)
      })

      app_object <- reactive({
        list(res=obj$res,
             rld=obj$rld,
             all_rld=obj$all_rld,
             dds_mapping=obj$dds_mapping)
      })

      hmap_coldata <- reactiveValues(all=NULL, current=NULL)

      data_loaded <- reactiveValues(flag=0)

      # reactiveValues to save heatmap data being plotted
      hmap_plot_data <- reactiveValues(all=NULL, plotted=NULL)

      # update from reactive config
      observeEvent(config(), {
        updateNumericInput(session, 'max_ngenes',
                           value=config()$server$de_analysis$heatmap$max_ngenes)
        updateNumericInput(session, 'fontsize_row',
                           value=config()$ui$de_analysis$heatmap$fontsize_row)
        updateNumericInput(session, 'fontsize_col',
                           value=config()$ui$de_analysis$heatmap$fontsize_col)

      })

      # observer to update heatmap comparisons if metadata is changed
      observeEvent(coldata.all()$curr, {
        validate(
          need(!is.null(coldata.all()$curr), 'Waiting for selection')
        )

        updateSelectInput(session, 'hmap_rld',
                          choices=names(coldata.all()$curr))
      }) # observeEvent

      # observer to update heatmap data
      observeEvent(c(app_object()$rld, input$hmap_rld, coldata.all()$curr), {
        validate(
            need(!is.null(app_object()$rld) & !is.null(app_object()$dds_mapping) & input$hmap_rld != '' & input$hmap_rld %in% names(coldata.all()$curr), 'Waiting for selection')
        )

        # get data to plot
        cdata <- coldata.all()$curr[[ input$hmap_rld ]]
        if(input$hmap_rld == 'all_samples'){
            rld <- app_object()$all_rld
            colData(rld) <- cdata
        } else {
            rld <- app_object()$rld[[ input$hmap_rld ]]
            colData(rld) <- cdata
        }

        # update heatmap coldata
        hmap_coldata$all <- cdata
        hmap_coldata$current <- cdata

        column.names <- setdiff(colnames(cdata), config()$server$cols.to.drop)

        if(is.null(input$hmap_cols)) selected <- column.names[1]
        else if(input$hmap_cols %in% column.names) selected <- input$hmap_cols
        else selected <- column.names[1]
        updateSelectInput(session, 'hmap_cols',
                          choices=column.names,
                          selected=selected)

        all_levels <- unique(cdata[, column.names[1]])
        validate(
            need(length(all_levels) > 0, 'Waiting for selection')
        )

        if(all(all_levels %in% input$hmap_col_levels)){
            selected <- input$hmap_col_levels
        } else {
            selected <- all_levels
        }
        updateSelectizeInput(session,
                             'hmap_col_levels',
                             choices=all_levels,
                             selected=selected)

        # set column names
        if(is.null(input$hmap_colnames)) selected <- column.names[1]
        else if(input$hmap_colnames %in% column.names) selected <- input$hmap_cols
        else selected <- column.names[1]
        updateSelectizeInput(session, 'hmap_colnames',
                             choices=column.names,
                             selected=selected)

        # get normalized counts
        all.dat <- assay(rld)

        # check if rownames are NAs and if yes, repl with
        # rownames(rld)
        if(any(is.na(rownames(all.dat))))
            rownames(all.dat) <- rownames(rld)

        hmap_plot_data$all <- all.dat

        data_loaded$flag <- data_loaded$flag + 1
      }) # observeEvent

      # update comparison menu
      observeEvent(app_object()$res, {
        updateSelectInput(session, 'de_comp',
                          choices=names(app_object()$res))
      }) # observeEvent

      thresholds <- reactiveValues(fdr.thres=NULL, fc.thres=NULL)

      # update FDR & FC threshold
      observeEvent(c(plot_args()$fdr.thres, plot_args()$fc.thres), {

        fc.thres <- ifelse(plot_args()$fc.thres == '' | is.na(plot_args()$fc.thres), 0, plot_args()$fc.thres)
        fdr.thres <- ifelse(plot_args()$fdr.thres == '' | is.na(plot_args()$fdr.thres), 0.1, plot_args()$fdr.thres)

        thresholds$fdr.thres <- fdr.thres
        thresholds$fc.thres <- fc.thres
      }) # observeEvent

      upset_data <- reactiveValues(genes=NULL, labels=NULL)
      gene_scratchpad <- reactiveValues(genes=NULL)

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

      observeEvent(gene_scratchpad(), {
        gene_scratchpad$genes <- gene_scratchpad()
      })

      # observer to get heatmap genes
      get_heatmap_genes <- eventReactive(c(input$geneset_type,
                                           input$de_comp,
                                           thresholds$fdr.thres,
                                           thresholds$fc.thres,
                                           input$hmap_type,

                                           upset_data$genes,
                                           input$upset_intersect,

                                           gene_scratchpad$genes,
                                           input$max_gene_num,
                                           input$hmap_rank), {
        validate(
          need(!is.null(app_object()$res), 'Waiting for selection')
        )

        # hard limit on genes plotted
        max_ngenes <- input$max_gene_num

        # scratchpad genes
        scratch <- gene_scratchpad$genes

        if(input$geneset_type == 'de_genes'){
          validate(
            need(!is.null(input$hmap_type), '')
          )

          fdr.thres <- thresholds$fdr.thres
          fc.thres <- thresholds$fc.thres

          # get results df
          res <- app_object()$res[[ input$de_comp ]]

          # order by ranking metric
          if(input$hmap_rank == 'padj'){
              res <- res[order(res[,input$hmap_rank]),]
          } else if(input$hmap_rank == 'log2FoldChange'){
              if(input$hmap_type == 'changed'){
                  res <- res[order(abs(res[, input$hmap_rank]), decreasing=TRUE),]
              } else if(input$hmap_type == 'upregulated'){
                  res <- res[order(res[, input$hmap_rank], decreasing=TRUE),]
              } else if(input$hmap_type == 'downregulated'){
                  res <- res[order(res[, input$hmap_rank]),]
              }
          }

          # order by ranking metric &  get DE gene indices
          idx <- res$padj < fdr.thres & !is.na(res$padj)
          if(input$hmap_type == 'upregulated'){
              idx_fc <- idx & res$log2FoldChange > fc.thres
          } else if(input$hmap_type == 'downregulated'){
              idx_fc <- idx & res$log2FoldChange < -fc.thres
          } else {
              idx_fc <- idx & abs(res$log2FoldChange) > fc.thres
          }

          idx_fc <- which(idx_fc)

          # get gene names
          if('symbol' %in% colnames(res)){
              s <- res$symbol[idx_fc]
          } else if('SYMBOL' %in% colnames(res)){
              s <- res$SYMBOL[idx_fc]
          }

          # if number of genes > max allowed, plot top genes
          if(length(s) > max_ngenes){
            showNotification(
                paste('Too many genes to plot. Showing top',
                      max_ngenes, 'genes based on', input$hmap_rank)
            )
            s <- s[seq_len(max_ngenes)]
          }

        } else if(input$geneset_type == 'upset_intersections'){

          validate(
            need(input$upset_intersect != '',
                 'Please select intersection to label')
          )
          s <- upset_data$genes[[ input$upset_intersect ]]

        } else if(input$geneset_type == 'gene_scratchpad'){
          s <- scratch
        }

        validate(
            need(length(s) > 1, 'Too few genes to plot. Try adjusting filters')
        )

        # return genes
        list(main=s,
             labels=scratch)

      }) # eventReactive get_heatmap_genes

      # observer to update heatmap col levels
      observeEvent(c(input$hmap_cols, hmap_coldata$all, hmap_coldata$current), {
        validate(
            need(!is.null(hmap_coldata$current) & input$hmap_cols != '' & input$hmap_cols %in% colnames(hmap_coldata$current), 'Waiting for selection')
        )

        all_levels <- unique(hmap_coldata$all[, input$hmap_cols])
        curr_levels <- unique(hmap_coldata$current[, input$hmap_cols])

        updateSelectizeInput(session,
                             'hmap_col_levels',
                             choices=all_levels,
                             selected=curr_levels)
      }) # observeEvent

      # observer to update heatmap metadata if col levels are changed
      # NOTE: levels can be removed, reordered or added back
      observeEvent(input$hmap_col_levels, {
        validate(
            need(!is.null(input$hmap_col_levels), 'Waiting for selection')
        )

        # get full metadata and current metadata
        df_all <- hmap_coldata$all
        df_current <- hmap_coldata$current

        # if all levels are present in current metadata
        if(all(input$hmap_col_levels %in% df_current[,input$hmap_cols])){
            # reorder factor order
            df_current[, input$hmap_cols] <- factor(df_current[, input$hmap_cols], levels=input$hmap_col_levels)

            # subset metadata to only keep col levels since
            # some levels in df_current might be missing
            idx <- df_current[,input$hmap_cols] %in% input$hmap_col_levels
            df_current <- df_current[idx,]

            # update coldata reactive
            hmap_coldata$current <- df_current[order(df_current[, input$hmap_cols]),]
        } else if(all(input$hmap_col_levels %in% df_all[,input$hmap_cols])){
            # reorder factor order
            df_all[, input$hmap_cols] <- factor(df_all[, input$hmap_cols], levels=input$hmap_col_levels)

            # get metadata from full corresponding to col levels
            idx <- df_all[,input$hmap_cols] %in% input$hmap_col_levels
            df_all <- df_all[idx,]

            # update coldata reactive
            hmap_coldata$current <- df_all[order(df_all[, input$hmap_cols]),]
        } else if(input$hmap_col_levels == ''){
            hmap_coldata$current <- NULL
        }
      }, ignoreNULL=FALSE) # observeEvent

      # observer to reset metadata
      observeEvent(input$reset_hmap_cols, {
        hmap_coldata$current <- hmap_coldata$all
      })

      # this reactive triggers a redraw when fc/fdr are updated
      # but only when geneset_type is 'de_genes' or 'upset_intersections'
      replot <- reactive({
        # the isolate makes sure choosing different geneset type
        # doesn't trigger redraw, only different fc/fdr thres does
        isolate({
          validate(
            need(input$geneset_type %in% c('de_genes', 'upset_intersections'), ''),
          )
        })
        list(
          plot_args()$fdr.thres,
          plot_args()$fc.thres
        )
      })

      # reactive to plot heatmap
      make_heatmap <- eventReactive(c(replot(),
                                      data_loaded$flag,
                                      gene_scratchpad$genes,
                                      input$plot_do), {

        validate(
          need(!is.null(app_object()$res), 'Waiting for selection')
        )

        validate(
            need(!is.null(hmap_coldata$current), 'Waiting for selection')
        )

        validate(
            need(input$hmap_colnames != '', 'Please select at least one labeling column')
        )

        # get heatmap data
        all.dat <- hmap_plot_data$all
        hmap_genes <- get_heatmap_genes()
        s <- hmap_genes$main

        # get genes that are present
        s <- intersect(s, rownames(all.dat))
        mat <- all.dat[s, ]

        ## if number of genes > max allowed, plot top most variable subset
        ## NOTE: this only happens when geneset is 'upset_intersections'
        if(length(s) > input$max_gene_num){
          showNotification(
              paste('Too many genes to plot. Showing top ',
                    input$max_gene_num, ' most variable genes')
          )

          mat <- mat[order(rowVars(mat), decreasing=TRUE),]
          mat <- mat[seq_len(input$max_gene_num),]
        }

        validate(
            need(!is.null(mat), 'Waiting for selection')
        )

        # subset data based on current coldata
        cols.to.keep <- rownames(hmap_coldata$current)
        if(is.null(cols.to.keep)){
            cols.to.keep <- hmap_coldata$current[, 'samplename']
        }
        mat <- mat[, cols.to.keep]

        hmap_plot_data$plotted <- mat

        # check for and remove duplicate rows
        if(anyDuplicated(rownames(mat)) != 0){
            nodup <- which(duplicated(rownames(mat)) == 0)
            showNotification(
                paste0('Removing ', nrow(mat) - length(nodup), ' duplicated rows')
            )
            mat <- mat[nodup, ]
        }

        # make row_side_colors vector if gene scratchpad is non-empty
        if(any(hmap_genes$labels != '')){
            g <- hmap_genes$labels
            g.keep <- g[g %in% rownames(mat)]
            showNotification(
                paste('Labeling', length(g.keep), 'genes found in current selection')
            )
            # get colors from config
            hmap_colors <- config()$server$de_analysis$heatmap$row_side_colors
            row_side_colors <- rep(hmap_colors[['not_labeled']], nrow(mat))
            row_side_colors[rownames(mat) %in% g.keep] <- hmap_colors[['labeled']]
            row_side_colors <- data.frame('label'=row_side_colors,
                                          row.names=rownames(mat))

            # set 3 subplot widths if row clustering is TRUE
            if(input$hmap_clust %in% c('row', 'both')){
                subplot_widths <- config()$server$de_analysis$heatmap$subplot_widths$with_row_clustering
            } else {
                subplot_widths <- config()$server$de_analysis$heatmap$subplot_widths$without_row_clustering
            }
        } else {
            row_side_colors <- NULL
            subplot_widths <- NULL
        }

        rowv <- ifelse(input$hmap_clust %in% c('row', 'both'),
                       TRUE, FALSE)
        colv <- ifelse(input$hmap_clust %in% c('column', 'both'),
                       TRUE, FALSE)

        if(length(input$hmap_colnames) > 1){
          col_labels <- apply(as.data.frame(hmap_coldata$current[, input$hmap_colnames]),
                              1, function(x) paste(as.character(x), collapse='_'))
        } else {
          col_labels <- hmap_coldata$current[, input$hmap_colnames]
        }
        # Determine colormap based on selection
        colormap <- if(input$hmap_colormap == 'blue-white-red') {
          colorRampPalette(c('blue', 'white', 'red'))(256)
        } else if(input$hmap_colormap %in% c('viridis', 'magma', 'inferno', 'plasma')){
          viridis(256, alpha=1, begin=0, end=1, option=input$hmap_colormap)
        } else if(input$hmap_colormap %in% c('RdBu', 'YlOrRd', 'YlGnBu')){
          brewer.pal(256, input$hmap_colormap)
        }

        p <- tryCatch(
                 heatmaply::heatmaply(mat,
                           Rowv=rowv, Colv=colv,
                           scale=input$hmap_scale,
                           labCol=col_labels,
                           RowSideColors=row_side_colors,
                           subplot_widths=subplot_widths,
                           fontsize_row=input$hmap_fontsize_row,
                           fontsize_col=input$hmap_fontsize_col,
                           colors=colormap),
                 error=function(e){ e })

        # check error message
        if(!heatmaply::is.plotly(p)){
            if(grepl('NA', p$message)){
                e <- 'error: multiple identical rows or rows with all identical values in the data. try setting "scale" to "none"'
            } else {
                e <- paste0('error: ', p$message)
            }
        } else {
            e <- ''
        }

        validate(
            need(is.plotly(p), e)
        )
        p
      }) # eventReactive make_heatmap

      output$heatmap <- renderPlotly({
        make_heatmap()
      })

      ################### buttons #######################

      # heatmap controls help
      helpButtonServer('hmap_geneset_help', size='l')
      helpButtonServer('de_hmap_help', size='l')
      downloadButtonServer('hmap_download', make_heatmap, 'heatmap')

    }
  )
}
