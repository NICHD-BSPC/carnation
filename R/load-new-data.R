#' Load data module UI
#'
#' @param id ID string used to match the ID used to call the module server function
#'
#' @export
loadDataUI <- function(id){
  ns <- NS(id)

  # TODO: add multiple dds
  tagList(
    span('Step 1: ',
         style='font-style: italic;'),
    actionButton(ns('add_counts'), 'Add counts'),
    textOutput(ns('counts_summary')),
    br(), br(),

    span('Step 2: ',
         style='font-style: italic;'),
    actionButton(ns('add_res'), 'Add DE results'),
    textOutput(ns('de_summary')),
    br(), br(),

    span('Step 3: ',
         style='font-style: italic;'),
    actionButton(ns('add_func'), 'Add functional enrichment results'),
    span('(Optional)',
         style='font-style: italic;'),
    textOutput(ns('func_summary')),
    br(), br(),

    actionButton(ns('create_new'), 'Create data set',
                 class='btn-primary')
  ) # tagList
}

#' Load data module server function
#'
#' @param id ID string used to match the ID used to call the module UI function
#' @param username user name
#' @param rds Object to be edited
#'
#' @export
loadDataServer <- function(id, username, rds=NULL){
  moduleServer(
    id,

    function(input, output, session){

      ns <- NS(id)

      config <- get_config()

      # flag to indicate whether to reload parent or not
      reload_parent <- reactiveVal(FALSE)

      # reactiveValues to hold new object
      new_obj <- reactiveValues(res_list=NULL,
                                dds_list=NULL,
                                rld_list=NULL,
                                enrich_list=NULL,
                                degpatterns=NULL,
                                genetonic=NULL)

      #################### edit data ####################

      edit_obj <- reactive({
        rds$obj
      })

      observeEvent(edit_obj(), {
        # match obj names to required patterns
        obj_names <- names(edit_obj())
        dds_name <- setdiff(grep('dds', obj_names),
                            c(grep('all_dds', obj_names), grep('dds_mapping', obj_names)))
        rld_name <- setdiff(grep('rld', obj_names), grep('all_rld', obj_names))
        enrich_name <- grep('enrich', obj_names)
        genetonic_name <- grep('genetonic', obj_names)
        degpatterns_name <-  grep('degpatterns', obj_names)

        # add required elements
        new_obj$res_list <- edit_obj()[[ grep('res', obj_names) ]]
        new_obj$dds_list <- edit_obj()[[ dds_name ]]
        new_obj$rld_list <- edit_obj()[[ rld_name ]]

        # add optional elements
        if(length(enrich_name) > 0) new_obj$enrich_list <- edit_obj()[[ enrich_name ]]
        if(length(genetonic_name) > 0) new_obj$genetonic <- edit_obj()[[ genetonic_name ]]
        if(length(degpatterns_name) > 0) new_obj$degpatterns <- edit_obj()[[ degpatterns_name ]]
      })

      #################### counts & metadata ####################

      observeEvent(input$add_counts, {
        showModal(
          modalDialog(
            fileInput(ns('cnts_file'),
                      label='Counts file'
            ), # fileInput
            span('Tab-delimited txt file (TSV) with gene IDs in the first column.',
                 style='font-style: italic'),
            br(), br(),

            fileInput(ns('cnts_mdata'),
                      label='Metadata file'
            ), # fileInput
            span('Tab-delimited txt file (TSV) with sample names in the first column. Samples should match those in the counts file.',
                 style='font-style: italic;'),
            br(), br(),

            textInput(ns('cnts_name'),
                      label='Counts name',
                      value='main'),
            span('Unique string without white-space or commas to identify counts. Use this to specify sample subsets used in the analysis.',
                 style='font-style: italic;'),
            footer=tagList(
                     actionButton(ns('add_counts_do'), 'OK'),
                     modalButton('Cancel')
                   ),
            easyClose=TRUE
          ) # modalDialog
        ) # showModal
      })

      observeEvent(input$add_counts_do, {

        # check that both counts table & metadata are added
        for(name in c('cnts_file', 'cnts_mdata')){
          if(is.null(input[[ name ]])){
            showNotification(
              'Both counts table & metadata must be added!',
              type='error'
            )

            req(input[[ name ]])
          }
        }

        # check that counts name is not empty
        if(input$cnts_name == ''){
          showNotification(
            'Counts name cannot be empty!',
            type='warning'
          )
        }

        validate(
          need(input$cnts_name != '', '')
        )

        # read counts
        counts <- read.table(input$cnts_file$datapath,
                             sep='\t', header=TRUE)

        # perform checks
        #
        # 1. 1st column of counts (gene IDs) should be character
        # 2. all remaining columns must be numeric
        # 3. 1st column of counts should not have duplicate values
        #
        col1_char <- is.character(counts[,1])
        cols_num <- apply(counts[, 2:ncol(counts)], 2, is.numeric)
        if(!col1_char){
          showNotification(
            'Counts file must contain gene IDs in the first column',
            type='error'
          )

          validate(
            need(col1_char, '')
          )
        } else if(!all(cols_num)){

          showNotification(
            'Counts file must contain gene counts in the 2nd column onwards',
            type='error'
          )

          validate(
            need(all(cols_num), '')
          )
        } else {
          # check that first column does not contain duplicate values
          if(sum(duplicated(counts[,1])) > 0){
            showNotification(
              'Duplicate gene IDs in 1st column!',
              type='error'
            )

            validate(
              need(sum(duplicated(counts[,1])) == 0, '')
            )
          }

          # set rownames & drop 1st column
          rownames(counts) <- counts[,1]
          counts <- counts[, 2:ncol(counts)]
        }

        # read metadata
        coldata <- read.table(input$cnts_mdata$datapath,
                              sep='\t', header=TRUE)

        # perform checks
        # 1st column (samplenames) should not have duplicates
        if(sum(duplicated(coldata[,1])) > 0){
          showNotification(
            'Duplicate values in 1st column!',
            type='error'
          )

          validate(
            need(sum(duplicated(coldata[,1])) == 0, '')
          )
        }

        # set rownames
        rownames(coldata) <- coldata[,1]

        # perform some checks
        # 1. ncol counts should be == nrow metadata
        # 2. all rownames of coldata must be present in colnames of counts
        if(nrow(coldata) != ncol(counts)){
          showNotification(
            paste0('# of samples in counts (n = ', ncol(counts),
                   ') and metadata files (n = ', nrow(coldata), ') must match!'),
            type='error'
          )

          validate(
            need(nrow(coldata) == ncol(counts), '')
          )
        } else if(!all(rownames(coldata) %in% colnames(counts)) | !all(colnames(counts) %in% rownames(coldata))){
          showNotification(
            'Sample names in metadata and counts file must match. Please retry',
            type='error'
          )

          validate(
            need(all(rownames(coldata) %in% colnames(counts)) &
                 all(colnames(counts) %in% rownames(coldata)), '')
          )
        } else {
          showModal(
            modalDialog(
              'Building count object & normalizing data',
              footer=NULL
            )
          )

          # order counts columns by coldata rownames
          counts <- counts[, rownames(coldata)]

          dds <- DESeqDataSetFromMatrix(counts,
                                        colData=coldata,
                                        design=~1)

          rld <- varianceStabilizingTransformation(dds, blind=TRUE)

          if(is.null(new_obj$dds_list)){
            new_obj$dds_list <- setNames(list(dds), input$cnts_name)
            new_obj$rld_list <- setNames(list(rld), input$cnts_name)
          } else {
            new_obj$dds_list[[ input$cnts_name ]] <- dds
            new_obj$rld_list[[ input$cnts_name ]] <- rld
          }

        }
        removeModal()
      })

      output$counts_summary <- renderText({
        if(!is.null(new_obj$dds_list)){
            nsets <- length(new_obj$dds_list)
            nsamples <- unlist(lapply(new_obj$dds_list, function(x){
                            nrow(colData(x))
            }))

            if(nsets > 1){
              msg <- paste(nsets, 'sample groups:',
                           paste(nsamples, collapse=', '),
                           'samples, respectively')
            } else {
              msg <- paste(nsets, 'sample group:',
                           nsamples, 'samples')
            }
            msg
        }
      })

      #################### DE results ####################

      observeEvent(input$add_res, {
        if(is.null(new_obj$dds_list)){
          showNotification(
            'Must add at least one counts table before uploading DE results!',
            type='error'
          )
        }

        validate(
          need(!is.null(new_obj$dds_list), 'Must have counts')
        )

        counts_choices <- names(new_obj$dds_list)
        showModal(
          modalDialog(
            fileInput(ns('res_file'),
                      label='DE results file(s)',
                      multiple=TRUE),
            tags$div(
              span('*Can add multiple files here'),
              br(), br(),

              span('These should be tab-delimited text files (TSV) generated by DE analysis tools:'),
              tags$ul(
                tags$li('Gene IDs should be in a "gene" column.'),
                tags$li('Gene symbols can also be present in a column named "symbol".'),
                tags$li('Only "DESeq2" results are supported at the moment.')
              ),
              style='font-style: italic;'
            ), # div

            footer=tagList(
                     actionButton(ns('add_res_files'), 'OK'),
                     modalButton('Cancel')
                   ),
            easyClose=TRUE
          ) # modalDialog
        ) # showModal
      })

      observeEvent(input$add_res_files, {

        if(is.null(input$res_file)){
          showNotification(
            'No DE results uploaded!',
            type='error'
          )
        }

        req(input$res_file)

        # check that at least 1 DE result is uploaded
        counts_choices <- names(new_obj$dds_list)
        tag <- tagList(
                 fluidRow(
                   column(3, strong('Name')),
                   column(4, strong('Description')),
                   column(3, strong('Count dataset'))
                 ) # fluidRow
               ) # tagList

        for(i in 1:nrow(input$res_file)){
          # get placeholder
          tmp_id <- tools::file_path_sans_ext(basename(input$res_file$name[i]))

          tag <- tagAppendChildren(
                   tag,
                   fluidRow(
                     column(3,
                       textInput(ns(paste0('res_id', i)),
                                 label=NULL,
                                 value=tmp_id)
                     ), # column

                     column(4,
                       textInput(ns(paste0('res_label', i)),
                                 label=NULL,
                                 value=tmp_id)
                     ), # column

                     column(3,
                       selectInput(ns(paste0('res_counts', i)),
                                   label=NULL,
                                   choices=counts_choices,
                                   selected=NULL)
                     ) # column
                   )
                 )
        }

        tag <- tagAppendChildren(
                 tag,
                 tags$div(
                   tags$ul(
                     tags$li('Name: Unique name for comparison without white-space or commas'),
                     tags$li('Description: Short description of comparison'),
                     tags$li('Count dataset: Counts data corresponding to results')),
                 style='font-style: italic;')
               )

        # show modal for multiple DE files
        showModal(
          modalDialog(
            tag,
            footer=tagList(
                     actionButton(ns('add_res_do'), 'OK'),
                     modalButton('Cancel')
                   ),
            easyClose=TRUE
          ) # modalDialog
        ) # showModal
      })

      observeEvent(input$add_res_do, {
        req(input$res_file)

        # get number of uploaded DE results
        nres <- nrow(input$res_file)
        all_names <- paste0(c('res_id', 'res_label'), 1:nres)

        # check that counts name & label are not empty
        for(name in c(all_names)){
          if(input[[ name ]] == ''){
            showNotification(
              'DE comparison name or description cannot be empty!',
              type='warning'
            )
          }

          validate(
            need(input[[ name ]] != '', '')
          )
        }

        # read DE files & build res_list
        for(i in 1:nres){
          res <- read.table(input$res_file$datapath[i],
                            sep='\t', header=TRUE)

          res_id <- input[[ paste0('res_id', i) ]]
          res_counts <- input[[ paste0('res_counts', i) ]]
          res_label <- input[[ paste0('res_label', i) ]]

          # check for DESeq2 columns
          deseq2_cols <- c('baseMean', 'log2FoldChange', 'padj', 'pvalue')

          if(!all(deseq2_cols  %in% colnames(res))){
            missing <- setdiff(deseq2_cols, colnames(res))
            showNotification(
              paste0('DE results table does not match DESeq2 format. Missing columns: ',
                     paste(missing, collapse=',')),
              type='error'
            )

            validate(
              need(all(deseq2_cols %in% colnames(res)), '')
            )
          }

          # check for 'gene' column
          if(!'gene' %in% tolower(colnames(res))){
            showNotification(
              'DE results table must contain "gene" column!',
              type='error'
            )

            validate(
              need('gene' %in% tolower(colnames(res)), '')
            )
          } else {
            # set gene column name to exactly 'gene'
            gcol <- grep('gene', tolower(colnames(res)))
            colnames(res)[gcol] <- 'gene'

            rownames(res) <- res$gene
            res <- res %>% relocate('gene')
          }

          # check for 'symbol' column
          # - if absent, copy 'gene' column instead
          if(!'symbol' %in% tolower(colnames(res))){
            showNotification(
              'DE results table does not contain "symbol" column! Using "gene" column instead',
              type='warning'
            )
            res[[ 'symbol' ]] <- res[[ 'gene' ]]
          } else {
            # set symbol column name to exactly 'symbol'
            gcol <- grep('symbol', tolower(colnames(res)))
            colnames(res)[gcol] <- 'symbol'
          }

          res_list <- list(res=res,
                           dds=res_counts,
                           label=res_label)

          if(is.null(new_obj$res_list)){
            new_obj$res_list <- list(res_list)
            names(new_obj$res_list) <- res_id
          } else {

            if(res_id %in% names(new_obj$res_list)){
              showNotification(
                paste0('DE results named "', res_id,
                       '" already present in object. Please choose different name'),
                type='error'
              )

              validate(
                need(!res_id %in% names(new_obj$res_list), '')
              )
            } else {
              new_obj$res_list[[ res_id ]] <- res_list
            }
          }
        }
        removeModal()
      })

      output$de_summary <- renderText({
        if(!is.null(new_obj$res_list)){
            nsets <- length(new_obj$res_list)

            if(nsets > 1){
              msg <- paste(nsets, 'DE results')
            } else {
              msg <- paste(nsets, 'DE result')
            }
            msg
        }
      })

      #################### FE results ####################

      observeEvent(input$add_func, {
        if(is.null(new_obj$dds_list)){
          showNotification(
            'Must add at least one counts table before uploading FE results!',
            type='error'
          )
        }

        validate(
          need(!is.null(new_obj$dds_list), 'Must have counts')
        )

        if(is.null(new_obj$res_list)){
          showNotification(
            'Must add at least one DE results table before uploading FE results!',
            type='error'
          )
        }

        validate(
          need(!is.null(new_obj$res_list), 'Must have DE results')
        )

        res_choices <- c('Choose one'='', names(new_obj$res_list))
        if(!is.null(new_obj$enrich_list)){
          func_choices <- c('Choose one or enter new'='',
                            names(new_obj$enrich_list))
        } else {
          func_choices <- ''
        }

        pathway_choices <- unlist(config$server$functional_enrichment$pathways)

        showModal(
          modalDialog(
            fileInput(ns('func_file'),
                      label='Functional enrichment results file(s)',
                      multiple=TRUE),

            tags$div(
              span('*Can add multiple files here'),
              br(), br(),

              span('These should be tab-delimited text files (TSV) generated by functional enrichment analysis tools.'),
              tags$ul(
                tags$li('Only "clusterProfiler" output is supported at the moment.')
              ),
              style='font-style: italic;'
            ), # div

            footer=tagList(
                     actionButton(ns('add_func_files'), 'OK'),
                     modalButton('Cancel')
                   ),
            easyClose=TRUE
          ) # modalDialog
        ) # showModal
      })


      observeEvent(input$add_func_files, {
        # check that at least 1 FE result is uploaded
        if(is.null(input$func_file)){
          showNotification(
            'No FE results uploaded!',
            type='error'
          )
        }

        req(input$func_file)

        res_choices <- c('Choose one'='', names(new_obj$res_list))
        if(!is.null(new_obj$enrich_list)){
          func_choices <- c('Choose one or enter new'='',
                            names(new_obj$enrich_list))
        } else {
          func_choices <- ''
        }

        pathway_choices <- unlist(config$server$functional_enrichment$pathways)

        tag <- tagList(
                 fluidRow(
                   column(3, strong('Name')),
                   column(3, strong('Comparison name')),
                   column(3, strong('Effect class')),
                   column(3, strong('Pathway'))
                 ) # fluidRow
               ) # tagList

        for(i in 1:nrow(input$func_file)){

          # get placeholder names
          tmp_id <- tools::file_path_sans_ext(basename(input$func_file$name[i]))

          tag <- tagAppendChildren(
                   tag,
                   fluidRow(
                     column(3,
                       textInput(ns(paste0('func_id', i)),
                                 label=NULL,
                                 value=tmp_id)
                     ), # column

                     column(3,
                       selectizeInput(ns(paste0('func_res_id', i)),
                                      label=NULL,
                                      choices=res_choices)
                     ), # column

                     column(3,
                       selectInput(ns(paste0('func_effect', i)),
                                   label=NULL,
                                   choices=c('changed', 'up', 'down', 'none'),
                                   selected=NULL)
                     ), # column

                     column(3,
                       selectizeInput(ns(paste0('func_pathway', i)),
                                      label=NULL,
                                      choices=pathway_choices,
                                      options=list(create=TRUE))
                     ) # column
                   )
                 )
        }

        tag <- tagAppendChildren(
                 tag,
                 tags$div(
                   tags$ul(
                     tags$li('Name: Unique name for FE result without white-space or commas.'),
                     tags$li('Comparison name: Name of DE result corresponding to FE result.'),
                     tags$li('Effect class: Direction of change. This can be "changed" (up & down), "up", "down" or "none" for gene lists made using other criteria, e.g. upset plot intersections.'),
                     tags$li('Pathway: Name of pathway/ontology. Choose from the options or create new.')),
                 style='font-style: italic;')
               )

        # show modal for multiple DE files
        showModal(
          modalDialog(
            tag,
            footer=tagList(
                     actionButton(ns('add_func_do'), 'OK'),
                     modalButton('Cancel')
                   ),
            size='l',
            easyClose=TRUE
          ) # modalDialog
        ) # showModal

      })

      output$func_summary <- renderText({
        if(!is.null(new_obj$enrich_list)){
          nsets <- 0
          for(fid in names(new_obj$enrich_list)){
            for(d in names(new_obj$enrich_list[[ fid ]])){
              for(ont in names(new_obj$enrich_list[[ fid ]][[ d ]])){
                nsets <- nsets + 1
              }
            }
          }

          if(nsets > 1){
            msg <- paste(nsets, 'FE results')
          } else {
            msg <- paste(nsets, 'FE result')
          }
          msg
        }
      })

      observeEvent(input$add_func_do, {
        req(input$func_file)

        neres <- nrow(input$func_file)
        all_names <- paste0(c('func_id', 'func_res_id', 'func_pathway'), 1:neres)

        # check that FE name, DE comp name & pathway names are not empty
        for(name in all_names){
          if(input[[ name ]] == ''){
            showNotification(
              'FE results, DE comparison or pathway names cannot be empty!',
              type='warning'
            )
          }

          validate(
            need(input[[ name ]] != '', '')
          )
        }

        for(i in 1:neres){
          # read FE table
          eres <- read.table(input$func_file$datapath[i],
                            sep='\t', header=TRUE)

          func_id <- input[[ paste0('func_id', i) ]]
          func_res_id <- input[[ paste0('func_res_id', i) ]]
          func_effect <- input[[ paste0('func_effect', i) ]]
          func_pathway <- input[[ paste0('func_pathway', i) ]]

          # TODO: move to config?
          # OR columns
          cprof_or_cols <- c('ID', 'Description', 'GeneRatio', 'BgRatio',
                             'pvalue', 'p.adjust', 'qvalue', 'geneID', 'Count')

          # GSEA columns
          cprof_gsea_cols <- c('ID', 'Description', 'core_enrichment', 'setSize',
                               'pvalue', 'p,adjust', 'qvalue', 'NES', 'setSize')

          # check for clusterProfiler columns
          if(all(cprof_or_cols %in% colnames(eres))){
            rownames(eres) <- eres$ID
            e <- makeEnrichResult(eres)
          } else if(all(cprof_gsea_cols) %in% colnames(eres)){
            rownames(eres) <- eres$ID
            e <- makeEnrichResult(eres, type='gseaResult')
          } else {
            showNotification(
              'FE results table does not match clusterProfiler format!',
              type='error'
            )

            validate(
              need(all(cprof_or_cols %in% colnames(eres)) |
                   all(cprof_gsea_cols %in% colnames(eres)), '')
            )
          }

          # get de results
          res <- new_obj$res_list[[ func_res_id ]][[ 'res' ]]

          # make genetonic object
          gt <- enrich_to_genetonic(e, res)

          e_list <- setNames(list(eres), func_pathway)
          gt_list <- setNames(list(gt), func_pathway)

          if(is.null(new_obj$enrich_list)){
            new_obj$enrich_list <-
              setNames(
                list(setNames(list(e_list), func_effect)),
                func_id
              )
            new_obj$genetonic <-
              setNames(
                list(setNames(list(gt_list), func_effect)),
                func_id
              )

            # save res_list key under 'res'
            new_obj$enrich_list[[ func_id ]][[ 'res' ]] <- func_res_id
          } else if(!func_id %in% names(new_obj$enrich_list)){
            new_obj$enrich_list[[ func_id ]] <-
              setNames(list(e_list), func_effect)
            new_obj$genetonic[[ func_id ]] <-
              setNames(list(gt_list), func_effect)

            # save res_list key under 'res'
            new_obj$enrich_list[[ func_id ]][[ 'res' ]] <- func_res_id
          } else if(!func_effect %in% names(new_obj$enrich_list[[ func_id ]])){
            new_obj$enrich_list[[ func_id ]][[ func_effect ]] <- e_list
            new_obj$genetonic[[ func_id ]][[ func_effect ]] <- gt_list
          } else {
            showNotification(
              paste0('FE result already exists: ',
                     'Name: "', func_id, '", ',
                     'Effect class: "', func_effect, '", ',
                     'Pathway: "', func_pathway, '"'),
              type='error'
            )

            validate(
              need(!(func_id %in% names(new_obj$enrich_list) &
                     func_effect %in% names(new_obj$enrich_list[[ func_id ]]) &
                     func_pathway %in% names(new_obj$enrich_list[[ func_id ]][[ func_effect ]])), '')
            )
          }
        }
        removeModal()
      })

      #################### save object ####################

      observeEvent(input$create_new, {

        if(is.null(new_obj$dds_list) & is.null(new_obj$res_list)){
          showNotification(
            "Must add at least 1 counts table & DE results table!",
            type='warning'
          )
        }

        validate(
          need(!is.null(new_obj$dds_list) & !is.null(new_obj$res_list), '')
        )

        showModal(
          modalDialog(
            'Saving object ...',
            footer=NULL
          )
        )

        showModal(
          obj_dir_modal(username())
        )
      })

      observeEvent(input$dir_new_do, {
        if(input$dir_new == ''){
          showNotification(
            'Please select or enter new data before saving file',
            type='warning'
          )
        }

        validate(
          need(input$dir_new != '', 'Waiting for input')
        )

        odir <- file.path(input$dir_new, input$proj_new)
        msg <- tryCatch(
                 dir.exists(odir),
                 error=function(e){ e }
               )

        if(!is.logical(msg)){
          showModal(
            modalDialog(
              span('Error saving to directory. Please retry',
                   style='color: red;'),
              footer=NULL
            )
          )

          Sys.sleep(3)
          showModal(
            obj_dir_modal(username())
          )
        } else {
          msg <- tryCatch(
                    dir.create(odir, recursive=TRUE),
                    warning=function(w){ w },
                    error=function(e){ e }
                 )

          if(!is.logical(msg)){
            if(!grepl('already exists', msg$message)){
              showNotification(
                paste('Error creating directory: ', msg$message,
                       'Please retry'),
                type='error'
              )
            } else {
              dir_ok <- TRUE
            }
          } else {
            dir_ok <- TRUE
          }

          validate(
            need(dir_ok, '')
          )

          # TODO: get suffix 'rnaseq.rds' from pattern
          ofile <- file.path(odir,
                             paste0(input$analysis_label,
                                   '.rds')
                            )

          if(!file.exists(ofile) | as.logical(input$overwrite)){
            # save file, add data area, reload
            showModal(
              modalDialog(
                'Saving RDS file',
                footer=NULL
              )
            )

            # create carnation obj & save
            combined <- list(res_list=new_obj$res_list,
                             dds_list=new_obj$dds_list,
                             rld_list=new_obj$rld_list,
                             enrich_list=new_obj$enrich_list,
                             genetonic=new_obj$genetonic,
                             degpatterns=new_obj$degpatterns)

            combined_final <- make_final_object(combined)

            # NOTE: remove .Environment attributes from @design slots of obj$dds
            #   elements & obj$all_dds
            # - this prevents saved object from becoming very large if another
            #   object has been previously loaded
            combined_final$dds <- lapply(combined_final$dds, function(x){
                                    attr(x@design, '.Environment') <- NULL
                                    x
                                  })

            attr(combined_final$all_dds@design, '.Environment') <- NULL

            saveRDS(combined_final, ofile,
                    compress=as.logical(input$compress))

            # add data area
            # get access yaml and add data area
            y <- read_access_yaml()
            if(is.null(username())) ug <- config$server$admin_group
            else ug <- input$user_group

            # check for empty user group
            if(ug == ''){
              showNotification(
                paste0('User group is empty. Using "', config$server$admin_group,
                       '" by default'),
                type='warning'
              )
              ug <- config$server$admin_group
            }

            # check for existence of user_group & add if new
            if(!ug %in% names(y$data_area)){
              y$data_area <- append(y$data_area,
                                setNames(as.list(input$dir_new), ug))
              showModal(
                modalDialog('Adding data area')
              )
            } else {
              if(!input$dir_new %in% y$data_area[[ug]]){
                y$data_area[[ug]] <- c(y$data_area[[ug]],
                                       input$dir_new)
              }
            }
            save_access_yaml(y)

            # reload
            showModal(
              modalDialog(
                'Reloading the app. Please wait',
                footer=NULL
              )
            )
            reload_parent(TRUE)
          } else {
            showNotification(
              'File already exists! Please retry with different analysis label or set "Force overwrite?" to TRUE',
              type='warning'
            )

            validate(
              need(!file.exists(ofile), '')
            )
          }
        }
      })

      obj_dir_modal <- function(username){
        # get data areas for user
        y <- read_access_yaml()

        # if username is NULL (single-user mode), user default user
        if(is.null(username)) username <- config$server$default_user
        ug <- y$user_group[[ username ]]
        da <- y$data_area[[ ug ]]

        tags <- tagList(
                  selectizeInput(ns('dir_new'), 'Choose data area to save object',
                                 choices=c('Select or create new'='',
                                           unname(unlist(da))),
                                 width='100%',
                                 options=list(create=TRUE)),

                  textInput(ns('proj_new'), 'Project name',
                            value='my-project'),
                  span('This will show up in "Available projects" menu',
                       style='font-style: italic;'), br(), br(),

                  textInput(ns('analysis_label'), 'Analysis label',
                            value='main'),
                  span('This will show up in "Available analyses" menu',
                       style='font-style: italic;'), br(), br(),
                )

        # don't show user group input in single-user mode
        if(username != config$server$default_user){
          tags <- tagAppendChildren(
                    tags,
                    textInput(ns('user_group'), 'User group',
                              value=config$server$group_admin),
                    span('User group that will have access to this data',
                         style='font-style: italic;'), br()
                  )
        }

        modalDialog(
          tags,
          selectInput(ns('compress'), label='Compress RDS?',
                      choices=c(FALSE, TRUE)),

          span('Uncompressed files load faster, but are larger on disk.',
               style='font-style: italic;'), br(), br(),

          selectInput(ns('overwrite'), label='Force overwrite?',
                      choices=c(FALSE, TRUE)),
          span('Caution: Setting this to TRUE will overwrite RDS file on disk',
               style='font-style: italic; color: red;'), br(), br(),

          footer=tagList(
                   actionButton(ns('dir_new_do'), 'OK'),
                   modalButton('Cancel')
                 )
        )
      }


      return(reactive(reload_parent()))

    } # function
  ) # moduleServer
} # function

