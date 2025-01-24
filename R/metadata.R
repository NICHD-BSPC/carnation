#' Metadata module ui
#'
#' This generates the metadata tab that allows users to
#' view the metadata attached to the objects being viewed.
#'
#' @param id Input id
#' @param panel context for generating ui elements ('sidebar' or 'main')
#'
metadataUI <- function(id, panel){
    ns <- NS(id)

    if(panel == 'sidebar'){
        tag <- tagList(
                 selectizeInput(ns('meta_key'),
                                label=NULL,
                                choices=NULL,
                                selected=NULL),

                 fluidRow(
                   column(6, strong('Edit metadata')),
                   column(6, align='right',
                     helpButtonUI(ns('edit_metadata_help'))
                   ) # column
                 ), # fluidRow

                 fluidRow(
                   column(12, align='center',
                     actionButton(ns('meta_add_dup'),
                       'Duplicate column',
                       style='margin-bottom: 10px')
                   ),
                   column(12, align='center',
                     actionButton(ns('meta_add_empty'),
                       'Add empty column',
                       style='margin-bottom: 10px')
                   ),

                   column(12, align='center',
                     actionButton(ns('meta_rm'),
                       'Remove column',
                       style='margin-bottom: 10px')
                   ) # column
                 ), # fluidRow

                 fluidRow(
                   column(6, align='right',
                       actionButton(ns('meta_apply'), 'Apply changes',
                             class='btn-primary',
                             style='margin-top: 10px;')
                   ),
                   column(6, align='center',
                       actionButton(ns('meta_reset'), 'Reset',
                             style='margin-top: 10px;')
                   ) # column
                 ) # fluidRow
               )  # tagList
    } else if(panel == 'main'){
        tag <- tagList(
                    fluidRow(
                       column(10,
                            withSpinner(
                                DTOutput(ns('metadata'))
                            ) # withSpinner
                        ) # column
                    ) # fluidRow
                ) # tagList
    }

    return(tag)
}

#' Metadata module server
#'
#' Server code for settings module
#'
#' @param id Input id
#' @param obj internal app object
#' @param cols.to.drop columns to hide from table
#'
metadataServer <- function(id, obj, cols.to.drop){
    moduleServer(
        id,

        function(input, output, session){

            ns <- NS(id)

            # reactive values to hold metadata
            coldata_all <- reactiveValues(init=NULL,
                                          curr=NULL,
                                          staging=NULL)

            # set/reset all metadata
            observeEvent(obj$dds, {
              validate(
                  need(!is.null(obj$dds) & !is.null(obj$dds_mapping),
                       'Waiting for data')
              )

              if(!is.null(obj$all_dds)){
                  all.cdata <- colData(obj$all_dds)
                  if('samplename' %in% colnames(all.cdata)){
                      all.list <- list('all_samples'=all.cdata)
                  } else {
                      rnames <- rownames(all.cdata)
                      all.df <- cbind(samplename=rnames, all.cdata)
                      rownames(all.df) <- rnames
                      all.list <- list('all_samples'=all.df)
                  }
                  coldata_all$init <- all.list
              } else {
                  coldata_all$init <- NULL
              }

              if(length(obj$dds) > 1){
                  clist <- lapply(obj$dds, function(x){
                                      df <- colData(x)

                                      # if 'samplename' column is missing make one with rownames
                                      if(!'samplename' %in% colnames(df)){
                                          rnames <- rownames(df)
                                          df <- cbind(samplename=rownames(df), df)
                                          rownames(df) <- rnames
                                      }
                                      df
                                  })

                  if(!is.null(coldata_all$init)){
                      coldata_all$init <- append(coldata_all$init, clist)
                  } else {
                      coldata_all$init <- clist
                  }
              }

              # copy to staging & current area
              coldata_all$curr <- coldata_all$init
              coldata_all$staging <- coldata_all$init

              # update input
              updateSelectInput(session, 'meta_key',
                                choices=names(coldata_all$init))
            })

            observeEvent(input$meta_add_dup, {
              validate(
                  need(!is.null(coldata_all$staging) & input$meta_key %in% names(coldata_all$staging), 'Waiting for data')
              )

              df <- coldata_all$staging[[ input$meta_key ]]

              cols.to.keep <- colnames(df)[!colnames(df) %in% c('samplename', cols.to.drop)]
              if(length(cols.to.keep) > 1){
                  df <- df[, cols.to.keep]
              }

              showModal(
                  modalDialog(
                      title='Duplicate column',
                      fluidRow(
                          column(6, 'Select column'),
                          column(6,
                              selectInput(ns('meta_dup_col'), label=NULL,
                                      choices=c('', colnames(df))))),
                      fluidRow(
                          column(6, 'Column label'),
                          column(6,
                              textInput(ns('meta_dup_col_name'),
                                        label=NULL, value=''))),
                      span('If left empty, column will be named "<selected column>_dup"', style='font-style: italic;'),
                      footer=tagList(
                          modalButton('Cancel'),
                          actionButton(ns('meta_do_dup'), 'OK')
                        ),
                      easyClose=TRUE
                  )
              )
            })

            observeEvent(input$meta_rm, {
              validate(
                  need(!is.null(coldata_all$staging) & input$meta_key %in% names(coldata_all$staging), 'Waiting for data')
              )

              column_choices <- setdiff(colnames(coldata_all$staging[[ input$meta_key ]]), c(cols.to.drop, 'samplename'))

              showModal(
                  modalDialog(
                      title='Delete column',
                      fluidRow(
                          column(6,  'Select column(s)'),
                          column(6,
                              selectizeInput(ns('meta_rm_col'), label=NULL,
                              choices=c('', column_choices),
                              multiple=TRUE))),
                      span('You can choose multiple columns to delete',
                           style='font-style: italic;'),
                      footer=tagList(
                          modalButton('Cancel'),
                          actionButton(ns('meta_do_rm'), 'OK')
                        ),
                      easyClose=TRUE
                  )
              )
            })

            colexistsModal <- function(){
              modalDialog(
                  div(tags$b('Column name already exists! Please choose different name for new column',
                       style='color: red;')),
                  footer=tagList(
                              modalButton('OK')
                              ),
                  easyClose=TRUE
              )
            }

            # duplicate selected column
            observeEvent(input$meta_do_dup, {
              validate(
                  need(!is.null(coldata_all$staging) & input$meta_key %in% names(coldata_all$staging) & input$meta_dup_col != '', 'Waiting for data')
              )

              df <- coldata_all$staging[[ input$meta_key ]]

              if(input$meta_dup_col_name == ''){
                  # add new user column
                  user_col <- paste0(input$meta_dup_col, '_dup')
              } else {
                  user_col <- input$meta_dup_col_name
              }

              if(user_col %in% colnames(df)){
                  showModal(colexistsModal())
              } else {
                  df[, user_col] <- as.character(df[, input$meta_dup_col])

                  coldata_all$staging[[ input$meta_key ]] <- df

                  removeModal()
              }
            })

            # add empty column
            observeEvent(input$meta_add_empty, {
              validate(
                  need(!is.null(coldata_all$staging) & input$meta_key %in% names(coldata_all$staging), 'Waiting for data')
              )

              df <- coldata_all$staging[[ input$meta_key ]]

              # get number of columns with prefix 'user_col'
              user_colnum <- length(grep('^user_col', colnames(df)))

              # new default empty column name
              user_col <- paste0('user_col', user_colnum + 1)

              showModal(
                  modalDialog(
                      title='Add empty column',
                      fluidRow(
                          column(6, 'Column label'),
                          column(6,
                              textInput(ns('meta_empty_col_name'),
                                        label=NULL, value=''))),
                      span(paste0('If left empty, column will be named "',
                                  user_col, '"'), style='font-style: italic;'),
                      footer=tagList(
                          modalButton('Cancel'),
                          actionButton(ns('meta_do_empty'), 'OK')
                        ),
                      easyClose=TRUE
                  )
              )
            })

            observeEvent(input$meta_do_empty, {

              df <- coldata_all$staging[[ input$meta_key ]]

              if(input$meta_empty_col_name == ''){
                  # get number of columns with prefix 'user_col'
                  user_colnum <- length(grep('^user_col', colnames(df)))

                  # add new user column
                  user_col <- paste0('user_col', user_colnum + 1)
              } else {
                  user_col <- input$meta_empty_col_name
              }

              if(user_col %in% colnames(df)){
                  showModal(colexistsModal())
              } else {
                  df[, user_col] <- rep('', nrow(df))

                  coldata_all$staging[[ input$meta_key ]] <- df

                  removeModal()
              }
            })

            # remove selected column
            observeEvent(input$meta_do_rm, {
              validate(
                  need(!is.null(coldata_all$staging) & input$meta_key %in% names(coldata_all$staging) & input$meta_rm_col != '', 'Waiting for data')
              )

              df <- coldata_all$staging[[ input$meta_key ]]

              df <- df[, setdiff(colnames(df), input$meta_rm_col)]
              coldata_all$staging[[ input$meta_key ]] <- df

              removeModal()

            })

            # get metadata table
            get_metadata_tbl <- eventReactive(c(coldata_all$staging, input$meta_key), {
              validate(
                  need(!is.null(input$meta_key), 'Waiting for data')
              )

              df <- coldata_all$staging[[ input$meta_key ]]
              cols.to.keep <- setdiff(colnames(df), cols.to.drop)
              if(length(cols.to.keep) > 1){
                  df <- df[, cols.to.keep]
              }

              # get user defined columns
              user_cols <- setdiff(colnames(df),
                                   colnames(coldata_all$init[[ input$meta_key ]]))

              # non-user defined columns are fixed
              fixed_cols <- setdiff(colnames(coldata_all$staging[[ input$meta_key]]),
                                    user_cols)
              fixed_cols_idx <- which(colnames(df) %in% fixed_cols)

              as.data.frame(df) %>%
                  datatable(editable=list(target='column',
                                          disable=list(columns=c(0, fixed_cols_idx))))
            })

            output$metadata <- renderDT({
              get_metadata_tbl()
            })

            # reactive to edit metadata
            observeEvent(input$metadata_cell_edit, {
              validate(
                  need(!is.null(coldata_all$staging) & input$meta_key %in% names(coldata_all$staging), 'Waiting for data')
              )

              # edited coldata df
              edit_df <- input$metadata_cell_edit

              # which cols have edits
              edit_cols <- unique(edit_df$col)

              # staged coldata df
              meta_df <- coldata_all$staging[[ input$meta_key ]]
              cols.to.keep <- setdiff(colnames(meta_df), cols.to.drop)
              meta_df <- meta_df[, cols.to.keep]

              # update values in staged df
              for(i in 1:nrow(edit_df)){
                  r <- edit_df$row[i]
                  c <- edit_df$col[i]
                  v <- edit_df$value[i]
                  meta_df[r, c] <- v
              }

              coldata_all$staging[[ input$meta_key ]] <- meta_df
            })

            # save edited metadata
            observeEvent(input$meta_apply, {
              showModal(
                  modalDialog(
                      div(tags$b('This will apply any changes to the metadata to the loaded dataset. Proceed?', style='color: red;')),
                      footer=tagList(
                          modalButton('Cancel'),
                          actionButton(ns('meta_apply_ok'), 'Yes'))
                      )
                  )
            })

            observeEvent(input$meta_apply_ok, {
              coldata_all$curr <- coldata_all$staging
              removeModal()
            })

            # reset metadata table
            observeEvent(input$meta_reset, {

              # get user columns
              all_user_cols <- setdiff(colnames(coldata_all$staging[[ input$meta_key ]]),
                                       colnames(coldata_all$init[[ input$meta_key ]]))

              showModal(
                  modalDialog(
                      div(tags$b('This will reset any changes made to the metadata. Proceed?', style='color: red;')),
                      footer=tagList(
                          modalButton('Cancel'),
                          actionButton(ns('meta_reset_ok'), 'Yes'))
                      )
                  )
            })

            observeEvent(input$meta_reset_ok, {
              coldata_all$curr <- coldata_all$init
              coldata_all$staging <- coldata_all$init
              removeModal()
            })

            coldata <- eventReactive(c(coldata_all$init,
                                    coldata_all$curr,
                                    coldata_all$staging), {
                coldata_all
            })

            helpButtonServer('edit_metadata_help', size='l')

            return(reactive(coldata()))

        } # function
    ) # moduleServer
}
