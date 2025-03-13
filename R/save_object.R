#' Save object module UI
#'
#' @param id ID string used to match the ID used to call the module server function
#'
#' @export
saveUI <- function(id){
  ns <- NS(id)

  actionButton(ns('save_rds'), label='Save object', icon=icon('save'),
               class='btn-primary')

}

#' Save object module server function
#'
#' @param id ID string used to match the ID used to call the module UI function
#' @param original original carnation object
#' @param current current carnation object
#' @param coldata reactiveValues object containing object metadata
#' @param pattern regex pattern for finding carnation data
#'
#' @export
saveServer <- function(id, original, current, coldata, pattern){
  moduleServer(
    id,

    function(input, output, session){
        ns <- NS(id)

        reload_parent <- reactiveValues(flag=FALSE)
        save_flag <- reactiveValues(l=FALSE)

        observeEvent(input$save_rds, {
          if(is.null(original$obj)){
            showNotification(
              'Error saving object: Object not loaded yet',
              type='error'
            )

            validate(
              need(!is.null(original$obj) == '', '')
            )
          } else {
            showModal(saveModal())
          }
        })

        saveModal <- function() {
          current_assay_label <- sub(paste0(pattern, '\\.rds$'), '',
                                     basename(original$path),
                                     ignore.case=TRUE)

          modalDialog(
            title = 'Save RDS object',
            textInput(ns('rds_label'), label='Analysis label',
                      value=current_assay_label,
                      width='100%'),
            span('This is the current analysis label by default.',
                 style='font-style: italic;'), br(), br(),

            textInput(ns('rds_path'), label='Destination directory',
                      value=dirname(original$path),
                      width='100%'),
            span('Folder where this file will be saved. This defaults to the current RDS location.',
                 style='font-style: italic;'), br(), br(),

            selectInput(ns('compress'), label='Compress RDS?',
                        choices=c(FALSE, TRUE)),

            span('Uncompressed files load faster, but are larger on disk.',
                 style='font-style: italic;'), br(), br(),

            selectInput(ns('overwrite'), label='Force overwrite?',
                        choices=c(FALSE, TRUE)),
            span('Caution: Setting this to TRUE will overwrite RDS file on disk',
                 style='font-style: italic; color: red;'), br(), br(),

            footer = tagList(
              modalButton("Cancel"),
              actionButton(ns("save_ok"), "OK")
            ),
            easyClose=TRUE
          )
        }

        # shows a 'saving' message that cannot be dismissed
        savingModal <- function() {
          modalDialog(
            span('Saving object, please wait'),

            footer = NULL
          )
        }

        # observer for save RDS button
        observeEvent(input$save_ok, {
          validate(
            need(!is.null(original$path), '')
          )

          if (file.exists(original$path)) {
            if(!dir.exists(input$rds_path)){
              showNotification(
                paste0('Destination directory "', input$rds_path, '" does not exist!'),
                type='error'
              )

              validate(
                need(dir.exists(input$rds_path), '')
              )
            }

            # get matched string from filename
            if(pattern == ''){
              match <- ''
            } else {
              rg <- regexpr(pattern, original$path)
              match <- substr(original$path, rg, rg + attr(rg, 'match.length'))
            }
            destpath <- file.path(input$rds_path,
                                  paste0(input$rds_label, match, '.rds'))

            # check if destination is a symlink
            # - if yes, don't overwrite
            if(file.exists(destpath)){
              link_to <- Sys.readlink(destpath)
              if(link_to == ''){
                # if force overwrite is TRUE, save object
                if(as.logical(input$overwrite)){
                  save_flag$l <- TRUE
                } else {

                  showNotification(
                    'File already exists! Please retry with different analysis label or set "Force overwrite?" to TRUE',
                    type='warning'
                  )

                  validate(
                    need(as.logical(input$overwrite), '')
                  )
                }
              } else if(link_to != destpath){
                showNotification(
                  paste('Destination path is a symlink and cannot be overwritten.',
                        'Please choose different assay label'),
                  type='error'
                )
                validate(
                  need(link_to == '', '')
                )
              } else {
                showNotification(
                  'Error saving object: unknown. Please check destination & source paths',
                  type='error'
                )

                validate(
                  need(link_to == '', '')
                )
              }
            } else {
              save_flag$l <- TRUE
            }
          } else {
            # this can happen if connection to server is lost after data is loaded
            showNotification(
              paste0('Error saving object: File "',
                     original$path, '"does not exist on disk!'),
              type='error'
            )

            validate(
              need(!file.exists(original$path), '')
            )
          }
        })

        # observer to save RDS file
        observeEvent(save_flag$l, {
          validate(
            need(!is.null(original$path), '')
          )

          # Check that data object & file path exists
          if (file.exists(original$path)) {
            # get matched string from filename
            if(pattern == ''){
              match <- ''
            } else {
              rg <- regexpr(pattern, original$path)
              match <- substr(original$path, rg, rg + attr(rg, 'match.length'))
            }
            destpath <- file.path(input$rds_path,
                                  paste0(input$rds_label, match, '.rds'))

            # check if user has write access
            write_access <- file.access(dirname(destpath), mode=2)

            if(write_access != '0'){
              showNotification(
                paste0('Error saving object: Please check to see if you read-write access to save location'),
                type='error'
              )

              validate(
                need(write_access == '0', '')
              )
            }

            showModal(savingModal())

            # recreate object from current app obj before saving
            obj <- reactiveValuesToList(current)

            # attach edited metadata
            # if dds has one element, metadata slot is called 'all_samples'
            # and it won't match dds slot name
            if(length(current$dds) == 1){
              colData(obj[[ 'dds' ]][[1]]) <- coldata[[1]]
              colData(obj[[ 'rld' ]][[1]]) <- coldata[[1]]

              if(!is.null(obj[[ 'all_dds' ]])) colData(obj[[ 'all_dds' ]]) <- coldata[[1]]
              if(!is.null(obj[[ 'all_rld' ]])) colData(obj[[ 'all_rld' ]]) <- coldata[[1]]
            } else {
              # otherwise, the metadata should have slots corresponding
              # to each dds slot in addition to 'all_samples'
              for(name in names(obj[[ 'dds' ]])){
                colData(obj[[ 'dds' ]][[name]]) <- coldata[[name]]
                colData(obj[[ 'rld' ]][[name]]) <- coldata[[name]]
              }
            }

            # NOTE: remove .Environment attributes from @design slots of obj$dds
            #   elements & obj$all_dds
            # - this prevents saved object from becoming very large if another
            #   object has been previously loaded
            obj$dds <- lapply(obj$dds, function(x){
                         attr(x@design, '.Environment') <- NULL
                         x
                       })

            attr(obj$all_dds@design, '.Environment') <- NULL

            saveRDS(obj, destpath, compress=as.logical(input$compress))

            removeModal()
            reload_parent$flag <- TRUE
          }
        })

        trigger <- eventReactive(reload_parent$flag, {
          if(reload_parent$flag){
              list(reload_parent=TRUE)
          }
        })

        return(reactive(trigger()))

    } # function
  ) # moduleServer
}
