#' Help button module
#'
#' @description
#' Module UI & server for help buttons.
#'
#' @param id Module id. This also doubles as prefixes for help text files.
#' @param ... other params passed to helpModal()
#'
#' @returns
#' UI returns tagList with help button UI.
#' Server invisibly returns NULL (used for side effects).
#'
#' @examples
#' library(shiny)
#'
#' # app with a single help button to show DE summary table details
#' if(interactive()){
#'   shinyApp(
#'     ui = fluidPage(
#'            helpButtonUI('de_summary_help')
#'          ),
#'     server = function(input, output, session){
#'                helpButtonServer('de_summary_help')
#'              }
#'   )
#' }
#'
#' @rdname helpmod
#' @name helpmod
NULL

#' @rdname helpmod
#' @export
helpButtonUI <- function(id){
  ns <- NS(id)

  config <- get_config()

  actionButton(ns('help_btn'),
               label=NULL,
               icon=icon('question'),
               style=config$style$help_buttons)
}

#' @rdname helpmod
#' @export
helpButtonServer <- function(id, ...){
  moduleServer(
    id,

    function(input, output, session){

      ns <- NS(id)

      observeEvent(input$help_btn, {
        # convert id to full system path
        mdfile <- system.file('extdata', 'help',
                              paste0(sub('_help', '', id), '.md'),
                              package=packageName())

        showModal(
          helpModal(mdfile=mdfile, ...)
        )
      })

    }   # function
  ) # moduleServer
}

#' Help modal
#'
#' This generates a modal dialog that includes text
#' from a markdown file.
#'
#' @param mdfile path to markdown file
#' @param title Title of modal dialog
#' @param ... other params passed to modalDialog()
#'
#' @return Modal dialog with help documentation.
#'
helpModal <- function(mdfile, title=NULL, ...){
  modalDialog(
      title=title,
      includeMarkdown(mdfile),
      footer=tagList(
          modalButton('OK')
      ),
      easyClose=TRUE,
      ...
  )
}

#' Download button module
#'
#' @description
#' Module UI & server for download buttons.
#'
#' @param id Module id
#' @param outplot reactive plot handle
#' @param plot_type reactive/static value used for output filename
#'
#' @returns
#' UI returns tagList with download button UI.
#' Server invisibly returns NULL (used for side effects).
#'
#' @examples
#' library(shiny)
#' library(ggplot2)
#'
#' # get example object
#' obj <- make_example_carnation_object()
#' res <- as.data.frame(obj$res[[1]])
#'
#' # make MA plot
#' p <- ggplot(res, aes(x=baseMean, y=log2foldChange)) +
#'        geom_point(color='black', alpha=0.5)
#'
#' outplot <- reactive({ p })
#'
#' # app with a single button to download a plot
#' if(interactive()){
#'   shinyApp(
#'     ui = fluidPage(
#'            downloadButtonUI('p')
#'          ),
#'     server = function(input, output, session){
#'                downloadButtonServer('p', outplot, 'maplot')
#'              }
#'   )
#' }
#'
#' @rdname dlmod
#' @name dlmod
NULL

#' @rdname dlmod
#' @export
downloadButtonUI <- function(id){
  ns <- NS(id)

  config <- get_config()

  actionButton(ns('dload_btn'),
               label='Download',
               icon=icon('download'),
               style=config$style$dload_buttons)
}

#' @rdname dlmod
#' @export
downloadButtonServer <- function(id, outplot, plot_type){
  moduleServer(
    id,

    function(input, output, session){

      ns <- session$ns

      config <- get_config()

      observeEvent(input$dload_btn, {
        dims <- config$server$de_analysis$pdf

        showModal(
          modalDialog(
            title='Save plot to PDF?',
            tagList(
              fluidRow(
                column(6, 'height (in inches)'),
                column(6,
                  numericInput(ns('plot_ht'),
                               label=NULL,
                               value=dims$height)
                ) # column
              ), # fluidRow
              fluidRow(
                column(6, 'width (in inches)'),
                column(6,
                  numericInput(ns('plot_wd'),
                               label=NULL,
                               value=dims$width)
                ) # column
              ) # fluidRow
            ), # tagList
            footer=tagList(
                downloadButton(ns('download'), label='OK'),
                modalButton('Cancel')
            ),
            easyClose=TRUE
          )
        ) # showModal
      })

      # reactive to return plot filename & handles
      # both reactive and static values
      plot_filename <- reactive({
          if(is.reactive(plot_type)) plot_type()
          else plot_type
      })

      output$download <- downloadHandler(
        filename = function(){
          paste0(plot_filename(), '.pdf')
        },
        content = function(file){
          # dendrogram & upset plot use base R graphics
          if(plot_filename() == 'dendrogram'){
            pdf(file, width=input$plot_wd, height=input$plot_ht)

            old.par <- par(mar = c(0, 0, 1, 25), cex=1.25)
            plot(outplot(), horiz=TRUE)
            par(old.par)

            dev.off()
          } else if(plot_filename() == 'upsetplot'){
            pdf(file, width=input$plot_wd, height=input$plot_ht)
            print(outplot())
            dev.off()
          } else if(inherits(outplot(), 'plotly')){
            # NOTE: this needs 'kaleido' module in python to be
            #       available for 'reticulate'
            ppi <- config$server$de_analysis$heatmap$pdf_res

            # NOTE: turn off mathjax in kaleido to prevent "Loading MathJax ..." box in saved plot
            # - solution from https://stackoverflow.com/questions/79464233/loading-mathjax-extensions-mathmenu-js-box-visible-in-pdf-when-running-plotl/
            k <- plotly::kaleido()
            k$scope$mathjax <- FALSE

            # add chromium flag if running in docker to prevent kaleido gpu crash error
            # see: https://github.com/plotly/Kaleido/issues/74
            if(Sys.getenv('IN_DOCKER') == 'true'){
              k$scope$chromium_args <- c(k$scope$chromium_args, '--single-process')
            }

            k$transform(outplot(),
                        file=file,
                        width=input$plot_wd*ppi,
                        height=input$plot_ht*ppi)
          } else if(inherits(outplot(), 'ggplot')){
            ggsave(file, plot = outplot(),
                   device='pdf',
                   width=input$plot_wd, height=input$plot_ht)
          } else if(inherits(outplot(), 'igraph')){
            pdf(file, width=input$plot_wd, height=input$plot_ht)

            plot(outplot())

            dev.off()
          } else {
            showNotification(
              paste0('Plot type: "', class(outplot()), '" not supported'),
              type='error'
            )
          }
          removeModal()
        }
      )

    }   # function
  ) # moduleServer
}

