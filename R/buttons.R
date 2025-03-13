#' Help button module ui
#'
#' UI for help button module
#'
#' @param id Input id
#'
helpButtonUI <- function(id){
  ns <- NS(id)

  config <- get_config()

  actionButton(ns('help_btn'),
               label=NULL,
               icon=icon('question'),
               style=config$style$help_buttons)
}

#' Help button module server
#'
#' Server for download button module
#'
#' @param id Input id
#' @param ... other params passed to helpModal()
#'
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

#' Download button module ui
#'
#' UI for download button module
#'
#' @param id Input id
#'
downloadButtonUI <- function(id){
  ns <- NS(id)

  config <- get_config()

  actionButton(ns('dload_btn'),
               label='Download',
               icon=icon('download'),
               style=config$style$dload_buttons)
}

#' Download button module server
#'
#' Server for download button module
#'
#' @param id Input id
#' @param outplot reactive plot handle
#' @param plot_type reactive/static value used for output filename
#'
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
            save_image(outplot(),
                       file=file,
                       width=input$plot_wd*ppi,
                       height=input$plot_ht*ppi)
          } else {
            ggsave(file, plot = outplot(),
                   device='pdf',
                   width=input$plot_wd, height=input$plot_ht)
          }
          removeModal()
        }
      )

    }   # function
  ) # moduleServer
}

