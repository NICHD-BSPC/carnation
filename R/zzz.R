# global variable to store carnation python env
carnation <- "r-carnation"

.onAttach <- function(libname, pkgname) {
  shiny::addResourcePath('help',
                         system.file('extdata', 'help',
                                     package=packageName()))

  shiny::addResourcePath("sbs",
                         system.file("www", package = "shinyBS"))
}

.onLoad <- function(...){
  reticulate::use_virtualenv(carnation, required=FALSE)
}


#' Create carnation python environment
#'
#' This function installs 'plotly' and 'kaleido' python packages
#' in an environment to allow PDF downloads from plotly plots.
#'
#' @param envname name of the python environment
#' @param ... parameters passed to reticulate::py_install
#'
#' @return \code{NULL}, invisibly. The function is called for its side effects.
#'
#' @examples
#' if(interactive()){
#'     install_carnation()
#' }
#'
#' @export
install_carnation <- function(envname, ...) {
  if(missing(envname)) envname <- carnation
  reticulate::py_install(c("plotly", "kaleido"),
                         envname = envname, ...)
}
