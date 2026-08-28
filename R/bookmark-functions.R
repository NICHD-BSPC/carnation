get_bookmarks_root_path <- function() {
    if (Sys.getenv('CARNATION_BOOKMARKS_ROOT') != '') {
      path <- Sys.getenv('CARNATION_BOOKMARKS_ROOT')
      if (!dir.exists(path)) {
        stop(
          'Environment variable "CARNATION_BOOKMARKs_ROOT" exists',
          'but specified location does not exist on disk:', path
        )
      }
    } else {
      path <- path.expand('~/.carnation-bookmarks/')
    }
    path
}

#' Builds path for setting the bookmark save directory
#' to.
#'
#' @param project A string containing lab/project, such as current$proj
#' @param analysis The current analysis name, such as current$analysis
#' @param bookmark_root The root of the bookmark directory you want to use
#'
#' @return Generates a path of bookmarks_root/lab/project/analysis
#'
get_current_bookmarks_path <- function(lab_project, analysis, bookmarks_root) {
    lab_project_split <- strsplit(
    lab_project,
    .Platform$file.sep,
    fixed = TRUE
    )[[1]]

    if (length(lab_project_split) != 2)  {
      stop (
        'Unexpected project path length'
      )
    }

    lab <- lab_project_split[1]
    project <- lab_project_split[2]

    file.path(bookmarks_root, lab, project, analysis)
}
