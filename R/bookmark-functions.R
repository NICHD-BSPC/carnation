bookmarkButtonUI <- function(id) {
  ns <- NS(id)

  actionButton(
    ns('bookmark'),
    label = 'Bookmark',
    class = 'btn-primary'
  )
}

savedBookmarksUI <- function(id) {
  ns <- NS(id)

  uiOutput(ns('saved_bookmarks_ui'))
}

bookmarksDescriptionUI <- function(id) {
  ns <- NS(id)

  DTOutput(ns('bookmarks_desc'))
}

bookmarkRestoreUI <- function(id) {
  ns <- NS(id)

  tags$script(
    HTML(
      paste0(
        "
          Shiny.addCustomMessageHandler(
            'bookmark_restore_probe',
            function(message) {
              window.requestAnimationFrame(function() {
                Shiny.setInputValue(
                  '",
        ns('bookmark_restore_ack'),
        "',
                  message,
                  { priority: 'event' }
                );
              });
            }
          );
        "
      )
    )
  )
}

################# Helper Functions #######################

#' Checks whether the given project and analysis are accessible
#' in the assay list and returns the file path if so.
#'
#' @param assay_list The named list of available assays per project.
#' @param project The project identifier.
#' @param analysis The analysis name.
#'
#' @return The assay path if access is valid, otherwise NULL.
#'
resolve_bookmark_access <- function(assay_list, project, analysis) {
  project_assays <- assay_list[[project]]
  if (is.null(project_assays)) return(NULL)

  matches <- which(names(project_assays) == analysis)
  if (length(matches) != 1) return(NULL)

  unname(project_assays[matches])
}

#' Checks if environment variable CARNATION_BOOKMARKS_ROOT is set,
#' if it is, it sets the root of the bookmarks storage to that
#' directory. If it is not, it sets the bookmark directory to
#' .carnation-bookmarks. If using the environment variable,
#' that directory must already exist.
#'
get_bookmarks_root_path <- function() {
  env_path <- Sys.getenv('CARNATION_BOOKMARKS_ROOT')
  if (env_path != '') {
    if (dir.exists(env_path)) {
      return(env_path)
    }
    message(
      'Environment variable "CARNATION_BOOKMARKS_ROOT" exists ',
      'but specified location does not exist on disk: ',
      env_path, '. Falling back to default.'
    )
  }
  path.expand('~/.carnation-bookmarks/')
}

#' Builds path to setting the bookmark save directory
#' to.
#'
#' @param lab_project A string containing lab/project, such as current$proj
#' @param analysis The current analysis name, such as current$analysis
#' @param bookmarks_root The root of the bookmark directory you want to use
#'
#' @return Generates a path of bookmarks_root/lab/project/analysis
#'
get_bookmarks_path <- function(lab_project, analysis, bookmarks_root) {
  lab_project_split <- strsplit(
    lab_project,
    .Platform$file.sep,
    fixed = TRUE
  )[[1]]

  if (length(lab_project_split) != 2) {
    message('Unexpected project path format: ', lab_project)
    return(NULL)
  }

  lab <- lab_project_split[1]
  project <- lab_project_split[2]

  file.path(bookmarks_root, lab, project, analysis)
}

#' Lists directories containing Shiny bookmark input files.
#'
#' @param root The directory to search from.
#' @param recursive Whether to search below immediate child directories.
#'
#' @return A character vector of bookmark directory paths.
#'
list_bookmark_dirs <- function(root, recursive = FALSE) {
  if (!dir.exists(root)) {
    return(character())
  }

  bookmark_dirs <- list.dirs(
    root,
    recursive = recursive,
    full.names = TRUE
  )

  bookmark_dirs[
    file.exists(file.path(bookmark_dirs, 'input.rds'))
  ]
}

#' Searches recursively from the root to find the
#' directory for the bookmark id provided.
#'
#' @param id The unique id for a bookmark.
#' @param root The root directory to seach from. Can be used from
#' the bookmark root level or the project root level. But with a large
#' number of bookmarks a smaller search area will be better.
#'
#' @return Gives the path to the requested bookmark directory.
#'
find_bookmark_dir <- function(id, root) {
  if (
    length(id) != 1 ||
      is.na(id) ||
      !grepl("^[A-Za-z0-9]+$", id)
  ) {
    message('Saved view is unavailable.')
    return()
  }

  bookmark_dirs <- list_bookmark_dirs(root, recursive = TRUE)
  matches <- bookmark_dirs[basename(bookmark_dirs) == id]

  if (length(matches) != 1) {
    stop('Saved view is unavailable.')
  }

  matches[[1]]
}

#' Extracts metadata from values.rds and returns it as a one-row data frame.
#' Does not contain the data itself, just the names of certain objects
#' and other metadata.
#'
#' @param dir The bookmark folder, named for the unique bookmark ID.
#'
#' @return A one-row data frame of normalized bookmark metadata, or `NULL`
#' when the metadata cannot be read.
#'
read_bookmark_metadata <- function(dir) {
  values_file <- file.path(dir, 'values.rds')
  if (!file.exists(values_file)) {
    return(NULL)
  }

  values <- tryCatch(readRDS(values_file), error = function(e) NULL)
  if (is.null(values)) {
    return(NULL)
  }

  data.frame(
    id = basename(dir),
    project = values$project %||% values$dds %||% NA_character_,
    analysis = values$assay_analysis %||% NA_character_,
    label = values$label %||% NA_character_,
    mode = values$mode %||% NA_character_,
    comp = values$comp %||% '',
    saved_at = values$saved_at %||% '',
    stringsAsFactors = FALSE
  )
}

list_saved_bookmarks <- function(bookmarks_dir, project, analysis) {
  empty_bookmarks <- data.frame(
    id = character(),
    label = character(),
    saved_at = character(),
    comp = character(),
    mode = character()
  )

  bookmark_dirs <- list_bookmark_dirs(bookmarks_dir)
  records <- Filter(
    Negate(is.null),
    lapply(bookmark_dirs, read_bookmark_metadata)
  )

  if (length(records) == 0) {
    return(empty_bookmarks)
  }

  bookmarks_df <- do.call(rbind, records)
  bookmarks_df <- bookmarks_df[
    !is.na(bookmarks_df$project) &
      bookmarks_df$project == project &
      !is.na(bookmarks_df$analysis) &
      bookmarks_df$analysis == analysis,
    ,
    drop = FALSE
  ]

  if (nrow(bookmarks_df) == 0) {
    return(empty_bookmarks)
  }

  missing_label <- is.na(bookmarks_df$label) |
    trimws(bookmarks_df$label) == ''
  bookmarks_df$label[missing_label] <- bookmarks_df$saved_at[missing_label]

  bookmarks_df <- bookmarks_df[,
    c('id', 'label', 'saved_at', 'comp', 'mode'),
    drop = FALSE
  ]
  bookmarks_df <- bookmarks_df[
    order(bookmarks_df$saved_at, decreasing = TRUE),
    ,
    drop = FALSE
  ]
  rownames(bookmarks_df) <- NULL

  bookmarks_df
}


# Takes the stored inputs and reapplies them to each
# type of field.
apply_restored_inputs <- function(saved_inputs, session) {
  if (is.null(saved_inputs) || length(saved_inputs) == 0) {
    return()
  }

  updateSelectizeInput(
    session,
    'gene.to.plot',
    choices = saved_inputs$gene.to.plot,
    selected = saved_inputs$gene.to.plot,
    server = TRUE
  )

  for (id in names(saved_inputs)) {
    if (is.null(saved_inputs[[id]])) {
      next
    }
    session$sendInputMessage(
      id,
      list(value = saved_inputs[[id]])
    )
  }
}
