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
get_bookmarks_path <- function(lab_project, analysis, bookmarks_root) {
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

find_bookmark_dir <- function(id, bookmarks_root) {
    if (
      length(id) != 1 ||
      is.na(id) ||
      !grepl("^[A-Za-z0-9]+$", id)
    ) {
      stop('Saved view is unavailable.')
    }

    input_files <- list.files(
      bookmarks_root,
      pattern = '^input\\.rds$',
      recursive = TRUE,
      full.names = TRUE
    )

    bookmark_dirs <- dirname(input_files)
    matches <- bookmark_dirs[basename(bookmark_dirs) == id]

    if (length(matches) != 1) {
      stop('Saved view is unavailable.')
    }

    matches[[1]]
}

read_bookmark_metadata <- function(dir) {
    values_file <- file.path(dir, 'values.rds')
    if (!file.exists(values_file)) {
      return(NULL)
    }

    values <- tryCatch(readRDS(values_file), error = function(e) NULL)
    if (is.null(values)) {
      return(NULL)
    }

    list(
      dds = values$dds %||% NA_character_,
      assay_analysis = values$assay_analysis %||% NA_character_,
      label = values$label %||% NA_character_,
      mode = values$mode %||% NA_character_,
      comp = values$comp %||% '',
      saved_at = values$saved_at %||% ''
    )
}

get_saved_bookmark_choices <- function(bookmarks_dir, project, analysis) {
    all_ids <- list.dirs(
      bookmarks_dir,
      recursive = FALSE,
      full.names = FALSE
    )

    matches <- Filter(
      Negate(is.null),
      lapply(all_ids, function(id) {
        meta <- read_bookmark_metadata(file.path(bookmarks_dir, id))
        if (
          is.null(meta) ||
          !identical(meta$dds, project) ||
          !identical(meta$assay_analysis, analysis)
        ) {
          return(NULL)
        }

        list(id = id, meta = meta)
      })
    )

    labels <- vapply(
      matches,
      function(match) {
        name <- if (!is.na(match$meta$label)) {
          match$meta$label
        } else {
          match$meta$saved_at
        }

        paste0(
          name,
          if (nzchar(match$meta$mode)) {
            paste0(' [', match$meta$mode, ']')
          } else {
            ''
          }
        )
      },
      character(1)
    )

    ids <- vapply(matches, function(match) match$id, character(1))

    list(
      has_bookmarks = length(all_ids) > 0,
      choices = setNames(ids, labels)
    )
}
