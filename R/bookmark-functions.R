
##### Storage #####

#' Gets the root folder for bookmark storage
#'
#' Uses the CARNATION_BOOKMARKS_ROOT environment variable if it points to an
#' existing folder, otherwise falls back to ~/.carnation-bookmarks/.
#'
#' @return The bookmarks root path.
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

#' Builds the bookmark folder path for an assay
#'
#' @param lab_project A string containing lab/project, such as current$proj
#' @param analysis The analysis name, such as current$analysis
#' @param bookmarks_root The bookmarks root folder
#'
#' @return bookmarks_root/lab/project/analysis, or NULL if lab_project isn't
#' in lab/project form.
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
#' @return Gives the path to the requested bookmark directory, or `NULL`
#' when the id is invalid or no single matching bookmark is found.
#'
find_bookmark_dir <- function(id, root) {
  if (length(id) != 1 || is.na(id) || !grepl('^[A-Za-z0-9]+$', id)) {
    return(NULL)
  }

  bookmark_dirs <- list_bookmark_dirs(root, recursive = TRUE)
  matches <- bookmark_dirs[basename(bookmark_dirs) == id]

  if (length(matches) == 1) matches[[1]] else NULL
}

#' Reads the stored values for a bookmark
#'
#' @param dir The bookmark folder named for the unique bookmark ID.
#'
#' @return The list saved in values.rds or NULL if it can't be read.
#'
read_bookmark_values <- function(dir) {
  if (is.null(dir)) {
    return(NULL)
  }
  values_file <- file.path(dir, 'values.rds')
  if (!file.exists(values_file)) {
    return(NULL)
  }
  tryCatch(readRDS(values_file), error = function(e) NULL)
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
  values <- read_bookmark_values(dir)
  if (is.null(values)) {
    return(NULL)
  }

  data.frame(
    id = basename(dir),
    project = values$project %||% NA_character_,
    analysis = values$assay_analysis %||% NA_character_,
    label = values$label %||% NA_character_,
    mode = values$mode %||% NA_character_,
    comp = values$comp %||% '',
    modified_at = values$modified_at %||% '',
    created_by = values$user %||% '',
    stringsAsFactors = FALSE
  )
}

#' Lists the saved bookmarks for one analysis, gathering metadata
#' intended for display. Only looks at immediate folders and checks
#' that the project and analysis match the bookmark data. If a
#' bookmark without a label is found it displays the time it
#' was made.
#'
#' @param bookmarks_dir The assay's bookmark folder
#' @param project The lab or project name such as current$proj
#' @param analysis The analysis name such as current$analysis
#'
#' @return A data frame with the character columns id, label, modified_at
#' comp, and mode. One row per matching bookmark.
#'
list_saved_bookmarks <- function(bookmarks_dir, project, analysis) {
  empty_bookmarks <- data.frame(
    id = character(),
    label = character(),
    modified_at = character(),
    comp = character(),
    mode = character(),
    created_by = character()
  )

  bookmark_dirs <- list_bookmark_dirs(bookmarks_dir)
  records <- Filter(
    Negate(is.null),
    lapply(bookmark_dirs, function(dir) {
      tryCatch(
        read_bookmark_metadata(dir),
        error = function(e) {
          warning('Skipping unreadable bookmark: ', basename(dir), call. = FALSE )
          NULL
        })
    })
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
  bookmarks_df$label[missing_label] <- bookmarks_df$modified_at[missing_label]

  bookmarks_df <- bookmarks_df[,
    c('id', 'label', 'modified_at', 'created_by', 'comp', 'mode'),
    drop = FALSE
  ]
  bookmarks_df <- bookmarks_df[
    order(bookmarks_df$modified_at, decreasing = TRUE),
    ,
    drop = FALSE
  ]
  rownames(bookmarks_df) <- NULL

  bookmarks_df
}

#' Deletes a saved bookmark
#'
#' @param id The unique ID of the bookmark to delete.
#' @param bookmarks_dir The bookmark directory to search.
#'
#' @return `TRUE` when the bookmark was deleted, otherwise `FALSE`.
#'
delete_saved_bookmark <- function(id, bookmarks_dir) {
  dir <- find_bookmark_dir(id, bookmarks_dir)

  if (is.null(dir)) {
    return(FALSE)
  }

  unlink(dir, recursive = TRUE) == 0
}

#' Can user modify a bookmark?
#'
#' @param bookmark unique ID of the bookmark to delete.
#' @param username username of who wants to modify a bookmark
#' @param is_admin boolean is the user admin?
#'
#' @return `TRUE` when the bookmark was deleted, otherwise `FALSE`.
#'
can_modify_bookmark <- function(bookmark, username, is_admin) {
  if (isTRUE(is_admin)) return(TRUE)
  if (is.null(bookmark) || nrow(bookmark) != 1) return(FALSE)
  is.character(username) && nzchar(username) && isTRUE(bookmark$created_by == username)
}


##### Labels #####

#' Validates a label for a new bookmark
#'
#' @param label The proposed bookmark label.
#' @param project The project containing the assay.
#' @param analysis The assay analysis name.
#' @param bookmarks_root The root bookmark storage directory.
#' @param username username of whoe wants to modify a bookmark.
#' @param is_admin boolean is the user an admin?
#' @param overwrite boolean overwrite existing bookmark?
#'
#' @return A list containing the validation result, cleaned label, and
#' an error message when validation fails.
#'
validate_new_bookmark_label <- function(
  label,
  project,
  analysis,
  bookmarks_root,
  username,
  is_admin,
  overwrite=FALSE
) {
  # Trim and collapse whitespace; case is kept for display
  label <- gsub('[[:space:]]+', ' ', trimws(label %||% ''))

  if (!nzchar(label)) {
    return(list(
      valid = FALSE,
      label = label,
      message = 'Enter a bookmark label.'
    ))
  }

  bookmarks_dir <- get_bookmarks_path(
    project,
    analysis,
    bookmarks_root
  )

  if (is.null(bookmarks_dir)) {
    return(list(
      valid = FALSE,
      label = label,
      message = 'Unable to determine the bookmark location.'
    ))
  }

  existing_bookmarks <- list_saved_bookmarks(
    bookmarks_dir,
    project,
    analysis
  )

  # Duplicate check ignores case
  dup_row <- tolower(existing_bookmarks$label) == tolower(label)
  if (any(dup_row)) {
    existing <- existing_bookmarks[which(dup_row)[1], , drop = FALSE]
    if(!overwrite){
      return(list(
        valid = FALSE,
        label = label,
        message = paste(
          'Bookmark with that label already exists',
          'for this analysis! To overwrite set "Force overwrite?" to TRUE'
        )
      ))
    }
    if (!can_modify_bookmark(existing, username, is_admin)) {
      return(list(
        valid = FALSE, label = label,
        message = paste0('Bookmark: ',
          label,
          ' was made by: ',
          existing$created_by,
          '. Only they or an admin can overwrite it.')
        )
      )
    }
    return(list(
      valid = TRUE,
      label = label,
      message = paste0(
        'Overwriting bookmark "', label, '"'
      ),
      overwrite_id = existing_bookmarks$id[dup_row][1]
    ))

  }

  list(
    valid = TRUE,
    label = label,
    message = 'Bookmark label is valid...'
  )
}


##### Access Check #####

#' Checks whether the user can access a bookmark
#'
#' This is the access gate for bookmark restoration, nothing from
#' values.rds should be sent to the browser unless this returns status 'ok'.
#'
#' @param values The list stored in a bookmark's values.rds
#' @param assay_list The named list of assays that the main auth of
#' the program has already determined the user can access
#'
#' @return A list with 'status' ('ok', 'unavailable', or 'denied'). When the
#' status is 'ok' it also has 'label' and 'bookmark' (project, analysis,
#' allowed_path and inputs).
#'
check_bookmark_access <- function(values, assay_list) {
  if (length(values) == 0) {
    return(list(status = 'unavailable'))
  }

  # The bookmark's assay must be in the user's assay list. The path to load
  # comes from that list, never from values.rds.
  project_assays <- assay_list[[values$project]]
  matches <- which(names(project_assays) == values$assay_analysis)
  if (length(matches) != 1) {
    return(list(status = 'denied'))
  }

  list(
    status = 'ok',
    label = values$label,
    bookmark = list(
      project = values$project,
      analysis = values$assay_analysis,
      allowed_path = unname(project_assays[matches]),
      inputs = values$inputs
    )
  )
}


##### Restore #####

# Checks whether the requested project/analysis is already the one
# currently loaded, so bookmark restoration can skip reloading it.
is_assay_already_loaded <- function(current, project, analysis) {
  isTRUE(identical(isolate(current$proj), project)) &&
    isTRUE(identical(isolate(current$analysis), analysis))
}

# Sends each saved input value to the browser.
apply_restored_inputs <- function(context, session) {
  saved_inputs <- context$inputs
  if (is.null(saved_inputs) || length(saved_inputs) == 0) {
    return()
  }

  for (id in names(saved_inputs)) {
    if (is.null(saved_inputs[[id]]) ||
        id == 'gene.to.plot') {
      next
    }
    session$sendInputMessage(
      id,
      list(value = saved_inputs[[id]])
    )
  }
}

# Resets all inputs in the session
reset_inputs <- function(session) {
    session$sendCustomMessage('prune_resettables', list())
    shinyjs::reset()
}

#' @rdname bookmarkingmod
#' @export
bookmarkRestoreScript <- function(id) {
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
          Shiny.addCustomMessageHandler(
            'prune_resettables', function(message) {
              $('.shinyjs-resettable').filter(function() {
                var value = $(this).attr('data-shinyjs-resettable-value');
                return $(this).attr('data-shinyjs-resettable-type') === 'Select' &&
                  (value === '\"\"' || value === '[]');
              }).removeClass('shinyjs-resettable');
            }
          )
        "
      )
    )
  )
}

# After the next server flush, asks the browser to confirm that it has reached
# a render frame.
send_bookmark_restore_probe <- function(session, root_session, round) {
  root_session$onFlushed(
    function() {
      if (!root_session$isClosed()) {
        session$sendCustomMessage(
          'bookmark_restore_probe',
          list(round = round)
        )
      }
    },
    once = TRUE
  )
}
