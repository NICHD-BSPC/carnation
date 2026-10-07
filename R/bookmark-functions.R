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

  conditionalPanel(
    'output.has_bookmarks', ns = ns,
    selectInput(ns('bookmark_dropdown'), 'Saved Bookmarks', choices = NULL),
    fluidRow(style='margin-bottom: 10px',
      column(6, align='left',
        actionButton(ns('go_bookmark'), 'Load Bookmark', class = 'btn-primary')
      ), # column
      column(6, align='right',
        conditionalPanel(
          'output.is_admin', ns = ns,
          actionButton(ns('delete_bookmark'), 'Delete Bookmark', class = 'btn-primary')
        )
      ) # column
    ) # fluidRow
  )
}

bookmarksDescriptionUI <- function(id) {
  ns <- NS(id)

  DTOutput(ns('bookmarks_desc'))
}

bookmark_name_modal <- function(ns) {
  modalDialog(
    title = 'Make a bookmark of the current view',
    textInput(
      ns('bookmark_name'),
      'Bookmark name',
      placeholder = 'Enter a descriptive name'
    ),
    footer = tagList(
      modalButton('Cancel'),
      actionButton(
        ns('confirm_bookmark'),
        'Save Bookmark',
        class = 'btn-primary'
      )
    ),
    easyClose = TRUE
  )
}

bookmark_loading_modal <- function() {
  modalDialog(
    icon('spinner', class = 'fa-spin'),
    'Loading bookmark, please wait',
    footer = NULL,
    easyClose = FALSE
  )
}

bookmark_access_denied_modal <- function() {
  modalDialog(
    title = 'You do not have access to this data area.',
    easyClose = FALSE
  )
}

bookmark_delete_modal <- function(ns) {
  modalDialog(
    title = 'Delete bookmark',
    'Are you sure you want to delete this bookmark?',
    footer = tagList(
      modalButton('Cancel'),
      actionButton(
        ns('delete_bookmark_do'),
        'Delete Bookmark',
        class = 'btn-danger'
      )
    ),
    easyClose = TRUE
  )
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

# Builds the bookmark description table displayed in the saved-bookmarks UI.
bookmark_description_table <- function(bookmarks_df) {
  hidden_columns <- match(
    c('id', 'comp', 'mode'),
    names(bookmarks_df)
  ) - 1

  datatable(
    bookmarks_df,
    rownames = FALSE,
    selection = 'single',
    caption = tags$caption(
      style = 'font-weight: bold; font-size: 15px;',
      'Available bookmarks'
    ),
    options = list(
      dom = 'tlfip',
      columnDefs = list(
        list(
          visible = FALSE,
          targets = hidden_columns
        )
      )
    )
  )
}

# Registers the reactive outputs and observers used to keep the saved-bookmark
# table and dropdown synchronized.
setup_saved_bookmarks_ui <- function(
  input,
  output,
  session,
  bookmarks,
  is_admin,
  loaded_label = function() NULL
) {
  output$bookmarks_desc <- renderDT({
    bookmark_description_table(bookmarks())
  })

  bookmarks_desc_proxy <- dataTableProxy(
    'bookmarks_desc',
    session = session
  )

  table_bookmark_selected <- eventReactive(
    input$bookmarks_desc_rows_selected,
    {
      bookmarks_df <- bookmarks()
      selection <- input$bookmarks_desc_rows_selected

      selected_bookmark <- bookmarks_df[selection, ]
      selected_bookmark$id
    }
  )

  observeEvent(table_bookmark_selected(), {
    selection <- input$bookmarks_desc_rows_selected
    req(length(selection) == 1)
    updateSelectizeInput(
      session,
      inputId = 'bookmark_dropdown',
      selected = table_bookmark_selected()
    )
  })

  # Updates the bookmark table selection based on dropdown selection.
  observeEvent(
    list(input$bookmark_dropdown, input$bookmarks_desc_rows_all),
    {
      bookmarks_df <- bookmarks()
      row_idx <- match(input$bookmark_dropdown, bookmarks_df$id)

      req(!is.na(row_idx))

      displayed_position <- match(
        row_idx,
        input$bookmarks_desc_rows_all
      )

      if (is.na(displayed_position)) {
        bookmarks_desc_proxy %>% selectRows(NULL)
        return()
      }

      page_length <- input$bookmarks_desc_state$length
      if (
        is.null(page_length) ||
          !is.numeric(page_length) ||
          page_length < 1
      ) {
        page_length <- 10
      }
      page_idx <- floor((displayed_position - 1) / page_length) + 1

      bookmarks_desc_proxy %>%
        selectRows(NULL) %>%
        selectRows(row_idx) %>%
        selectPage(page_idx)
    }
  )

  output$has_bookmarks <- reactive(nrow(bookmarks()) > 0)
  output$is_admin <- is_admin
  outputOptions(output, 'has_bookmarks', suspendWhenHidden = FALSE)
  outputOptions(output, 'is_admin', suspendWhenHidden = FALSE)

  observe({
    bookmarks_df <- bookmarks()
    sel <- bookmarks_df$id[match(loaded_label(), bookmarks_df$label)]
    updateSelectInput(
      session,
      'bookmark_dropdown',
      choices = setNames(bookmarks_df$id, bookmarks_df$label),
      selected = if (length(sel) == 1 && !is.na(sel)) sel else NULL
    )
  })

  invisible(NULL)
}

################# Helper Functions #######################

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

#' Resolves an assay selected in the load-data interface
#'
#' @param assay_list The named list of available assays per project.
#' @param project The selected project identifier.
#' @param assay_path The selected assay path.
#'
#' @return A list containing the project, analysis, and path, or `NULL` when
#' the selection cannot be resolved.
#'
resolve_selected_assay <- function(assay_list, project, assay_path) {
  if (
    length(project) != 1L ||
      is.na(project) ||
      length(assay_path) != 1L ||
      is.na(assay_path)
  ) {
    return(NULL)
  }

  project_assays <- assay_list[[project]]
  if (is.null(project_assays)) {
    return(NULL)
  }

  idx <- match(assay_path, unname(project_assays))
  if (is.na(idx)) {
    return(NULL)
  }

  list(
    project = project,
    analysis = names(project_assays)[idx],
    path = unname(project_assays[idx])
  )
}

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

#' Canonicalizes whitespace in a bookmark label
#'
#' @param label A character vector of bookmark labels.
#'
#' @return Labels with surrounding whitespace removed and internal whitespace
#' collapsed to single spaces.
#'
canonicalize_bookmark_label <- function(label) {
  label <- trimws(label)
  gsub('[[:space:]]+', ' ', label)
}

#' Normalizes a bookmark label for comparison
#'
#' @param label A character vector of bookmark labels.
#'
#' @return The canonicalized, lowercase labels.
#'
normalize_bookmark_label <- function(label) {
  tolower(canonicalize_bookmark_label(label))
}

#' Checks whether a bookmark label is already in use
#'
#' @param label The proposed bookmark label.
#' @param bookmarks A bookmark data frame returned by
#'   `list_saved_bookmarks()`.
#'
#' @return `TRUE` when the normalized label already exists.
#'
bookmark_label_exists <- function(label, bookmarks) {
  normalized_label <- normalize_bookmark_label(label)
  existing_labels <- normalize_bookmark_label(bookmarks$label)

  any(existing_labels == normalized_label, na.rm = TRUE)
}

#' Validates a label for a new bookmark
#'
#' @param label The proposed bookmark label.
#' @param project The project containing the assay.
#' @param analysis The assay analysis name.
#' @param bookmarks_root The root bookmark storage directory.
#'
#' @return A list containing the validation result, canonicalized label, and
#' an error message when validation fails.
#'
validate_new_bookmark_label <- function(
  label,
  project,
  analysis,
  bookmarks_root
) {
  label <- canonicalize_bookmark_label(label %||% '')

  if (!nzchar(label)) {
    return(list(
      valid = FALSE,
      label = label,
      message = 'Enter a bookmark name.'
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

  if (bookmark_label_exists(label, existing_bookmarks)) {
    return(list(
      valid = FALSE,
      label = label,
      message = paste(
        'A bookmark with that name already exists',
        'for this assay.'
      )
    ))
  }

  list(
    valid = TRUE,
    label = label,
    message = NULL
  )
}

# Checks whether the requested project/analysis is already the one
# currently loaded, so bookmark restoration can skip reloading it.
is_assay_already_loaded <- function(current, project, analysis) {
  isTRUE(identical(isolate(current$proj), project)) &&
    isTRUE(identical(isolate(current$analysis), analysis))
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
