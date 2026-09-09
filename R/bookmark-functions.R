# TODO add documentation
register_bookmark_handlers <- function(
    input,
    output,
    session,
    assay.list,
    gene.id,
    load_assay,
    current,
    original,
    bookmarks_root,
    settings,
    restore_roundtrips = 7L
  ) {

    # NULL during normal operation. During bookmark restoration, this holds
    # the current stage and all data needed across Shiny's restore functions.
    #
    # The stages used in stage are 'shiny_native_restore', 'load_data',
    # 'replaying_initial_inputs', 'stabilizing_restore', then back to NULL
    restore_context <- reactiveVal(NULL)

    # Is returned into the main app and fed into plotting modules that use a
    # refresh button to gate updating their plots or have potential race conditions
    # when all of their inputs are updated at once. This is used as a reactive
    # dependency and updated/invalidated after a server flushing loop.
    bookmark_restore_tick <- reactiveVal(0L)

    refresh_bookmarks_tick <- reactiveVal(0L)

    bookmark_label <- reactiveVal(NULL)

    setBookmarkExclude(
      c(
        'assay_do',
        'bookmark',
        'bookmark_dropdown',
        'delete_bookmark',
        'delete_bookmark_do',
        'gene.to.plot',
        'confirm_bookmark',
        'bookmark_name',
        'bookmark_restore_ack'
      ),
      session = session
    )

    # Spawns a modal to take the bookmark name and confirmation.
    observeEvent(
      input$bookmark,
      {
        if (!is.null(isolate(restore_context()))) {
          return()
        }

        showModal(
          modalDialog(
            title = 'Make a bookmark of the current view',
            textInput(
              'bookmark_name',
              'Bookmark name',
              placeholder = 'Enter a descriptive name'
            ),
            footer = tagList(
              modalButton('Cancel'),
              actionButton(
                'confirm_bookmark',
                'Save Bookmark',
                class = 'btn-primary'
              )
            ),
            easyClose = TRUE
          )
        )
      },
      ignoreInit = TRUE
    )


    observeEvent(
      input$confirm_bookmark,
      {
        label <- trimws(input$bookmark_name %||% '')

        bookmark_label(
          if (nzchar(label)) label else NA_character_
        )

        removeModal()

        session$doBookmark()
      },
      ignoreInit = TRUE
    )

    # Updates the current bookmark save directory whenever a new dataset is loaded.
    observeEvent(list(current$proj, current$analysis), {
      req(current$proj, current$analysis)
      current_bookmarks_dir <- get_bookmarks_path(
        current$proj,
        current$analysis,
        bookmarks_root
      )
      shinyOptions(
        save.interface = function(id, callback) {
          dir <- file.path(current_bookmarks_dir, id)
          dir.create(dir, recursive = TRUE, showWarnings = FALSE)
          callback(dir)
        }
      )
    })

    # Runs just before bookmarking.
    onBookmark(
      function(state) {
        state$values$project <- current$proj
        state$values$assay <- original$path
        state$values$assay_analysis <- current$analysis
        state$values$label <- bookmark_label()
        state$values$mode <- input$mode
        state$values$comp <- paste(input$comp_all, collapse = ', ')
        state$values$saved_at <- as.character(Sys.time())
        state$values$gene.to.plot <- input$gene.to.plot
      },
      session = session
    )

    # Runs after the bookmark has been generated.
    onBookmarked(
      function(url) {
        updateQueryString(url, session = session)
      },
      session = session
    )
    # Runs at start of loading bookmark/state.
    onRestore(
      function(state) {
        restore_context(list(
          stage = 'shiny_native_restore',
          project = state$values$project %||% state$values$dds,
          analysis = state$values$assay_analysis,
          inputs = as.list(state$input),
          genes = as.list(state$values$gene.to.plot)
        ))
      },
      session = session
    )

    # At this point Shiny has loaded the bookmarked input state from the RDS
    # object. When the data is loaded, side effects will overwrite these
    # inputs, so advance the saved snapshot to the data-loading stage.
    onRestored(
      function(state) {
        context <- isolate(restore_context())
        req(context)

        context$stage <- 'load_data'
        restore_context(context)
      },
      session = session
    )

    # Data loading phase.
    observe({
      context <- restore_context()
      # Uses existing authentication to determine what data should be available.
      allowed_list <- assay.list$l

      req(context)
      req(identical(context$stage, 'load_data'))
      req(allowed_list)

      project_assays <- allowed_list[[context$project]]

      if (is.null(project_assays)) {
        return()
      }

      matches <- which(names(project_assays) == context$analysis)

      if (length(matches) != 1) {
        return()
      }

      allowed_path <- unname(project_assays[matches])

      load_assay(
        context$project,
        context$analysis,
        allowed_path
      )

      context$stage <- 'replaying_initial_inputs'
      restore_context(context)
    })

    # After the next server flush send a message asking
    # the browser to confirm it has reached a render frame.
    send_restore_probe <- function(token, round) {

      session$onFlushed(
        function() {
          if (!session$isClosed()) {
            session$sendCustomMessage(
              'bookmark_restore_probe',
              list(token = token, round = round)
            )
          }
        },
        once = TRUE
      )
    }


    # Reapplies the stored inputs after the data has been reloaded.
    observeEvent(
      restore_context(),
      {
        context <- restore_context()

        if (is.null(context) || !identical(context$stage, 'replaying_initial_inputs')) {
          return()
        }

        apply_restored_inputs(
          context$inputs,
          session,
          context$genes,
          gene.id$gene
        )

        # Is a count of which resoration we are on in this session.
        # Gives identity of which restoration in this session the custom
        # message belongs to. Currently this is mostly future proofing and defensive
        # because loading a bookmark the way it is now will make a new session
        # each time. Keeping it in for now because it costs little.
        restore_token <- as.character(
          isolate(bookmark_restore_tick()) + 1L
        )

        context$stage <- 'stabilizing_restore'
        context$restore_token <- restore_token
        context$restore_round <- 1L
        restore_context(context)

        send_restore_probe(restore_token, 1L)
      },
      #
      priority = -100
    )

    # Handles and advances the server-flush/browser update one
    # round trip at a time. Ignores stale, duplicate or out of order
    # acknowledgements.
    observeEvent(
      input$bookmark_restore_ack,
      {
        context <- isolate(restore_context())
        ack <- input$bookmark_restore_ack

        if (
          is.null(context) ||
          !identical(context$stage, 'stabilizing_restore') ||
          !is.list(ack) ||
          is.null(ack$token) ||
          !identical(ack$token, context$restore_token) ||
          !isTRUE(ack$round == context$restore_round)
        ) {
          return()
        }

        current_round <- as.integer(ack$round)

        # When the correct number of roundtrips is reached,
        # increment bookmark_restore tick to trigger all the
        # plot refreshes and other things that depend on it.
        # This completes the session restore
        if (current_round >= restore_roundtrips) {
          bookmark_restore_tick(
            isolate(bookmark_restore_tick()) + 1L
          )

          restore_context(NULL)
          return()
        }

        next_round <- current_round + 1L
        context$restore_round <- next_round
        restore_context(context)

        # On the second to last round, applies stored inputs again.
        # This allows to the system to reach settings that are gated
        # by other settings being needed to be selected first.
        # When and the number of times this happens may need to
        # be tuned.
        if (current_round == restore_roundtrips - 2L ||
            current_round == 3L) {
          apply_restored_inputs(
            context$inputs,
            session,
            context$genes,
            gene.id$gene
          )
        }

        send_restore_probe(context$restore_token, next_round)
      },
      ignoreInit = TRUE,
      priority = -100
    )

    ########### Bookmark loading UI ###############

    ui_selected_assay <- reactive({
      req(input$dds, input$assay)

      project_assays <- assay.list$l[[input$dds]]
      idx <- match(input$assay, unname(project_assays))

      req(!is.na(idx))

      list(
        project = input$dds,
        analysis = names(project_assays)[idx],
        path = unname(project_assays[idx])
      )
    })

    ui_selected_assay_bookmarks <- reactive({
      req(ui_selected_assay())
      refresh_bookmarks_tick()

      selected <- ui_selected_assay()

      selected_bookmarks_dir <- get_bookmarks_path(
        selected$project,
        selected$analysis,
        bookmarks_root
      )

      list_saved_bookmarks(
        selected_bookmarks_dir,
        selected$project,
        selected$analysis
      )
    })

    # UI for bookmarks description data table.
    output$bookmarks_desc <- renderDT({
      bookmarks_df <- ui_selected_assay_bookmarks()

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
          columnDefs = list(
            list(
              visible = FALSE,
              targets = hidden_columns
            )
          )
        )
      )
    })

    bookmarks_desc_proxy <- dataTableProxy('bookmarks_desc')


    table_bookmark_selected <- eventReactive(
      input$bookmarks_desc_rows_selected,
      {
        bookmarks_df <- ui_selected_assay_bookmarks()
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
        bookmarks_df <- ui_selected_assay_bookmarks()
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

    # UI for the dropdown bookmarks list.
    output$saved_bookmarks_ui <- renderUI({
      bookmarks_df <- ui_selected_assay_bookmarks()

      if (nrow(bookmarks_df) == 0) {
        return(tags$p('No saved bookmarks for this analysis yet.'))
      }

      choices <- setNames(bookmarks_df$id, bookmarks_df$label)

      tagList(
        selectInput(
          'bookmark_dropdown',
          'Saved Bookmarks',
          choices = choices
        ),
        actionButton('go_bookmark', 'Load Bookmark', class = 'btn-primary'),
        if (settings()$is_admin) {
          actionButton('delete_bookmark', 'Delete Bookmark', class = 'btn-primary')
        }
      )
    })

    observeEvent(input$go_bookmark, {
      req(input$bookmark_dropdown)
      updateQueryString(
        paste0('?_state_id_=', input$bookmark_dropdown),
        mode = 'push',
        session = session
      )
      session$reload()
    })

    observeEvent(input$delete_bookmark, {
      req(input$bookmark_dropdown)
      if (!settings()$is_admin) {
        return()
      }

      showModal(
        modalDialog(
          title = 'Delete bookmark',
          'Are you sure you want to delete this bookmark?',
          footer = tagList(
            modalButton('Cancel'),
            actionButton(
              'delete_bookmark_do',
              'Delete Bookmark',
              class = 'btn-danger'
            )
          ),
          easyClose = TRUE
        )
      )
    })

    observeEvent(input$delete_bookmark_do, {
      req(input$bookmark_dropdown)
      if (!settings()$is_admin) {
        return()
      }

      assay <- ui_selected_assay()
      bookmarks_dir <- get_bookmarks_path(
        assay$project,
        assay$analysis,
        bookmarks_root
      )
      dir <- find_bookmark_dir(input$bookmark_dropdown, bookmarks_dir)
      delete_status <- unlink(dir, recursive = TRUE)

      if (delete_status == 0) {
        removeModal()
        refresh_bookmarks_tick(
          isolate(refresh_bookmarks_tick()) + 1L
        )
        showNotification('Bookmark deleted.')
      } else {
        showNotification('Unable to delete bookmark.', type = 'error')
      }
    })

    # reavtiveVals returned and used by other parts of the
    # dashboard.
    list(
      restore_context = restore_context,
      restore_tick = bookmark_restore_tick
    )

}

################# Helper Functions #######################

#' Checks if environment variable CARNATION_BOOKMARKS_ROOT is set,
#' if it is, it sets the root of the bookmarks storage to that
#' directory. If it is not, it sets the bookmark directory to
#' .carnation-bookmarks. If using the environment variable,
#' that directory must already exist.
#'
get_bookmarks_root_path <- function() {
    if (Sys.getenv('CARNATION_BOOKMARKS_ROOT') != '') {
      path <- Sys.getenv('CARNATION_BOOKMARKS_ROOT')
      if (!dir.exists(path)) {
        stop(
          'Environment variable "CARNATION_BOOKMARKS_ROOT" exists',
          'but specified location does not exist on disk:', path
        )
      }
    } else {
      path <- path.expand('~/.carnation-bookmarks/')
    }
    path
}

#' Builds path to setting the bookmark save directory
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
      stop('Saved view is unavailable.')
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

    bookmarks_df <- bookmarks_df[
      ,
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
apply_restored_inputs <- function(saved_inputs, session, genes, gene_choices) {
  if (length(genes) > 0) {
    updateSelectizeInput(
      session,
      'gene.to.plot',
      choices = gene_choices,
      selected = genes,
      server = TRUE
    )
  }

  if (is.null(saved_inputs) || length(saved_inputs) == 0) {
    return()
  }

  skip_pat <- "(-help_btn|-dload_btn|_do$|^intro$|^bookmark|^go_bookmark|^load_bookmark|^assay_do|^dds$|^assay$|^data_type$|select_all|select_none)"

  replayable <- saved_inputs[
    !grepl(skip_pat, names(saved_inputs))
  ]

  for (id in names(replayable)) {
    if (is.null(replayable[[id]])) {
      next
    }

    session$sendInputMessage(
      id,
      list(value = replayable[[id]])
    )
  }
}
