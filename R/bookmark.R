#' Bookmarking module
#'
#' @description
#' UI & server for module to save, restore and delete bookmarks of the
#' dashboard state. Bookmarks are stored on the server, one set per
#' project/analysis. `bookmarkButtonUI()` saves a bookmark,
#' `bookmarksDescriptionUI()` shows a table of saved bookmarks along
#' with load and delete buttons, and `bookmarkRestoreScript()` is a
#' script that needs to be embedded in the browser for restoring.
#'
#' @param id Module id
#' @param root_input the app's top-level `input` (not the module's own); all
#' inputs are saved with the bookmark
#' @param assay.list reactiveValues object with `l`, the list of assays the
#' user can access
#' @param load_assay function called as
#' `load_assay(project, analysis, path, reset_tabs)` to load the bookmarked assay
#' @param current reactiveValues object with 'proj' (lab/project) & 'analysis'
#' of the loaded assay
#' @param bookmarks_root root folder for bookmark storage, usually from
#' `get_bookmarks_root_path()`
#' @param settings reactive list with user settings; 'is_admin' allows deleting
#' or overwriting any users bookmarks
#' @param gene_choices reactive containing the list of genes that are availble for
#' the gene scratchpad to choose from
#' @param username string containing the session username, allows for user specific control
#' @param input_applies number of times saved inputs are applied during a restore
#' @param rounds_per_apply number of browser round trips after each apply
#'
#' @details
#' The app must call `enableBookmarking("server")`, define `ui` as
#' `function(request)` and set `shinyOptions(load.interface = ...)` using
#' `find_bookmark_dir()`. `bookmarkRestoreScript()` must be included in the UI,
#' otherwise restores never finish.
#'
#' @returns
#' UI functions return tags with bookmarking UI.
#' Server returns a list of reactives: 'is_restoring' (TRUE while a restore
#' is running), 'bookmark_restore_tick' (counter incremented after a restore,
#' passed to plot modules) & 'loaded_bookmark_label' (label of the loaded
#' bookmark or NULL).
#'
#' @examplesIf interactive()
#' library(shiny)
#'
#' bookmarks_root <- tempdir()
#' shinyOptions(
#'   load.interface = function(id, callback) {
#'     dir <- find_bookmark_dir(id, bookmarks_root)
#'     if (!is.null(dir)) callback(dir)
#'   }
#' )
#'
#' ui <- function(request) {
#'   fluidPage(
#'     bookmarkRestoreScript('bookmarking'),
#'     bookmarkButtonUI('bookmarking'),
#'     bookmarksDescriptionUI('bookmarking')
#'   )
#' }
#'
#' server <- function(input, output, session) {
#'   # Simulate app state
#'   current <- reactiveValues(proj = 'lab/project', analysis = 'main')
#'   assay.list <- reactiveValues(l = NULL)
#'   settings <- reactive(list(is_admin = TRUE))
#'
#'   bookmarkingServer(
#'     'bookmarking', input,
#'     assay.list = assay.list,
#'     load_assay = function(...) NULL,
#'     current = current,
#'     bookmarks_root = bookmarks_root,
#'     settings = settings
#'   )
#' }
#'
#' shinyApp(ui, server, enableBookmarking = 'server')
#'
#' @name bookmarkingmod
#' @rdname bookmarkingmod
NULL

#' @rdname bookmarkingmod
#' @export
bookmarkingServer <- function(
  id,
  root_input,
  assay.list,
  load_assay,
  current,
  bookmarks_root,
  settings,
  gene_choices,
  username,
  input_applies = 3L,
  rounds_per_apply = 2L
) {
  moduleServer(
    id,
    function(input, output, session) {
      ##### Shared State #####

      root_session <- session$rootScope()

      # NULL when no restore is running. Otherwise holds the bookmark's project,
      # analysis, allowed_path and inputs, plus the current stage and round.
      # Stages, in order:
      #   'loading_data'     load the bookmarked assay unless it's already loaded
      #   'applying_inputs'  first replay of the saved inputs
      #   'stabilizing'      browser round trips, re-applying inputs every
      #                      rounds_per_apply rounds
      # then back to NULL.
      restore_context <- reactiveVal(NULL)

      # Inputs are applied input_applies times in total, each followed by
      # rounds_per_apply browser round trips.
      total_restore_rounds <- input_applies * rounds_per_apply

      # Label and copy of all inputs for the bookmark being saved; cleared in
      # onBookmarked. Inputs are stored in values.rds instead of Shiny's
      # input.rds so they only reach the browser after the access check.
      pending_save <- reactiveVal(NULL)

      # Incremented when a restore finishes; plot modules use it to refresh.
      bookmark_restore_tick <- reactiveVal(0L)

      # Incremented on save/delete to refresh the saved-bookmarks list.
      saved_list_tick <- reactiveVal(0L)

      # id, label, project and analysis of the bookmark last saved or restored.
      loaded_bookmark <- reactiveVal(NULL)

      ##### Selected Assay #####

      # The assay picked in the load-data dropdowns (not necessarily the one
      # loaded). Only resolves for assays present in the user's assay.list$l.
      assay_selection <- reactive({
        req(root_input$dds, root_input$assay)

        project_assays <- assay.list$l[[root_input$dds]]
        idx <- match(root_input$assay, unname(project_assays))
        req(!is.na(idx))

        list(project = root_input$dds, analysis = names(project_assays)[idx])
      })

      # Bookmark folder for the selected assay
      # (bookmarks_root/lab/project/analysis).
      selected_bookmarks_dir <- reactive({
        selected <- assay_selection()
        dir <- get_bookmarks_path(
          selected$project,
          selected$analysis,
          bookmarks_root
        )
        req(dir)
        dir
      })

      # Saved bookmarks directly inside the selected assay's folder whose
      # values.rds matches that assay. Refreshes on save/delete.
      selected_assay_bookmarks <- reactive({
        selected <- assay_selection()
        saved_list_tick()

        list_saved_bookmarks(
          selected_bookmarks_dir(),
          selected$project,
          selected$analysis
        )
      })

      ##### Saving Bookmarks #####

      # Updates the current bookmark save directory whenever a new dataset is loaded.
      observeEvent(list(current$proj, current$analysis), {
        req(current$proj, current$analysis)
        current_bookmarks_dir <- get_bookmarks_path(
          current$proj,
          current$analysis,
          bookmarks_root
        )
        req(current_bookmarks_dir)
        shinyOptions(
          save.interface = function(id, callback) {
            dir <- file.path(current_bookmarks_dir, id)
            dir.create(dir, recursive = TRUE, showWarnings = FALSE)
            callback(dir)
          }
        )
      })

      # Spawns a modal to take the bookmark name and confirmation.
      observeEvent(
        input$bookmark,
        {
          if (!is.null(restore_context())) {
            return()
          }
          showModal(bookmark_name_modal(
            session,
            loaded_bookmark()$label %||% ''
          ))
        },
        ignoreInit = TRUE
      )

      observeEvent(
        input$confirm_bookmark,
        {
          req(current$proj, current$analysis)

          # Checks for unique label, and if the user is allowed to overwrite.
          label_validation <- validate_new_bookmark_label(
            input$bookmark_name,
            current$proj,
            current$analysis,
            bookmarks_root,
            username(),
            is_admin(),
            overwrite=as.logical(input$overwrite)
          )

          if (!label_validation$valid) {
            showNotification(
              label_validation$message,
              type = 'error'
            )
            return()
          } else {
              showNotification(
                label_validation$message,
              )
          }


          removeModal()

          # Excludes all inputs from Shiny's input.rds; onBookmark stores them
          # in values.rds instead. This module's own inputs (such as table
          # selection) are left out so a restore doesn't overwrite the
          # loaded-bookmark selection with stale values.
          all_inputs <- reactiveValuesToList(root_input, all.names = TRUE)
          own_inputs <- startsWith(names(all_inputs), session$ns(''))
          pending_save(list(
            label = label_validation$label,
            inputs = all_inputs[!own_inputs],
            overwrite_id = label_validation$overwrite_id
          ))
          setBookmarkExclude(
            names(all_inputs),
            session = root_session
          )

          root_session$doBookmark()
        },
        ignoreInit = TRUE
      )

      onBookmark(
        function(state) {
          state$values$project <- current$proj
          state$values$assay_analysis <- current$analysis
          state$values$label <- pending_save()$label
          state$values$mode <- root_input$mode
          state$values$comp <- paste(root_input$comp_all, collapse = ', ')
          state$values$saved_at <- format(Sys.time(), '%y-%m-%d, %H:%M')
          state$values$inputs <- pending_save()$inputs
          state$values$user <- username()
        },
        session = root_session
      )

      # Runs after the bookmark has been generated.
      onBookmarked(
        function(url) {
          saved_label <- pending_save()$label

          updateQueryString(url)
          loaded_bookmark(list(
            id = parseQueryString(sub('^[^?]*', '', url))[['_state_id_']],
            label = saved_label,
            project = isolate(current$proj),
            analysis = isolate(current$analysis)
          ))
          showNotification(
            paste0('Bookmark "', saved_label, '" saved.'),
            type = 'message'
          )
          showModal(show_bookmark_url_modal(url))

          analysis_bookmark_dir <- get_bookmarks_path(isolate(current$proj),
            isolate(current$analysis), bookmarks_root)
          old_id <- pending_save()$overwrite_id
          if (!is.null(old_id) && !delete_saved_bookmark(old_id, analysis_bookmark_dir)) {
            showNotification('Saved, but could not remove the old bookmark.', type = 'warning')
          }
          saved_list_tick(isolate(saved_list_tick()) + 1L)
          pending_save(NULL)
        },
        session = root_session
      )

      ##### Restoring Bookmarks #####

      # Single entry point for every bookmark restore, whether from a URL
      # (onRestored) or from the saved bookmarks UI within the current session.
      start_bookmark_restore <- function(values, id) {
        access <- check_bookmark_access(values, isolate(assay.list$l))
        switch(
          access$status,
          unavailable = showNotification(
            'Bookmark is unavailable.',
            type = 'warning'
          ),
          denied = showModal(bookmark_access_denied_modal()),
          ok = {
            loaded_bookmark(list(
              id = id,
              label = access$label,
              project = access$bookmark$project,
              analysis = access$bookmark$analysis
            ))
            restore_context(c(
              access$bookmark,
              list(stage = 'loading_data', round = 0L)
            ))
          }
        )
      }

      # Runs once after a URL load or page reload. Hands Shiny's restored values
      # to start_bookmark_restore(), which checks access before anything is applied.
      onRestored(
        function(state) {
          # state$dir is the bookmark folder, named for its id. It is NULL
          # when load.interface couldn't find the bookmark.
          start_bookmark_restore(
            state$values,
            if (!is.null(state$dir)) basename(state$dir)
          )
        },
        session = root_session
      )

      # Loads the bookmarked assay, unless it's already loaded.
      observeEvent(restore_context(), {
        req(restore_context()$stage == 'loading_data')
        context <- restore_context()

        if (!is_assay_already_loaded(current, context$project, context$analysis)) {
          load_assay(
            context$project,
            context$analysis,
            context$allowed_path,
            reset_tabs = FALSE
          )
        }

        showModal(bookmark_loading_modal())
        context$stage <- 'applying_inputs'
        restore_context(context)
      })

      # Applies the saved inputs for the first time, then starts the
      # browser round trips.
      observeEvent(restore_context(), {
        req(restore_context()$stage == 'applying_inputs')
        context <- restore_context()

        updateSelectizeInput(
          root_session,
          'gene.to.plot',
          choices = gene_choices(),
          selected = context$inputs$gene.to.plot,
          server = TRUE
        )

        apply_restored_inputs(context, root_session)
        context$stage <- 'stabilizing'
        context$round <- 1L
        restore_context(context)
        send_bookmark_restore_probe(session, root_session, context$round)
      })

      # Each ack means the browser has rendered the previous flush. Ignores acks
      # that are stale, duplicated, or arrive when no restore is running.
      observeEvent(
        input$bookmark_restore_ack,
        {
          context <- restore_context()
          ack <- input$bookmark_restore_ack

          if (
            !identical(context$stage, 'stabilizing') ||
              !is.list(ack) ||
              !isTRUE(ack$round == context$round)
          ) {
            return()
          }

          # Done: tell plot modules to refresh and end the restore.
          if (context$round >= total_restore_rounds) {
            bookmark_restore_tick(bookmark_restore_tick() + 1L)
            restore_context(NULL)
            removeModal()
            return()
          }

          # Re-apply so inputs whose choices depend on earlier inputs can take
          # their saved values.
          if (context$round %% rounds_per_apply == 0) {
            apply_restored_inputs(context, root_session)
          }

          context$round <- context$round + 1L
          restore_context(context)
          send_bookmark_restore_probe(session, root_session, context$round)
        },
        ignoreInit = TRUE
      )

      ##### Save/Delete Actions #####

      can_modify_selected <- reactive(
        can_modify_bookmark(selected_bookmark(), username(), is_admin())
      )

      # reactive and output with info about admin user
      is_admin <- reactive(isTRUE(settings()$is_admin))
      output$can_modify_selected <- can_modify_selected
      outputOptions(output, 'can_modify_selected', suspendWhenHidden = FALSE)

      # output flag used to toggle bookmark table if bookmarks present
      output$has_bookmarks <- reactive({
        req(root_input$dds, root_input$assay, assay.list$l)
        bookmarks <- selected_assay_bookmarks()
        nrow(bookmarks) > 0
      })
      outputOptions(output, 'has_bookmarks', suspendWhenHidden = FALSE)


      # Save Bookmark Modal Machinery
      observeEvent(input$save_bookmarks_desc_rows_selected, {
        row <- input$save_bookmarks_desc_rows_selected
        updateTextInput(session, 'bookmark_name',
          value = current_assay_bookmarks()$label[row])
      })

      # Show bookmarks for what is currently loaded. As opposed
      # to selected_assay_bookmarks which shows bookmarks for what is
      # selected in the bookmarks description table.
      current_assay_bookmarks <- reactive({
        saved_list_tick()
        list_saved_bookmarks(
          get_bookmarks_path(current$proj, current$analysis, bookmarks_root),
          current$proj,
          current$analysis)
      })

      output$save_bookmarks_desc <- renderDT({
        bookmarks_df <- current_assay_bookmarks()

        hidden_columns <- match(
          c('id', 'comp', 'mode'),
          names(bookmarks_df)
        ) - 1

        selected <- which(bookmarks_df$id %in% loaded_bookmark()$id)

        datatable(
          bookmarks_df,
          rownames = FALSE,
          selection = list(mode = 'single', selected = selected),
          caption = tags$caption(
            style = 'font-weight: bold; font-size 15px;',
              'Existing bookmarks for loaded analysis'
          ),
          options = list(
            dom = 'tp',
            columnDefs = list(
              list(
                visible = FALSE,
                targets = hidden_columns
              )
            ),
            displayStart = if (length(selected)) (selected - 1) %/% 10 * 10 else 0
          )
        )
      })


      # bookmark description table displayed in the saved-bookmarks UI.
      output$bookmarks_desc <- renderDT({
        bookmarks_df <- selected_assay_bookmarks()

        hidden_columns <- match(
          c('id', 'comp', 'mode'),
          names(bookmarks_df)
        ) - 1

        selected <- which(bookmarks_df$id %in% loaded_bookmark()$id)

        datatable(
          bookmarks_df,
          rownames = FALSE,
          selection = list(mode = 'single', selected = selected),
          caption = tags$caption(
            style = 'font-weight: bold; font-size: 15px;',
            'Bookmarks for selected analysis'
          ),
          options = list(
            dom = 'tp',
            columnDefs = list(
              list(
                visible = FALSE,
                targets = hidden_columns
              )
            ),
            displayStart = if (length(selected)) (selected - 1) %/% 10 * 10 else 0
          )
        )
      })

      selected_bookmark <- reactive({
        row <- input$bookmarks_desc_rows_selected
        bookmarks <- selected_assay_bookmarks()
        if (length(row) != 1 || row > nrow(bookmarks)) return(NULL)
        bookmarks[row, , drop = FALSE]
      })

      observeEvent(input$show_bookmark_url, {
        id <- selected_bookmark()$id

        if (is.null(id)) {
          showNotification('Please select a bookmark first.')
          return()
        }
        client_data <- session$clientData
        url <- paste0(client_data$url_protocol, '//', client_data$url_hostname,
          if (nzchar(client_data$url_port)) paste0(':', client_data$url_port),
          client_data$url_pathname, '?_state_id_=', id
        )
        showModal(show_bookmark_url_modal(url))
      })

      # Loads a bookmark picked from the UI. Restores in-session when the
      # bookmarked assay is already loaded, otherwise reloads via the URL path.
      # No access checks here: both paths go through start_bookmark_restore().
      observeEvent(input$go_bookmark, {
        id <- selected_bookmark()$id
        req(id, is.null(restore_context()))

        values <- read_bookmark_values(
          find_bookmark_dir(id, selected_bookmarks_dir())
        )

        updateQueryString(paste0('?_state_id_=', id), mode = 'push')

        if (!is.null(values) &&
            is_assay_already_loaded(current, values$project, values$assay_analysis)) {
          start_bookmark_restore(values, id)
        } else {
          session$reload()
        }
      })

      observeEvent(input$delete_bookmark, {
        req(selected_bookmark()$id)
        showModal(bookmark_delete_modal(session$ns))
      })

      # Deleting is a separate rule from restoring: admin only, and limited to
      # the selected assay's bookmark folder.
      observeEvent(input$delete_bookmark_do, {
        id <- selected_bookmark()$id
        req(id, can_modify_selected())


        removeModal()
        if (!delete_saved_bookmark(id, selected_bookmarks_dir())) {
          showNotification('Unable to delete bookmark.', type = 'error')
          return()
        }
        saved_list_tick(saved_list_tick() + 1L)
        showNotification('Bookmark deleted.')
      })

      ##### Loaded Bookmark Label #####

      # Label for the 'currently loaded' card (returned to app.R). Only shown
      # while the loaded assay is still the one the bookmark belongs to, since
      # loaded_bookmark isn't cleared when another assay is loaded.
      loaded_bookmark_label <- reactive({
        bookmark <- loaded_bookmark()
        if (!is.null(bookmark) &&
            identical(bookmark$project, current$proj) &&
            identical(bookmark$analysis, current$analysis)) bookmark$label else NULL
      })

      helpButtonServer('bookmarking_help')

      ##### Returned Values #####

      # Reactives used by other parts of the dashboard.
      list(
        is_restoring = reactive(!is.null(restore_context())),
        bookmark_restore_tick = bookmark_restore_tick,
        loaded_bookmark_label = loaded_bookmark_label
      )
    }
  )
}


######## UI #######################

#' @rdname bookmarkingmod
#' @export
bookmarkButtonUI <- function(id) {
  ns <- NS(id)

  actionButton(
    ns('bookmark'),
    label = 'Bookmark',
    class = 'btn-primary'
  )
}


#' @rdname bookmarkingmod
#' @export
bookmarksDescriptionUI <- function(id) {
  ns <- NS(id)
  tagList(
    DTOutput(ns('bookmarks_desc')),
    br(),
    conditionalPanel('output.has_bookmarks', ns = ns,
      fluidRow(style='margin-bottom: 10px',
            column(1,
              actionButton(ns('go_bookmark'), 'Load', class = 'btn-primary')
            ),
            column(1,
              conditionalPanel(
                'output.can_modify_selected', ns = ns,
                actionButton(ns('delete_bookmark'), 'Delete', class = 'btn-primary')
              )
            ), # column
            column(7, span()),
            column(1,
              actionButton(ns('show_bookmark_url'), 'Show URL', class = 'btn-primary'),
            ),
            column(1, offset = 1,
              helpButtonUI(ns('bookmarking_help'))
            ), # column
      ) # fluidRow
    ) # conditionalPanel
  ) # tagList
}

##### Modals #####

bookmark_name_modal <- function(session, loaded_bookmark) {
  modalDialog(
    title = 'Save current view',
    textInput(
        session$ns('bookmark_name'),
        'Bookmark name',
        value = loaded_bookmark,
        placeholder = 'Enter a descriptive name',
        width = '100%'
    ),
    DTOutput(session$ns('save_bookmarks_desc')),
    selectInput(
      session$ns('overwrite'),
      label='Force overwrite?',
      choices=c(FALSE, TRUE)
    ),
    footer = tagList(
      modalButton('Cancel'),
      actionButton(
        session$ns('confirm_bookmark'),
        'Save Bookmark',
        class = 'btn-primary'
      )
    ),
    easyClose = TRUE,
    size = 'l'
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

show_bookmark_url_modal <- function(url) {
  modalDialog(
    fluidRow(
      column(12,
        paste0("Share with: "),
        br(),
        strong(url),
        br(),
        em('The link will only load for people who have access to the data.')
      ) # column
    ), # fluidRow
    title = 'Bookmark Link',
    easyClose = TRUE
  )
}
