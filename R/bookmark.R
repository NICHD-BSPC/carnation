bookmarkingServer <- function(
  id,
  root_input,
  assay.list,
  load_assay,
  current,
  bookmarks_root,
  settings,
  restore_applies = 3L,
  rounds_per_apply = 2L
) {
  moduleServer(
    id,
    function(input, output, session) {
      ##### Shared State #####

      root_session <- session$rootScope()
      # NULL during normal operation. During bookmark restoration, this holds
      # the data needed across Shiny's restore functions.
      restore_context <- reactiveVal(NULL)

      # The stages held in restore_context$stage are 'load_data',
      # 'replaying_initial_inputs', 'stabilizing_restore', then back to NULL


      # On the bookmark saving side of this process, stores all the inputs so
      # that they can be excluded from being saved into the inputs.RDS object
      # and therefore from being automatically reapplied on restore. They are
      # saved into the state values object instead so we can control what is
      # sent to the browser as we load them manually after auth.
      save_context <- reactiveVal(NULL)

      # Is returned into the main app and fed into plotting modules that use a
      # refresh button to gate updating their plots or have potential race conditions
      # when all of their inputs are updated at once. This is used as a reactive
      # dependency and updated/invalidated after a server flushing loop.
      bookmark_restore_tick <- reactiveVal(0L)

      # Is used to repopulate saved bookmark lists when the assay is changed
      # or one is deleted.
      refresh_bookmarks_tick <- reactiveVal(0L)

      # Used to hold the bookmark label when saving.
      bookmark_label <- reactiveVal(NULL)

      # Holds the label and assay of the bookmark restored this session.
      loaded_bookmark <- reactiveVal(NULL)

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
          if (!is.null(isolate(restore_context()))) {
            return()
          }

          showModal(bookmark_name_modal(session$ns))
        },
        ignoreInit = TRUE
      )

      # Is triggered by the bookmark modal and triggers the main cascade.
      observeEvent(
        input$confirm_bookmark,
        {
          req(current$proj, current$analysis)

          label_validation <- validate_new_bookmark_label(
            input$bookmark_name,
            current$proj,
            current$analysis,
            bookmarks_root
          )

          if (!label_validation$valid) {
            showNotification(
              label_validation$message,
              type = 'error'
            )
            return()
          }

          bookmark_label(label_validation$label)
          removeModal()

          # Prevents inputs from being saved by shiny bookmarking in the default fashion.
          save_context(isolate(reactiveValuesToList(
            root_input,
            all.names = TRUE
          )))
          input_names <- names(save_context())
          setBookmarkExclude(
            input_names,
            session = root_session
          )

          root_session$doBookmark()
        },
        ignoreInit = TRUE
      )


      onBookmark(
        function(state) {
          # Stores metadata that isnt well captured by the inputs
          # snapshot.
          state$values$project <- current$proj
          state$values$assay_analysis <- current$analysis
          state$values$label <- bookmark_label()
          state$values$mode <- root_input$mode
          state$values$comp <- paste(root_input$comp_all, collapse = ', ')
          state$values$saved_at <- format(Sys.time(), '%Y-%m-%d %H:%M:%S')

          state$values$inputs <- save_context()
        },
        session = root_session
      )

      # Runs after the bookmark has been generated.
      onBookmarked(
        function(url) {
          saved_label <- bookmark_label()

          updateQueryString(url)
          loaded_bookmark(list(
            label = saved_label,
            project = isolate(current$proj),
            analysis = isolate(current$analysis)
          ))
          refresh_bookmarks_tick(
            isolate(refresh_bookmarks_tick()) + 1L
          )
          showNotification(
            paste0('Bookmark "', saved_label, '" saved.'),
            type = 'message'
          )
          save_context(NULL)
          bookmark_label(NULL)
        },
        session = root_session
      )

      ##### Restoring Bookmarks #####

      # Runs a the end of the initial session load. Captures
      # the server values loaded by the built-in portion of shiny's
      # bookmarking so it can be applied and stablized
      # after auth check and data is loaded.
      onRestored(
        function(state) {
          # Empty values mean load.interface couldn't find the bookmark.
          if (length(state$values) == 0) {
            showNotification('Bookmark is unavailable.', type = 'warning')
            return()
          }

          allowed_path <- resolve_bookmark_access(
            assay.list$l,
            state$values$project,
            state$values$assay_analysis
          )

          if (is.null(allowed_path)) {
            showModal(bookmark_access_denied_modal())
            return()
          }

          loaded_bookmark(list(
            label = state$values$label,
            project = state$values$project,
            analysis = state$values$assay_analysis
          ))

          restore_context(list(
            stage = 'load_data',
            project = state$values$project,
            analysis = state$values$assay_analysis,
            allowed_path = allowed_path,
            inputs = state$values$inputs,
            restore_loop_i = 0L
          ))
        },
        session = root_session
      )

      # Data loading phase.
      observe({
        req(restore_context()$stage ==  'load_data')

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

        context <- restore_context()
        context$stage <- 'replaying_initial_inputs'
        restore_context(context)
      })

      # Reapplies the stored inputs for the first time after the data has been reloaded,
      # then launches the flushing and input reapplying restore loop.
      observeEvent(
        restore_context(),
        {
          req(restore_context()$stage == 'replaying_initial_inputs')
          apply_restored_inputs(
            restore_context()$inputs,
            root_session
          )
          context <- restore_context()
          context$stage <- 'stabilizing_restore'
          context$restore_loop_i <- 1L
          restore_context(context)
          send_bookmark_restore_probe(
            session,
            root_session,
            context$restore_loop_i
          )
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
              !isTRUE(ack$round == context$restore_loop_i)
          ) {
            return()
          }

          current_round <- as.integer(ack$round)

          # When the correct number of roundtrips is reached,
          # increment bookmark_restore tick to trigger all the
          # plot refreshes and other things that depend on it.
          # This completes the session restore
          restore_roundtrips <- (restore_applies - 1L) * rounds_per_apply

          if (current_round >= restore_roundtrips) {
            bookmark_restore_tick(
              isolate(bookmark_restore_tick()) + 1L
            )

            restore_context(NULL)
            removeModal()
            return()
          }

          next_round <- current_round + 1L
          context$restore_loop_i <- next_round
          restore_context(context)

          # This allows to the system to reach settings that are gated
          # by other settings being needed to be selected first by reapplying
          # inputs after allowing earlier settings to be selected first.
          # When and the number of times this happens may need to
          # be tuned.
          if (current_round %% rounds_per_apply == 0) {
            apply_restored_inputs(
              context$inputs,
              root_session
            )
          }

          send_bookmark_restore_probe(
            session,
            root_session,
            next_round
          )
        },
        ignoreInit = TRUE,
        priority = -100
      )

      ##### Managing Saved Bookmarks #####

      assay_selection <- reactive({
        req(root_input$dds, root_input$assay)

        selected <- resolve_selected_assay(
          assay.list$l,
          root_input$dds,
          root_input$assay
        )

        req(selected)
        selected
      })

      selected_assay_bookmarks <- reactive({
        selected <- assay_selection()
        refresh_bookmarks_tick()

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

      is_admin <- reactive(isTRUE(settings()$is_admin))

      # Label for the 'currently loaded' card. Only shown while the loaded
      # assay is still the one the bookmark was restored into.
      loaded_bookmark_label <- reactive({
        b <- loaded_bookmark()
        if (!is.null(b) &&
            identical(b$project, current$proj) &&
            identical(b$analysis, current$analysis)) b$label else NULL
      })

      # Label for the picker. Only used when the dropdowns point at the
      # assay the bookmark was restored into.
      picker_bookmark_label <- reactive({
        b <- loaded_bookmark()
        sel <- assay_selection()
        if (!is.null(b) &&
            identical(b$project, sel$project) &&
            identical(b$analysis, sel$analysis)) b$label else NULL
      })

      setup_saved_bookmarks_ui(
        input = input,
        output = output,
        session = session,
        bookmarks = selected_assay_bookmarks,
        is_admin = is_admin,
        loaded_label = picker_bookmark_label
      )

      observeEvent(input$go_bookmark, {
        req(input$bookmark_dropdown)
        updateQueryString(
          paste0('?_state_id_=', input$bookmark_dropdown),
          mode = 'push'
        )
        session$reload()
      })

      observeEvent(input$delete_bookmark, {
        req(input$bookmark_dropdown)
        showModal(bookmark_delete_modal(session$ns))
      })

      observeEvent(input$delete_bookmark_do, {
        req(input$bookmark_dropdown, is_admin())

        assay <- assay_selection()
        bookmarks_dir <- get_bookmarks_path(
          assay$project,
          assay$analysis,
          bookmarks_root
        )
        deleted <- delete_saved_bookmark(
          input$bookmark_dropdown,
          bookmarks_dir
        )

        removeModal()
        if (!deleted) {
          showNotification('Unable to delete bookmark.', type = 'error')
          return()
        }
        refresh_bookmarks_tick(isolate(refresh_bookmarks_tick()) + 1L)
        showNotification('Bookmark deleted.')
      })

      ##### Public Module Values #####

      # reactiveVals returned and used by other parts of the dashboard.
      list(
        restore_context = restore_context,
        restore_tick = bookmark_restore_tick,
        loaded_label = loaded_bookmark_label
      )
    }
  )
}
