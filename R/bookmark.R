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

      onBookmark(
        function(state) {
          # Stores metadata that isnt well captured by the inputs
          # snapshot.
          state$values$project <- current$proj
          state$values$assay_analysis <- current$analysis
          state$values$label <- bookmark_label()
          state$values$mode <- root_input$mode
          state$values$comp <- paste(root_input$comp_all, collapse = ', ')
          state$values$saved_at <- as.character(Sys.time())

          state$values$inputs <- save_context()
        },
        session = root_session
      )

      # Runs after the bookmark has been generated.
      onBookmarked(
        function(url) {
          updateQueryString(url)
          refresh_bookmarks_tick(
            isolate(refresh_bookmarks_tick()) + 1L
          )
          save_context(NULL)
        },
        session = root_session
      )

      # Runs a the end of the initial session load. Captures
      # the server values loaded by the built-in portion of shiny's
      # bookmarking so it can be applied and stablized
      # after auth check and data is loaded.
      onRestored(
        function(state) {
          allowed_path <- resolve_bookmark_access(
            assay.list$l,
            state$values$project,
            state$values$assay_analysis
          )

          if (is.null(allowed_path)) {
            showModal(
              modalDialog(
                title='You do not have access to this data area.',
                easyClose = FALSE
              )
            )
            return()
          }

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

        load_assay(
          context$project,
          context$analysis,
          context$allowed_path,
          reset_tabs = FALSE
        )

        showModal(
          modalDialog(
            icon('spinner', class = 'fa-spin'),
            'Loading bookmark, please wait',
            footer = NULL,
            easyClose = FALSE
          )
        )

        context <- restore_context()
        context$stage <- 'replaying_initial_inputs'
        restore_context(context)
      })

      # After the next server flush send a message asking
      # the browser to confirm it has reached a render frame.
      send_restore_probe <- function(round) {
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
          send_restore_probe(
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

          send_restore_probe(next_round)
        },
        ignoreInit = TRUE,
        priority = -100
      )

      ########### Bookmark loading UI ###############

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
                session$ns('bookmark_name'),
                'Bookmark name',
                placeholder = 'Enter a descriptive name'
              ),
              footer = tagList(
                modalButton('Cancel'),
                actionButton(
                  session$ns('confirm_bookmark'),
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

          # Prevents inputs from being saved by shiny bookmarking as default.
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

      ui_selected_assay <- reactive({
        req(root_input$dds, root_input$assay)

        project_assays <- assay.list$l[[root_input$dds]]
        idx <- match(root_input$assay, unname(project_assays))

        req(!is.na(idx))

        list(
          project = root_input$dds,
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
        ) -
          1

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
          return(NULL)
        }

        choices <- setNames(bookmarks_df$id, bookmarks_df$label)

        tagList(
          selectInput(
            session$ns('bookmark_dropdown'),
            'Saved Bookmarks',
            choices = choices
          ),
          actionButton(
            session$ns('go_bookmark'),
            'Load Bookmark',
            class = 'btn-primary'
          ),
          if (settings()$is_admin) {
            actionButton(
              session$ns('delete_bookmark'),
              'Delete Bookmark',
              class = 'btn-primary'
            )
          }
        )
      })

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
                session$ns('delete_bookmark_do'),
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

      # reactiveVals returned and used by other parts of the
      # dashboard.
      list(
        restore_context = restore_context,
        restore_tick = bookmark_restore_tick
      )
    }
  )
}
