#' Settings module ui
#'
#' This generates the settings tab that allows users to
#' add data areas and app user/user groups.
#'
#' @param id Input id
#' @param panel context for generating ui elements ('sidebar' or 'main')
#' @param username user name
#'
#' @export
settingsUI <- function(id, panel, username){
  ns <- NS(id)

  if(!is.null(username())){
    user_settings <- tagList(
                       conditionalPanel(paste0('input["', ns('access'), '"] == "User list"'),
                         div(tags$b('Edit user list'),
                             style='text-align: center; margin-bottom: 20px; font-size: 15px;'),

                         fluidRow(
                           column(4, align='center',
                             actionButton(ns('add_user'), 'Add user')
                           ), # column
                           column(4,
                             actionButton(ns('edit_user'), 'Edit user'),
                           ), # column
                           column(4,
                             actionButton(ns('rm_user'), 'Delete user')
                           )  # column
                         )  # fluidRow
                       ) # conditionalPanel
                     ) # tagList

    user_main <- tabPanel('User list',
                   fluidRow(
                     column(10,
                       DTOutput(ns('user_list'))
                     )  # column
                   )  # fluidRow
                 ) # tabPanel

    msg <- paste('Add/edit', packageName(), 'data & user settings here')
  } else {
    user_settings <- tagList()
    user_main <- tabPanel(title=NULL)
    msg <- paste('Add/edit', packageName(), 'data settings here')
  }

  if(panel == 'sidebar'){
    tag <- tagList(
              br(),
              br(),
              span(msg, style='font-size: 15px;'),
              br(),
              br()
            )
  } else if(panel == 'main'){
    tag <- tagList(
             sidebarPanel(width=4,

               # style the sidebar
               HTML(
                 '<style type="text/css">
                 .well { background-color: #FFFFFF; border-width: 0px 3px 3px 1px; }
                 </style>'
               ), # HTML

               user_settings,

               conditionalPanel(paste0('input["', ns('access'), '"] == "Data areas"'),
                 div(tags$b('Edit data areas'),
                     style='text-align: center; margin-bottom: 20px; font-size: 15px;'),

                 fluidRow(
                   column(12,
                     align='center',
                     style='margin-bottom: 10px;',
                     actionButton(ns('add_grp'), 'Add area')
                   ),
                   column(12,
                     align='center',
                     style='margin-bottom: 10px;',
                     actionButton(ns('edit_grp'), 'Edit area')
                   ),
                   column(12,
                     align='center',
                     style='margin-bottom: 10px;',
                     actionButton(ns('rm_grp'), 'Remove area')
                   )  # column
                 )  # fluidRow
               ), # conditionalPanel

               fluidRow(
                 column(5, align='center',
                   actionButton(ns('reset_access'), 'Reset',
                                style='margin-top: 20px')
                 ),  #column
                 column(7, align='center',
                   actionButton(ns('save_access'), 'Save changes',
                                style='margin-top: 20px;',
                                class='btn-primary')
                 )  # column
               )  # fluidRow

             ), # sidebarPanel

             mainPanel(
               tabsetPanel(id=ns('access'),
                 tabPanel('Data areas',
                   fluidRow(
                     column(10,
                       DTOutput(ns('data_areas'))
                     ), # column
                     column(1, align='left',
                       helpButtonUI(ns('settings_help'))
                     ) # column
                   )    # fluidRow
                 ), # tabPanel
                 user_main
               ) # tabsetPanel
             ) # mainPanel
    ) # tagList
  }
  return(tag)
}

#' Settings module server
#'
#' Server code for settings module
#'
#' @param id Input id
#' @param details reactive list with user name & app location details
#' @param depth project name depth
#' @param end_offset project name end offset
#' @param assay_fun function to parse assay names from file path
#' @param config reactive list with config settings
#'
#' @export
settingsServer <- function(id, details, depth, end_offset, assay_fun, config){
  moduleServer(
    id,

    function(input, output, session){

    ns <- NS(id)

    # get vars from config
    pattern <- reactive({ config()$server$pattern })
    admin_group <- reactive({ config()$server$admin_group })
    staging_dir <- reactive({ config()$server$staging_dir })

    access_yaml <- reactiveValues(l=NULL)
    username <- reactive({
      details()$username
    })

    observe({
      # check if access yaml is present
      p <- tryCatch(
             read_access_yaml(),
             error=function(e){ e }
           ) # tryCatch

      if(!all(names(p) %in% c('user_group', 'data_area'))){
        if(grepl('Access yaml not found', p$message)){
          if(!is.null(username())){
            user_inputs <- tagList(
                             textInput(ns('user_name'), label='User name',
                                       value=username()),
                             textInput(ns('user_group'), label='User group',
                                       value=admin_group())
                            )
          } else {
            user_inputs <- tagList()
          }

          showModal(
            modalDialog(
              tags$b('Access settings not configured', style='color: red;'),
              br(), br(),
              span('Please enter details below before first run'),
              br(),
              tagList(
                user_inputs,
                textInput(ns('data_area'), label='Data area',
                          value=NULL,
                          placeholder='Folder where data is stored')
              ),
              footer=tagList(
                       actionButton(ns('create_access_do'), 'OK')
                     )
            ) # modalDialog
          ) # showModal
        } else {
          showModal(
            modalDialog(
              tags$b('Error loading access settings', style='color: red;'),
              br(), br(),
              span(p$message),
              br(),
              span('Please check file permissions of "', get_access_path(), '" and try again'),
              footer=NULL
            ) # modalDialog
          ) # showModal
        }
      } else {
        access_yaml$l <- p
      }
    }) # observe

    observeEvent(input$create_access_do, {
      if(is.null(input$data_area)){
        showNotification(
          'Data area cannot be empty!', type='warning'
        )

        req(input$data_area)
      } else if(!dir.exists(input$data_area)){
        showNotification(
          'Data area does not exist on disk!', type='warning'
        )

        validate(
          need(dir.exists(input$data_area), '')
        )
      }

      if(is.null(username())){
        # check that user name/group is not empty
        create_access_yaml(config()$server$default_user,
                           admin_group(),
                           input$data_area)
      } else {
        # check that user name is not empty
        if(is.null(input$user_name)){
          showNotification(
            'User name cannot be empty!', type='warning'
          )

          req(input$user_name)
        }

        # check that user name is not empty
        if(is.null(input$user_group)){
          showNotification(
            'User group cannot be empty!', type='warning'
          )

          req(input$user_group)
        }

        create_access_yaml(input$user_name,
                           input$user_group,
                           input$data_area)
      }

      access_yaml$l <- read_access_yaml()

      removeModal()
    })

    # reactive to keep list of users and data areas
    user_access <- reactiveValues(orig=NULL, users=NULL, areas=NULL)
    reload_parent <- reactiveValues(flag=FALSE)

    observeEvent(access_yaml$l, {
      if(is.null(username())){
        lst <- access_yaml$l
      } else {
        lst <- check_user_access(access_yaml$l, username(),
                                 admin=admin_group())
      }

      user_access$orig <- lst
      user_access$users <- lst$user_group
      user_access$areas <- lst$data_area
    })

    # format user list table
    get_user_df <- function(u){
      user_df <- do.call('rbind',
                         lapply(u,
                                function(x) paste(x[order(x)], collapse=', '))
                         )
      user_df <- data.frame(user=rownames(user_df),
                            user_group=user_df[,1])
      user_df <- user_df[order(user_df$user_group),]

      user_df
    }

    ######### User list #########

    observeEvent(user_access$users, {
      validate(
        need(!is.null(user_access$users), 'No users added!')
      )
      all_users <- names(user_access$users)

      user_df <- do.call('rbind',
                         lapply(user_access$users,
                                function(x) paste(x[order(x)], collapse=', '))
                         )
      user_df <- cbind(user=rownames(user_df),
                       user_group=user_df[,1])

      output$user_list <- renderDT({
          get_user_df(user_access$users) %>%
              datatable(rownames=FALSE,
                                selection='single')
      })
    })

    # add user
    observeEvent(input$add_user, {
      showModal(
        modalDialog(
            title='Add user',
            fluidRow(
                column(6, 'Enter user name'),
                column(6,
                    textInput(ns('add_user_name'), label=NULL,
                              value='')
                    )  # column
                ),  # fluidRow
            fluidRow(
                column(6, 'Select user group(s)'),
                column(6,
                    selectizeInput(ns('add_user_grps'),
                            label=NULL,
                            choices=c(admin_group(), names(user_access$areas)),
                            multiple=TRUE
                        )  # selectizeInput
                    )  # column
                ),  # fluidRow
            span('You can choose multiple user_groups',
                 style='font-style: italic;'),
            footer=tagList(
                modalButton('Cancel'),
                actionButton(ns('add_user_do'), 'OK')
              ),  # tagList
            easyClose=TRUE
        )  # modalDialog
      )  # showModal
    })

    observeEvent(input$add_user_do, {
      df <- user_access$users
      g <- input$add_user_grps
      u <- input$add_user_name

      if(length(g) == 0 | u == '' | is.null(u)){
        showNotification('User name or user groups cannot be empty')
      }

      validate(
        need(g != '' & !is.null(u) & u != '', 'Waiting for data')
      )

      if(u == config()$server$default_user){
        showNotification('User name is reserved for internal use!', type='warning')
      }

      validate(
        need(u != config()$server$default_user, 'User name reserved for internal use')
      )

      if(!all(g %in% df[[u]])){
          df[[ u ]] <- unique(c(df[[u]], g))
      } else if(!u %in% names(df)){
          df[[ u ]] <- g
      }
      user_access$users <- df

      removeModal()
    })

    # remove user
    observeEvent(input$rm_user, {
      # get user from selected row (if any)
      if(!is.null(input$user_list_rows_selected)){
        idx <- input$user_list_rows_selected
        user_df <- get_user_df(user_access$users)
        selected <- user_df[idx, 'user']
      } else {
        selected <- NULL
      }

      showModal(
        modalDialog(
            title='Remove user',
            fluidRow(
                column(6,  'Select user'),
                column(6,
                    selectizeInput(ns('rm_user_name'),
                                label=NULL,
                                choices=c('Choose one', names(user_access$users)),
                                selected=selected
                      )  # selectizeInput
                    )  # column
                ),  # fluidRow
            footer=tagList(
                modalButton('Cancel'),
                actionButton(ns('rm_user_do'), 'OK')
              ),  # tagList
            easyClose=TRUE
        )  # modalDialog
      )  # showModal
    })

    observeEvent(input$rm_user_do, {
      df <- user_access$users
      u <- input$rm_user_name

      if(u == '' | u == 'Choose one'){
        showNotification('User name cannot be empty')
      }

      validate(
        need(u != 'Choose one' & u != '' & u != username(), 'Waiting for data')
      )

      if(u == username()){
        showNotification('Cannot remove current user')
      } else {
          df[[ u ]] <- NULL
          user_access$users <- df
      }

      removeModal()
    })

    # edit user
    observeEvent(input$edit_user, {
        # get user from selected row (if any)
        if(!is.null(input$user_list_rows_selected)){
          idx <- input$user_list_rows_selected
          user_df <- get_user_df(user_access$users)

          selected_user <- user_df[idx, 'user']
          selected_grp <- user_df[idx, 'user_group']
        } else {
          selected_user <- NULL
          selected_grp <- NULL
        }

        showModal(
            modalDialog(
                title='Edit user',
                fluidRow(
                    column(6,  'Select user'),
                    column(6,
                        selectizeInput(ns('edit_user_name'),
                                    label=NULL,
                                    choices=c('Choose one', names(user_access$users)),
                                    selected=selected_user
                      )  # selectizeInput
                    )  # column
                ),  # fluidRow
                fluidRow(
                    column(6,  'Select user groups'),
                    column(6,
                        selectizeInput(ns('edit_user_grps'),
                                       label=NULL,
                                       choices=NULL,
                                       selected=selected_grp,
                                       multiple=TRUE,
                                       options=list(create=TRUE)
                        )  # selectizeInput
                    )  # column
                ),  # fluidRow
                span('You can choose multiple user_groups',
                     style='font-style: italic;'),
                footer=tagList(
                    modalButton('Cancel'),
                    actionButton(ns('edit_user_do'), 'OK')
                  ),
                easyClose=TRUE
            )  # modalDialog
        )  # showModal
    })

    observeEvent(c(input$edit_user_name, user_access$areas), {
        validate(
            need(input$edit_user_name != 'Choose one',
                 '')
        )

        user_groups <- user_access$users[[ input$edit_user_name ]]
        all_groups <- unique(c(admin_group(), names(user_access$areas)))

        updateSelectizeInput(session,
                             'edit_user_grps',
                             choices=all_groups,
                             selected=user_groups)
    })

    observeEvent(input$edit_user_do, {
        validate(
            need(input$edit_user_name != 'Choose one',
                 '')
        )

        user_access$users[[ input$edit_user_name ]] <- input$edit_user_grps
        removeModal()
    })

    ######### Data areas ##########

    observeEvent(user_access$areas, {

      output$data_areas <- renderDT({
          get_area_df(user_access$areas) %>%
              datatable(rownames=FALSE,
                        selection='single')
      })

    })

    # format data areas table
    get_area_df <- function(u){
      area_df <-  do.call('rbind',
                         lapply(u, function(x) paste(x, collapse=', '))
                         )

      area_df <- data.frame(user_group=rownames(area_df),
                            data_area=area_df[,1])
      if(length(area_df$user_group) > 0)
        area_df <- area_df[order(area_df$user_group),]

      area_df
    }

    # add user group
    observeEvent(input$add_grp, {
      showModal(
        modalDialog(
            title='Add data group',
            fluidRow(
                column(6, 'Enter data group name'),
                column(6,
                    textInput(ns('add_grp_name'), label=NULL,
                              value='')
                    )  # column
                ),  # fluidRow
            fluidRow(
                column(6, 'Enter data area(s)'),
                column(6,
                    selectizeInput(ns('add_grp_areas'),
                            label=NULL,
                            choices=c('', unname(unlist(user_access$areas))),
                            multiple=TRUE,
                            options=list(create=TRUE)
                        )  # selectizeInput
                    )  # column
                ),  # fluidRow
            span('You can enter multiple data areas',
                 style='font-style: italic;'),
            footer=tagList(
                modalButton('Cancel'),
                actionButton(ns('add_grp_do'), 'OK')
              ),  # tagList
            easyClose=TRUE
        )  # modalDialog
      )  # showModal
    })

    observeEvent(input$add_grp_do, {
      df <- user_access$areas
      g <- input$add_grp_name
      a <- input$add_grp_areas

      if(length(a) == 0 | g == '' | is.null(g)){
        showNotification('User group name or data areas cannot be empty')
      }

      validate(
        need(!is.null(g) & g != '' & !is.null(a), 'Waiting for data')
      )

      if(!all(file.exists(a))){
          showNotification('Error: Data area does not exist on file', type='error')

          # this removes entered path from options
          updateSelectizeInput(session,
                               'add_grp_areas',
                               choices=c('', unname(unlist(user_access$areas)))
            )  # selectizeInput
      } else if(g %in% names(df)){
          showNotification('Error: user group already exists. Try "Edit user group" instead or choose different group name')
      }

      # check that path exists
      validate(
        need(all(file.exists(a)) & !g %in% names(df), '')
      )

      df[[ g ]] <- a
      user_access$areas <- df

      removeModal()
    })

    # remove user group
    observeEvent(input$rm_grp, {
      if(!is.null(input$data_areas_rows_selected)){
        idx <- input$data_areas_rows_selected
        area_df <- get_area_df(user_access$areas)
        selected <- area_df[idx, 'user_group']
      } else {
        selected <- NULL
      }

      showModal(
        modalDialog(
            title='Remove user group',
            fluidRow(
                column(6,  'Select user group'),
                column(6,
                    selectizeInput(ns('rm_grp_name'),
                                label=NULL,
                                choices=c('Choose one', names(user_access$areas)),
                                selected=selected
                      )  # selectizeInput
                    )  # column
                ),  # fluidRow
            footer=tagList(
                modalButton('Cancel'),
                actionButton(ns('rm_grp_do'), 'OK')
              ),  # tagList
            easyClose=TRUE
        )  # modalDialog
      )  # showModal
    })

    observeEvent(input$rm_grp_do, {
      df <- user_access$areas
      g <- input$rm_grp_name

      if(length(g) == 1){
          if(g == '' | g == 'Choose one')
            showNotification('User group name cannot be empty')
      }

      validate(
        need(g != 'Choose one' & g != '', 'Waiting for data')
      )

      df <- df[!names(df) %in% g]
      user_access$areas <- df

      removeModal()
    })

    # edit user group
    observeEvent(input$edit_grp, {
        if(!is.null(input$data_areas_rows_selected)){
          idx <- input$data_areas_rows_selected
          area_df <- get_area_df(user_access$areas)
          selected_name <- area_df[idx, 'user_group']
          selected_area <- area_df[idx, 'data_area']
        } else {
          selected_name <- NULL
          selected_area <- NULL
        }

        showModal(
            modalDialog(
                title='Edit user group',
                fluidRow(
                    column(6,  'Select user group'),
                    column(6,
                        selectizeInput(ns('edit_grp_name'),
                                    label=NULL,
                                    choices=c('Choose one', names(user_access$areas)),
                                    selected=selected_name
                          )  # selectizeInput
                        )  # column
                    ),  # fluidRow
                fluidRow(
                    column(6,  'Select data areas'),
                    column(6,
                        selectizeInput(ns('edit_grp_areas'),
                                       label=NULL,
                                       choices=NULL,
                                       selected=selected_area,
                                       multiple=TRUE,
                                       options=list(create=TRUE)
                            )  # selectizeInput
                        )  # column
                    ),  # fluidRow
                span('You can choose multiple data areas',
                     style='font-style: italic;'),
                footer=tagList(
                    modalButton('Cancel'),
                    actionButton(ns('edit_grp_do'), 'OK')
                  ),
                easyClose=TRUE
            )  # modalDialog
        )  # showModal
    })

    observeEvent(c(input$edit_grp_name, input$edit_grp_do), {
        validate(
            need(input$edit_grp_name != 'Choose one',
                 '')
        )

        data_areas <- user_access$areas[[ input$edit_grp_name ]]
        all_areas <- unname(unique(unlist(user_access$areas)))

        updateSelectizeInput(session,
                             'edit_grp_areas',
                             choices=all_areas,
                             selected=data_areas)
    })

    observeEvent(input$edit_grp_do, {
        validate(
            need(input$edit_grp_name != 'Choose one',
                 '')
        )

        if(is.null(input$edit_grp_areas)){
            showModal(
              modalDialog(
                tags$b('Error: Data areas cannot be empty',
                       style='color: red;'),
                br(), br(),
                span('To delete user group, please use "Remove user group" button'),
                footer=modalButton('OK'),
                easyClose=TRUE
              )
            )
        } else {
            # check that data areas exist
            check_exists <- unlist(lapply(input$edit_grp_areas, dir.exists))
            if(any(!check_exists)){
              showNotification(
                paste('Error: Data areas must exist on disk! Following folders not found:',
                      paste(input$edit_grp_areas[!check_exists], collapse=', ')),
                type='error'
              )
            }
            validate(
              need(!any(!check_exists), '')
            )

            user_access$areas[[ input$edit_grp_name ]] <- input$edit_grp_areas
            removeModal()
        }

    })

    # save access changes
    observeEvent(input$save_access, {
        # check if user is admin
        if(is.null(username())) u <- config()$server$default_user
        else u <- username()
        is_admin <- admin_group() %in% user_access$orig$user_group[[ u ]]

        if(is_admin){
            showModal(
                modalDialog(
                    div(tags$b('This will save any changes to access settings & cannot be undone. Proceed?', style='color: red;')),
                    footer=tagList(
                        modalButton('Cancel'),
                        actionButton(ns('save_access_do'), 'Yes')
                    )  # tagList
                )  # modalDialog
            )  # showModal
        } else {
            reset_fun()
            showModal(
                modalDialog(
                    div(tags$b(
                      paste0('Forbidden: You must be member of "',
                             admin_group(), '" group to edit access settings!'),
                      style='color: red;'
                      ) # tags
                    ), # div
                    footer=tagList(
                        modalButton('OK')
                    )  # tagList
                )  # modalDialog
            )  # showModal
        }
    })

    observeEvent(input$save_access_do, {
      lst <- list(user_group=user_access$users,
                  data_area=user_access$areas)
      save_access_yaml(lst)

      reload_parent$flag <- TRUE

      removeModal()
    })

    # reset access changes
    observeEvent(input$reset_access, {
        showModal(
            modalDialog(
                div(tags$b('This will reset any changes to user access settings. Proceed?', style='color: red;')),
                footer=tagList(
                    modalButton('Cancel'),
                    actionButton(ns('reset_access_do'), 'Yes')
                )  # tagList
            )  # modalDialog
        )  # showModal
    })

    reset_fun <- function(){
      lst <- user_access$orig
      user_access$users <- lst$user_group
      user_access$areas <- lst$data_area
    }

    observeEvent(input$reset_access_do, {
      reset_fun()

      removeModal()
    })

    no_projects_modal <- function(){
      modalDialog(
        div(tags$b('No projects found!', style='color: red;')),
        br(),
        span('To add/edit a data area go to the "Settings" tab'),
        footer=modalButton('OK'),
        easyClose=TRUE
      ) # modalDialog
    }

    # update object menu
    settings <- eventReactive(c(access_yaml$l,
                                reload_parent$flag), {
      validate(
        need(!is.null(access_yaml$l), '')
      )

      if(is.null(username())) u <- config()$server$default_user
      else u <- username()
      is_admin <- admin_group() %in% user_access$orig$user_group[[ u ]]

      if(is.null(username())){
        d <- check_user_access(access_yaml$l,
                               config()$server$default_user,
                               admin=admin_group())
      } else {
        d <- check_user_access(access_yaml$l,
                               username(),
                               admin=admin_group())
      }

      # var to detect if in shinymanager admin view
      if(is.null(details()$where)) shinyadmin <- FALSE
      else if(details()$where != 'admin') shinyadmin <- FALSE
      else shinyadmin <- TRUE

      if(is.null(d)){
        # single-user mode
        if(is.null(u)){
          showModal(
            no_projects_modal()
          )
        } else {
          # don't show if in shinymanager admin view
          if(!shinyadmin){
            if(!is_admin){
              showModal(
                modalDialog(
                  div(tags$b('No access permissions!', style='color: red;')),
                  br(),
                  span('Please contact site administrators to enable access'),
                  footer=modalButton('OK')
                )
              )
            } else {
              showModal(
                no_projects_modal()
              )
            }
          }
        }
      }

      validate(
          need(!is.null(d), 'No access permissions')
      )

      l <- unlist(lapply(unique(d$data_area),
              function(x) list.files(x,
                              pattern=paste0(pattern(), '\\.rds$'),
                              ignore.case=TRUE,
                              recursive=TRUE,
                              full.names = TRUE)))

      if(is.null(l) | length(l) == 0){
        # don't show modal if in admin panel
        if(!shinyadmin){
          showModal(
            no_projects_modal()
          ) # showModal
        }

        ll <- list(assay_list=NULL,
                   reload_parent=reload_parent$flag,
                   is_admin=is_admin)
        return(ll)

      }

      # get project name
      proj <- unlist(lapply(l, function(x){
                get_project_name_from_path(x,
                                           depth=depth,
                                           end_offset=end_offset,
                                           staging_dir=staging_dir())
              }))

      # get assay name
      assays <- unlist(lapply(l, assay_fun))

      # list of projects, with associated assays
      alist <- list()
      for(p in unique(proj)){
        idx <- which(proj == p)
        alist[[p]] <- l[idx]
        names(alist[[p]]) <- assays[idx]
      }

      # if not admin, filter out 'dev' datasets
      # else, fix project & assay labels

      # for a dev dataset in pi/project/dev/main.rnaseq.rds
      # the project name has to be changed from
      # pi/project/dev (old_label) to pi/project (fixed_label)
      #
      # and the assay name has to be changed from
      # main to dev/main

      # get list of dev datasets
      dev_idx <- which(basename(names(alist)) == staging_dir())
      if(length(dev_idx) > 0){
        old_labels <- NULL
        if(!is_admin){
          # if not admin, remove dev data from assay list
          alist <- alist[-dev_idx]
        } else {
          # if admin, tweak project & assay name
          for(i in dev_idx){
              old_label <- names(alist)[i]
              fixed_label <- dirname(old_label)

              # add dev/ prefix to assay label for dev datasets
              names(alist[[old_label]]) <- paste0(staging_dir(),
                                                  '/',
                                                  names(alist[[old_label]]))

              # if fixed label already exists, add to existing list
              # otherwise create label
              if(fixed_label %in% names(alist)){
                  # if assay label is duplicated, skip dataset
                  if(any(names(alist[[old_label]]) %in% names(alist[[fixed_label]]))){
                      showNotification(
                        paste0('Duplicate assay label for dev dataset in "',
                               fixed_label, '" ... Skipping'),
                        type='warning')
                  } else {
                      alist[[fixed_label]] <- c(alist[[fixed_label]], alist[[old_label]])
                  }
              } else {
                  alist[[fixed_label]] <- alist[[old_label]]
              }
              old_labels <- c(old_labels, old_label)
          }
        }
        alist <- alist[!names(alist) %in% old_labels]
      }

      assay_list <- alist

      # find and read project descriptions
      project_descriptions <- list()
      for(name in names(assay_list)){
        # First parse project directory
        # 1. find match for project name
        # 2. get substring from assay_choices
        mm <- regexpr(name, assay_list[[name]][1])
        pd_path <- file.path(
                     substr(assay_list[[name]][1], 1, mm + attr(mm, 'match.length') - 1),
                     'project-description.yaml')

        # read project descriptions if file exists
        if(file.exists(pd_path)){
          tmp_desc <- read_yaml(pd_path)

          # if not admin, filter out staged data
          if(!is_admin){
            idx <- grep(staging_dir(), names(tmp_desc))
            tmp_desc <- tmp_desc[ -idx ]
          }

          project_descriptions[[ name ]] <- tmp_desc
        }
      }

      list(assay_list=assay_list,
           reload_parent=reload_parent$flag,
           is_admin=is_admin,
           project_descriptions=project_descriptions)
    })

    helpButtonServer('settings_help', size='l')

    return(reactive(settings()))

    }   # function
  ) # moduleServer
}

