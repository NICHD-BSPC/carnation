#' Get path to access yaml file
#'
#' This function checks for an environment variable 'CARNATION_ACCESS_YAML'
#' to specify directory to save access yaml. If env variable does not exist
#' uses home directory as save location.
#'
#' @return path to access yaml
#'
#' @export
get_access_path <- function(){
  if(Sys.getenv('CARNATION_ACCESS_YAML') != ''){
    path <- Sys.getenv('CARNATION_ACCESS_YAML')
    if(!dir.exists(path)){
      stop(
        paste('Environment variable "CARNATION_ACCESS_YAML" exists',
              'but specified location does not exist on disk:', path)
      )
    }
  } else {
    path <- path.expand('~')
    message(
      paste('Environment variable "CARNATION_ACCESS_YAML" not found.',
            'Using default location to save access yaml:', path)
    )
  }
  file.path(path, '.carnation-access.yaml')
}

#' Get data areas a user has access to
#'
#' This function takes a username and returns a
#' list with two elements:
#'
#' user_group: one element vector
#' data_area: vector of data areas
#'
#' @param al list with access settings; should have two elements - user_group & data_area
#' @param u user name
#' @param admin Admin user group
#'
#' @export
check_user_access <- function(al, u, admin='admin'){

  # lab of user
  idx <- which(names(al$user_group) == u)
  if(length(idx) == 0){
    return(NULL)
  } else {
    user_group <- al$user_group[idx]

    # if admin, give access to everything
    if(admin %in% unlist(user_group)){
      data_area <- al$data_area
    } else {
      idx <- which(names(al$data_area) %in% unlist(user_group))
      if(length(idx) == 0){
        return(NULL)
      } else {
        data_area <- al$data_area[ idx ]
      }
    }
    ll <- list(user_group=user_group, data_area=data_area)
  }

  return(ll)
}


