# Get data areas a user has access to

This function takes a username and returns a list with two elements:

## Usage

``` r
check_user_access(al, u, admin = "admin")
```

## Arguments

- al:

  list with access settings; should have two elements - user_group &
  data_area

- u:

  user name

- admin:

  Admin user group

## Value

list of user groups and data areas

## Details

user_group: one element vector data_area: vector of data areas

## Examples

``` r
# save access details to file
home <- Sys.getenv('HOME')

# create carnation data area if it doesn't exist
carnation_home <- file.path(home, 'carnation/data')
if(!dir.exists(carnation_home)) dir.create(carnation_home)
#> Warning: cannot create dir '/home/runner/carnation/data', reason 'No such file or directory'

create_access_yaml(user = 'admin',
                   user_group = 'admin',
                   data_area = carnation_home)
#> Environment variable "CARNATION_ACCESS_YAML" not found.Using default location to save access yaml:/home/runner

# get current user access details
al <- read_access_yaml()
#> Environment variable "CARNATION_ACCESS_YAML" not found.Using default location to save access yaml:/home/runner

lst <- check_user_access(al, u='admin')
```
