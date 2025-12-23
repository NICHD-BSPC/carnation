# Save access yaml to file

This function saves access details (user groups and data areas) to the
designated access yaml file.

## Usage

``` r
save_access_yaml(lst)
```

## Arguments

- lst:

  list of data frames with user_groups and data_areas

## Value

save access settings to yaml file

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

# read access yaml
lst <- read_access_yaml()
#> Environment variable "CARNATION_ACCESS_YAML" not found.Using default location to save access yaml:/home/runner

# add new user
lst$user_group$admin <- c(lst$user_group$admin, 'user1')

# save to access settings
save_access_yaml(lst)
#> Environment variable "CARNATION_ACCESS_YAML" not found.Using default location to save access yaml:/home/runner
```
