# Create access yaml

This function creates an access yaml file. This is primarily intended
for the first run.

## Usage

``` r
create_access_yaml(user, user_group, data_area)
```

## Arguments

- user:

  User name

- user_group:

  User group

- data_area:

  Path to data area containing RDS files

## Value

Invisibly returns `NULL`. This function is primarily used for its side
effect of saving a yaml file with access settings

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
```
