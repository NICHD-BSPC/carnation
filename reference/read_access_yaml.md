# Read access yaml with user groups and data areas

This function reads the access yaml file and returns user groups and
data areas as a list of data frames.

## Usage

``` r
read_access_yaml()
```

## Value

return carnation access settings from yaml file

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

al <- read_access_yaml()
#> Environment variable "CARNATION_ACCESS_YAML" not found.Using default location to save access yaml:/home/runner
```
