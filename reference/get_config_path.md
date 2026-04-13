# Get path to local config yaml file

This function checks for an environment variable `CARNATION_CONFIG_YAML`
to specify the local config yaml path. If the variable is not set, a
default path in the home directory is used.

## Usage

``` r
get_config_path()
```

## Value

path to local config yaml

## Examples

``` r
p <- get_config_path()
#> Error in get_config_path(): could not find function "get_config_path"
```
