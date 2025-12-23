# Get path to access yaml file

This function checks for an environment variable 'CARNATION_ACCESS_YAML'
to specify directory to save access yaml. If env variable does not exist
uses home directory as save location.

## Usage

``` r
get_access_path()
```

## Value

path to access yaml

## Examples

``` r
p <- get_access_path()
#> Environment variable "CARNATION_ACCESS_YAML" not found.Using default location to save access yaml:/home/runner
```
