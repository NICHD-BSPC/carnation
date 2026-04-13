# Get config

This function reads the bundled package config and returns it. If a
local config yaml exists, only supported user-editable settings are
merged into the returned config.

## Usage

``` r
get_config(config_path = NULL)
```

## Arguments

- config_path:

  optional path to a local config yaml. If `NULL`, uses the path
  returned by
  [`get_config_path()`](https://nichd-bspc.github.io/carnation/reference/get_config_path.md).

## Value

list containing config items

## Examples

``` r
cfg <- get_config()
```
