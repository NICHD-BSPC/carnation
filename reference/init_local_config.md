# Initialize local config

This function copies the bundled package config to a user-writable local
config yaml. This is intended for users who want to customize the
supported config settings without editing the installed package.

## Usage

``` r
init_local_config(config_path = get_config_path(), overwrite = FALSE)
```

## Arguments

- config_path:

  path to the local config yaml to create. Defaults to
  [`get_config_path()`](https://nichd-bspc.github.io/carnation/reference/get_config_path.md).

- overwrite:

  logical indicating whether to overwrite an existing file.

## Value

Path to the local config yaml, invisibly.

## Examples

``` r
cfg_out <- tempfile(fileext = ".yaml")
init_local_config(cfg_out)
#> Error in init_local_config(cfg_out): could not find function "init_local_config"
```
