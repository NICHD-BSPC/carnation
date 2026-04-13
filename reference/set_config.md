# Set config

This function updates a limited subset of the package config YAML. Only
stable user-facing settings are writable; style settings and other
internal options are intentionally left untouched.

## Usage

``` r
set_config(
  config_path = get_config_path(),
  de_analysis = NULL,
  fdr_threshold = NULL,
  log2fc_threshold = NULL,
  max_upload_size = NULL,
  cores = NULL
)
```

## Arguments

- config_path:

  character path to the config YAML file to update. Defaults to the
  local config returned by
  [`get_config_path()`](https://nichd-bspc.github.io/carnation/reference/get_config_path.md).
  If the file does not exist yet, it is initialized from the bundled
  package config.

- de_analysis:

  optional list with DE analysis config updates. Currently only
  `de_analysis$column_names` is supported, and the provided aliases are
  merged into the existing column-name mappings.

- fdr_threshold:

  optional numeric FDR threshold between 0 and 1.

- log2fc_threshold:

  optional numeric log2 fold-change threshold greater than or equal to
  0.

- max_upload_size:

  optional positive numeric upload limit in MB.

- cores:

  optional positive integer number of cores to use.

## Value

Updated config list, invisibly.

## Examples

``` r
cfg_out <- tempfile(fileext = ".yaml")

set_config(
  config_path = cfg_out,
  de_analysis = list(
    column_names = list(
      padj = "qvalue",
      log2FoldChange = c("logFC", "avg_log2FC")
    )
  ),
  fdr_threshold = 0.05,
  log2fc_threshold = 1,
  max_upload_size = 50,
  cores = 2
)
#> Error in set_config(config_path = cfg_out, de_analysis = list(column_names = list(padj = "qvalue",     log2FoldChange = c("logFC", "avg_log2FC"))), fdr_threshold = 0.05,     log2fc_threshold = 1, max_upload_size = 50, cores = 2): could not find function "set_config"
```
