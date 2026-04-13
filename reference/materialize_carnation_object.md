# Materialize expensive carnation object components

This function materializes expensive derived pieces for a validated
carnation object, including DESeqDataSet creation from raw count
matrices, variance-stabilized counts, and GeneTonic conversions.

## Usage

``` r
materialize_carnation_object(obj, config = NULL, cores = NULL)
```

## Arguments

- obj:

  A validated object returned by
  [`validate_carnation_object()`](https://nichd-bspc.github.io/carnation/reference/validate_carnation_object.md)
  or `validate_loaded_carnation_object()`.

- config:

  Optional config list. If NULL, will use
  [`get_config()`](https://nichd-bspc.github.io/carnation/reference/get_config.md).

- cores:

  Optional number of worker processes. If NULL, uses
  `config$server$cores`.

## Value

The input object with materialized `dds_list`, `rld_list`, and optional
`genetonic` slots.

## Examples

``` r
if (FALSE) { # interactive()
# Minimal example with DE results and counts
library(DESeq2)

# Create example data
dds <- makeExampleDESeqDataSet()
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "A", "B"))
rld <- varianceStabilizingTransformation(dds, blind = TRUE)

# Validate object inputs
obj <- validate_carnation_object(
  res_list = list(
    comp1 = list(
      res = as.data.frame(res),
      dds = "main",
      label = "A vs B"
    )
  ),
  dds_list = list(main = dds),
  rld_list = list(main = rld)
)

materialized <- materialize_carnation_object(obj, cores = 1)
}
```
