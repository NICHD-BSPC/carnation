# Validate Pattern Analysis Object Schema

Validate the schema for a single `degpatterns` analysis element used by
the pattern analysis module.

## Usage

``` r
is_valid_pattern_obj(pattern_obj, require_symbol = FALSE)
```

## Arguments

- pattern_obj:

  A single pattern analysis element. Must be either a `data.frame` or a
  list containing a `normalized` `data.frame`.

- require_symbol:

  Logical, if `TRUE` require a `symbol` column in the analysis table.

## Value

Returns `TRUE` when validation succeeds, otherwise returns `FALSE` after
emitting a message describing the issue.

## Examples

``` r
data(degpatterns_dex, package = "carnation")

is_valid_pattern_obj(degpatterns_dex)
#> Error in is_valid_pattern_obj(degpatterns_dex): could not find function "is_valid_pattern_obj"
```
