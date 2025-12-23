# Convert enrichResult to GeneTonic object

This function takes an enrichResult object and DE analysis results and
creates a GeneTonic object.

## Usage

``` r
enrich_to_genetonic(enrich, res)
```

## Arguments

- enrich:

  enrichResult object

- res:

  data frame with DE analysis results

## Value

GeneTonic object

## Examples

``` r
# get enrich & res objects
data(res_dex, package="carnation")
data(eres_dex, package="carnation")

# convert to GeneTonic object
gt <- enrich_to_genetonic(eres_dex, res_dex)
#> Found 2186 gene sets in `enrichResult` object, of which 2186 are significant.
#> Converting for usage in GeneTonic...

```
