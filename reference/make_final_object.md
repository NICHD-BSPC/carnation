# Make final object for internal use by the app

This function takes an uploaded object and sanitizes it to make sure it
is suitable for internal use along with other additions:

- adds a 'dds_mapping' element that maps dds_list keys to res_list
  objects.

- if there are multiple dds_list objects, it adds a 'all_dds' element
  combining all samples.

## Usage

``` r
make_final_object(obj)
```

## Arguments

- obj:

  list object containing lists of DE analysis results, functional
  enrichment objects, pattern analysis objects & raw and normalized
  counts objects.

## Value

final carnation object with additional pre-processing

## Examples

``` r
library(DESeq2)

# make example DESeq dataset
dds <- makeExampleDESeqDataSet()

# run DE analysis
dds <- DESeq(dds)
#> estimating size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing

# extract comparison of interest
res <- results(dds, contrast = c("condition", "A", "B"))

# perform VST normalization
rld <- varianceStabilizingTransformation(dds, blind = TRUE)

# build minimal object
obj <- list(
           res_list = list(
                          comp = list(
                              res = res,
                              dds = "main",
                              label = "A vs B"
                          )
                      ),
           dds_list = list(main = dds),
           rld_list = list(main = rld)
       )

# final object
final_obj <- make_final_object(obj)
```
