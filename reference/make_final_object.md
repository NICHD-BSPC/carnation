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
