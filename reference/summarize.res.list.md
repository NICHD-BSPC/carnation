# Combine everything in the results list into a single table

Combine everything in the results list into a single table

## Usage

``` r
summarize.res.list(
  res.list,
  dds.list,
  dds_mapping,
  alpha,
  lfc.thresh,
  labels = NULL
)
```

## Arguments

- res.list:

  Named list of lists, where each sublist contains the following names:
  c('res', 'dds', 'label'). "res" is a DESeqResults object, "dds" is
  either the indexing label for the dds.list object or the DESeq object,
  and "label" is a nicer-looking label to use. NOTE: backwards
  compatibility with older versions of lcdb-wf depends on no dds.list
  object being passed.

- dds.list:

  List of DESeqDataSet objects whose names are expected to match 'dds'
  slots in the 'res.list' object

- dds_mapping:

  List mapping names of dds.list to res.list elements

- alpha:

  false-discovery rate threshold

- lfc.thresh:

  log2FoldChange threshold

- labels:

  list of descriptions for res.list elements

  NOTE: this is edited to match the structure used in the shiny app.
  Specifically res.list and dds.list both are expected to have the same
  number and names of elements.

## Value

Dataframe
