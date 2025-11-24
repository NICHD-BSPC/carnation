# Get read counts for gene

This is a simple function to obtain read counts for a specified gene,
based on the DESeq2::plotCounts function.

## Usage

``` r
get_gene_counts(dds, gene, intgroup = "condition", norm_method = "libsize")
```

## Arguments

- dds:

  DESeqDataSet object

- gene:

  gene name vector

- intgroup:

  metadata variable to attach to counts

- norm_method:

  normalization method, can be 'libsize' (default) or 'vst'

## Value

data.frame with gene counts
