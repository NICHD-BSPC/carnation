# Combine everything in the results list into a single table

Combine everything in the results list into a single table

## Usage

``` r
summarize_res_list(
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

## Value

Dataframe

## Examples

``` r
n_genes <- 100

#  make mock dds list
dds_list <- list(main=DESeq2::makeExampleDESeqDataSet(n=n_genes))

# make mock results df
res1 <- data.frame(
          baseMean = runif(n_genes, 10, 1000),
          log2FoldChange = rnorm(n_genes, 0, 2),
          lfcSE = runif(n_genes, 0.1, 0.5),
          stat = rnorm(n_genes, 0, 3),
          pvalue = runif(n_genes, 0, 1),
          padj = runif(n_genes, 0, 1),
          symbol = paste0("GENE", 1:n_genes),
          row.names = paste0("gene", 1:n_genes)
        )

res2 <- data.frame(
          baseMean = runif(n_genes, 10, 1000),
          log2FoldChange = rnorm(n_genes, 0, 2),
          lfcSE = runif(n_genes, 0.1, 0.5),
          stat = rnorm(n_genes, 0, 3),
          pvalue = runif(n_genes, 0, 1),
          padj = runif(n_genes, 0, 1),
          symbol = paste0("GENE", 1:n_genes),
          row.names = paste0("gene", 1:n_genes)
        )

# make list of results
res_list <- list(
              comp1=res1,
              comp2=res2
            )

# make dds mapping
dds_mapping <- list(comp1='main', comp2='main')

# get summary
df <- summarize_res_list(res_list, dds_list, dds_mapping, alpha=0.1, lfc.thresh=0)
```
