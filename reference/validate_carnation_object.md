# Validate a carnation object

This function takes various input data types (DE results, counts,
enrichment, pattern analysis) and validates them according to
carnation's requirements, returning a normalized intermediate object.
Expensive derived-object creation steps such as variance-stabilized
counts and GeneTonic conversion are handled separately by
[`materialize_carnation_object()`](https://nichd-bspc.github.io/carnation/reference/materialize_carnation_object.md).

## Usage

``` r
validate_carnation_object(
  res_list,
  dds_list,
  rld_list = NULL,
  labels = NULL,
  enrich_list = NULL,
  degpatterns = NULL,
  metadata = NULL,
  dds_mapping = NULL,
  config = NULL
)
```

## Arguments

- res_list:

  Named list of DE results. Each element should be either:

  - A data frame with DE results containing gene, symbol, pvalue, padj,
    log2FoldChange, and baseMean columns (or tool-specific alternatives)

  - A list with slots: `res` (data frame), `dds` (name reference to
    dds_list element), `label` (comparison label)

- dds_list:

  Named list of count data. Each element should be either:

  - A `DESeqDataSet` object

  - A data frame or matrix of raw counts (first column=gene IDs,
    remaining=samples)

- rld_list:

  Optional named list of variance-stabilized count objects. If NULL,
  these can be generated later via
  [`materialize_carnation_object()`](https://nichd-bspc.github.io/carnation/reference/materialize_carnation_object.md).

- labels:

  Optional named list of comparison labels. If NULL and res_list
  contains nested structure with `label` slots, labels will be
  extracted.

- enrich_list:

  Optional named list of functional enrichment results. Should be
  structured as: `enrich_list[[func_id]][[effect]][[pathway]]`. Each
  enrichment result must be a data frame in clusterProfiler format:

  - Over-representation: ID, Description, GeneRatio, BgRatio, pvalue,
    p.adjust, qvalue, geneID, Count

  - GSEA: ID, Description, core_enrichment, setSize, pvalue, p.adjust,
    qvalue, NES

- degpatterns:

  Optional named list of pattern analysis results. Each element should
  be either a data frame or a list with `$normalized` slot containing a
  data frame with columns: genes, value, and either cluster or columns
  starting with "cutoff".

- metadata:

  Optional data frame with sample metadata. Required if `dds_list`
  contains count matrices instead of DESeqDataSet objects. First column
  should be sample names matching column names in count matrices.

- dds_mapping:

  Optional named list mapping `res_list` elements to `dds_list` objects.
  Required if `res_list` is a list of data frames.

- config:

  Optional config list. If NULL, will use
  [`get_config()`](https://nichd-bspc.github.io/carnation/reference/get_config.md),
  including any supported local config overrides.

## Value

A validated list with canonical slots `res_list`, `dds_list`, optional
`rld_list`, `labels`, `dds_mapping`, `enrich_list`, `degpatterns`, and
`metadata` when supplied.

A list containing normalized inputs with elements `res_list`,
`dds_list`, optional `rld_list`, `labels`, `dds_mapping`, and optional
`enrich_list`, `degpatterns`, and `metadata`.

## Details

This function performs comprehensive validation of all input data:

- DE results: Checks for required columns (with support for DESeq2,
  edgeR, limma), ensures gene and symbol columns exist

- Counts: Validates structure, checks sample name matching with metadata

- Enrichment: Validates clusterProfiler format (OR or GSEA)

- Pattern analysis: Checks for required columns (genes, value, cluster)

If validation fails, the function will stop with an informative error
message.

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
final_obj <- make_final_object(materialized)

# Save for use with carnation
saveRDS(final_obj, "my_analysis.rds")

# Alternative: start from count matrix and metadata
counts <- as.data.frame(counts(dds))
counts$gene <- rownames(counts)
counts <- counts[, c(ncol(counts), 1:(ncol(counts)-1))]
metadata <- as.data.frame(colData(dds))
metadata$sample <- rownames(metadata)
metadata <- metadata[, c(ncol(metadata), 1:(ncol(metadata)-1))]

obj <- validate_carnation_object(
  res_list = list(comp1 = as.data.frame(res)),
  dds_list = list(main = counts),
  metadata = metadata,
  dds_mapping = list(comp1 = "main")
)
}
```
