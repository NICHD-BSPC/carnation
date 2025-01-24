#### Summary table
------------------

Table summarizing differential expression (DE) analysis results for
all comparisons in an analysis.

Column key:
- `comparison`: unique id for the comparison performed.
- `up`, `down`: number of genes up or downregulated.
- `total.genes`: number of annotated genes that are
  included in the analysis.
- `total.nonzero`: number of genes that do not have zero counts
  in all samples.
- `outliers`: number of genes that are flagged
  as count outliers by DESeq2. Typically, these are genes
  that have widely varying counts between replicates.
- `low.counts`: number of genes with low counts across all
  replicates that are filtered by DESeq2 in order to improve
  the p-value adjustment for multiple testing.
- `design`: linear model used for the analysis with DESeq2.

