#### Distill enrichment
-----------------------

Here the top functional terms detected in the analysis are clustered based on similarities in
associated genes. This provides a way to look for similar/dissimilar sets of functional terms enriched
among genes of interest.
- The table consists of the following columns:
  - `cluster`: cluster number
  - `#_terms`: number of terms grouped in the cluster
  - `#_genes`: number of genes that are associated with the grouped terms
  - `most_significant_term`: this is the term with the lowest p-value in the cluster
  - `strongest_term`: this is the most highly connected member of the cluster
  - `genes`: genes present in the cluster
  - `term_description_list`: comma-separated list of descriptions of terms in the cluster
- The number of functional terms in this meta-analysis can be adjusted using the `# of terms` control.
  If this number is greater than the number of rows in the enrichment table, then this
  automatically defaults to the latter.
- The data shown in this table is visualized in the `emap_distill` plot.

