#### Fuzzy clustering
---------------------

Here we cluster top functional terms, using a fuzzy approach which means that
terms can appear in multiple clusters. This is a reimplementation of the algorithm
underlying DAVID (https://david.ncifcrf.gov/).
- The table shares several columns with the `Enrichment` table, but has the following
  additional columns:
  - `fuzzycluster`: cluster number
  - `cluster_status`: this can be either `Representative` or `Member`. The most significant
    term in the cluster is termed `Representative` while all others are `Members`.
  - `genes`: genes present in the cluster
  - `z_score`: measure of the direction of change. If `u` is the number of upregulated genes
    and `d` the number of downregulated genes, then `z = (u - d)/sqrt(u + d)`.
  - `aggr_score`: mean LFC of genes associated with the term
- The number of functional terms included in the fuzzy clustering can be adjusted using
  the `# of terms` control. As above, if this number is greater than the number of terms
  in the enrichment table, then the latter is used instead.
- The data shown in this table is visualized in the `emap_fuzzy` plot.

