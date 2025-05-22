#### Fuzzy Clustering
---------------------

Interactive table showing functional terms grouped using a fuzzy clustering approach.
This is a reimplementation of the algorithm underlying DAVID (https://david.ncifcrf.gov/), allowing terms to appear in multiple clusters.
This approach is particularly useful for understanding complex biological responses where pathways may participate in multiple processes.

**What it shows:**
- Functional terms organized into overlapping clusters
- Representative terms for each cluster
- Directional statistics for gene expression within clusters

**When to use it:**
- To identify related functional categories
- To allow terms to belong to multiple biological themes
- To explore directional changes within functional groups

**How to interpret:**
- `fuzzycluster`: Cluster identifier number
- `cluster_status`: Indicates whether a term is the `Representative` (most significant) or a `Member` of a cluster
- `genes`: Genes associated with the term in this cluster
- `z_score`: Measure of expression direction calculated as (u-d)/√(u+d), where:
  - u = number of upregulated genes
  - d = number of downregulated genes
  - Positive values indicate predominantly upregulated genes
  - Negative values indicate predominantly downregulated genes
- `aggr_score`: Mean log fold change of genes associated with the term

**Interactive features:**
- Filter by cluster number
- Sort by any column
- Search for specific genes or terms
- Adjust the number of terms included using the `# of terms` control
- View the corresponding visualization in the `emap_fuzzy` plot


