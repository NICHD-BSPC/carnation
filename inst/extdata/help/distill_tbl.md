#### Distill Enrichment
-----------------------

Interactive table showing clustered functional terms based on gene overlap similarity.
This analysis is helpful to identify the major biological themes in your enrichment results by grouping similar terms together.

**What it shows:**
- Functional terms grouped into clusters based on shared genes
- Summary statistics for each cluster
- Representative terms for each functional group

**When to use it:**
- To reduce redundancy in enrichment results
- To identify major functional themes in your data
- To simplify interpretation of complex enrichment results

**How to interpret:**
- `cluster`: Unique identifier for each functional term group
- `#_terms`: Number of related terms grouped in this cluster
  - If this exceeds the number of rows in the enrichment table,
    it automatically defaults to the latter.
- `#_genes`: Total number of genes associated with the grouped terms
- `most_significant_term`: Term with the lowest p-value in the cluster
- `strongest_term`: Most highly connected member of the cluster
- `genes`: Genes present in the cluster
- `term_description_list`: Comma-separated descriptions of all terms in the cluster

**Interactive features:**
- Sort clusters by any column
- Search for specific genes or terms
- Adjust the number of terms included using the `# of terms` control
- View the corresponding visualization in the `emap_distill` plot

