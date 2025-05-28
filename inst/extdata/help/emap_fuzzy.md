#### Fuzzy enrichment map
-------------------------

Enrichment map view of the `Fuzzy clustering` meta-analysis.
- Unlike traditional clustering where terms belong to only one cluster, fuzzy clustering allows terms to participate in multiple clusters.
- This approach is particularly useful for understanding complex biological responses where pathways may participate in multiple processes.

**What it shows:**
- Network visualization of fuzzy-clustered functional terms
- Representative terms and cluster members
- Overlapping cluster membership through network structure

**When to use it:**
- To visualize the fuzzy clusters identified in the Fuzzy clustering table
- To explore how terms belong to multiple biological themes
- To identify relationships between overlapping functional clusters

**How to interpret:**
- Nodes represent functional terms or pathways
- Colors represent different fuzzy clusters
- `Representative` terms in a cluster are outlined in black
- Node size represents the number of DE genes associated with the term
- Edge thickness represents the number of genes shared between terms
- Terms may appear in multiple clusters, reflecting biological complexity

**Plot options:**
- `# of terms`: Control how many terms appear in the plot (default: 30)