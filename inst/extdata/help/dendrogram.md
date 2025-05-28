#### Dendrogram
---------------

Tree visualization showing hierarchical clustering of functional terms.
- This plot provides an alternative to enrichment maps for grouping terms based on shared genes.
- The hierarchical structure can reveal relationships between terms at different levels of similarity, helping to identify both broad functional categories and specific sub-processes.

**What it shows:**
- Hierarchical relationships between functional terms based on shared genes
- Clustering of related biological processes
- Expression direction through leaf coloring

**When to use it:**
- To visualize hierarchical relationships between functional terms
- To identify groups of related biological processes
- As an alternative to enrichment maps for term clustering

**How to interpret:**
- Each leaf represents a functional term or pathway
- Branches represent clusters of related terms
- Color of leaves represents z-scores (direction of expression change)
- Size of leaves represents p-values (significance of enrichment)
- Color of branches represents different clusters
- Terms that cluster together share more genes

**Plot options:**
- `# of terms`: Control how many terms appear in the plot (default: 30)
- `max name length`: Limit the character length of term names for better display
