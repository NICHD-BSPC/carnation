#### Enrichment map
-------------------

Interactive network visualization of the enrichment analysis.
- This plot is particularly good for identifying similarities between enriched terms, since terms with shared genes will cluster together in the network visualization.
- This helps reveal the underlying structure of your functional enrichment results.

**What it shows:**
- Network of related functional terms connected by shared genes
- Significance of enrichment through node color
- Gene set size through node size
- Strength of term relationships through edge thickness

**When to use it:**
- To identify clusters of related biological processes
- To visualize relationships between enriched terms
- To discover functional modules in your enrichment results

**How to interpret:**
- Nodes represent functional terms or pathways
- Edges represent shared genes between terms
- Node color corresponds to p-values (darker = more significant)
- Node size represents the number of DE genes associated with the term
- Edge thickness represents the number of genes shared between terms
- Terms with many shared genes cluster together in the network

**Plot options:**
- `# of terms`: Control how many terms appear in the plot (default: 30)
