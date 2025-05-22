#### Alluvial
-------------

Interactive Sankey diagram showing functional terms and associated genes.
- Similar to the cnetplot, this visualization shows the association of genes with multiple terms, but in a flow-based layout.
- This can provide a clearer view of how genes are distributed across functional categories, especially when genes participate in multiple biological  
  processes.

**What it shows:**
- Flow of connections between genes and functional terms
- Distribution of genes across multiple functional categories
- Relative importance of terms through flow width

**When to use it:**
- To visualize how genes are distributed across functional categories
- To identify genes that participate in multiple biological processes
- To create an alternative visualization to the cnetplot

**How to interpret:**
- Each colored flow represents a functional term or pathway
- The width of flows indicates the number of associated genes
- Genes connecting to multiple terms appear as flow splits
- Incoming flow count indicates the number of terms associated with a gene
- Outgoing flow count indicates the number of genes associated with a functional term

**Plot options:**
- `# of terms`: Control how many terms appear in the plot (default: 4)
