#### Cnetplot
-------------

Network visualization showing the associations between genes and multiple functional terms or pathways.
- This plot is particularly useful for identifying the specific genes that drive enrichment results and for discovering genes that participate in multiple biological processes, potentially revealing key regulators in your system.

**What it shows:**
- Connections between genes and the functional categories they belong to
- Relative importance of terms through node size
- Genes that link multiple functional categories

**When to use it:**
- To identify which genes act as bridges between different biological processes
- To visualize the overlap between related functional categories
- To explore the specific genes driving your enrichment results

**How to interpret:**
- Larger nodes represent terms with more associated genes
- Genes connecting to multiple terms may be key regulators or multifunctional proteins
- Densely connected regions indicate related biological processes

**Plot options:**
- `# of terms`: Control how many terms appear in the plot (default: 5)
- `node labels`: Choose to display `all` labels, only `category` names, only `gene` names, or no labels
- `max name length`: Limit the character length of term names for better display (default: 50)
- `color edge by terms`: Color edges based on their connected terms
- `circular layout`: Arrange the network in a circular pattern instead of force-directed layout