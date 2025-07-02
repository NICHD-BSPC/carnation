#### Distill enrichment map
---------------------------

Enrichment map view of the `Distill enrichment` meta-analysis.
- This visualization complements the Distill enrichment table by showing the network structure of the identified clusters.
- The clustering helps reduce redundancy in enrichment results and highlights the major biological themes in your data.

**What it shows:**
- Network visualization of clustered functional terms
- Relationships between terms within and across clusters
- Significance and gene content through node properties

**When to use it:**
- To visualize the clusters identified in the Distill enrichment table
- To explore relationships between functional term clusters
- To identify higher-level biological themes in your results

**How to interpret:**
- Nodes represent functional terms or pathways
- Colors represent different clusters identified by the distill algorithm
- Node size represents the number of DE genes associated with the term
- Edge thickness represents the number of genes shared between terms
- Terms within the same cluster have the same color

**Plot options:**
- `# of terms`: Control how many terms appear in the plot (default: 30)
