#### MA Plot
------------

Interactive visualization showing the relationship between average expression and log fold change for all genes.

**What it shows:**
- Average expression level (`baseMean`) on the x-axis
- Log2 fold change (`log2FoldChange`) on the y-axis
- Statistical significance through color coding

**When to use it:**
- To identify differentially expressed genes
- To assess the distribution of expression changes
- To detect potential biases in your differential expression analysis

**How to interpret:**
- Red points represent statistically significant differentially expressed genes (based on FDR and LFC thresholds)
- Grey points represent non-significant genes
- Triangles at the top/bottom indicate genes with fold changes beyond the y-axis limits
- Genes from your "Gene scratchpad" are highlighted with dark circles

**Interactive features:**
- Hover over points to see gene names and exact values
- Adjust FDR and LFC thresholds using the sidebar controls
- Modify y-axis limits to focus on specific fold change ranges
- Download a publication-ready PDF version

**Note:** Changes to the FDR and LFC thresholds affect all DE analysis visualizations.

