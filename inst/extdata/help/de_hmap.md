#### Heatmap
------------

Hierarchical clustering visualization of expression patterns for top differentially expressed genes.

**What it shows:**
- Normalized expression values across samples and conditions
- Hierarchical clustering of genes with similar expression patterns
- Sample relationships based on expression profiles

**When to use it:**
- To visualize expression patterns across experimental conditions
- To identify co-regulated gene clusters
- To examine sample similarity based on DE gene expression

**How to interpret:**
- Each row represents a gene
- Each column represents a sample
- Color intensity indicates expression level (typically red for high, blue for low)
- The dendrogram shows hierarchical relationships between genes
- Genes found in the Gene scratchpad are labeled in red on the side of the heatmap

**Plot options:**
- Select different sample groups to visualize
- Choose clustering methods (row, column, both, or none)
- Scale data by row, column, or none
- Limit the number of genes shown (top genes based on significance)

**Note:** Changes to the FDR and LFC thresholds affect all DE analysis visualizations.
