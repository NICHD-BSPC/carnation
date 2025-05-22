#### Scatter Plot Comparisons
-----------------------------

Interactive visualization comparing gene expression changes between two different contrasts.

**What it shows:**
- Correlation between fold changes or p-values across two comparisons
- Genes with consistent or divergent responses between conditions
- Statistical significance in either or both comparisons

**When to use it:**
- To compare treatment responses across different conditions
- To identify genes with consistent or unique responses
- To discover potential interaction effects between experimental factors

**How to interpret:**
- Each point represents a gene
- Position indicates the gene's behavior in both comparisons
- Points in quadrants I and III show consistent direction of change
- Points in quadrants II and IV show opposite directions of change
- Color coding indicates significance categories:
  1. Significant in both comparisons with same direction of change
  2. Significant in both comparisons with opposite direction
  3. Significant only in comparison 1 (x-axis)
  4. Significant only in comparison 2 (y-axis)
  5. Not significant in either comparison

**Interactive features:**
- Switch between comparing log2FoldChange values or adjusted p-values
- Toggle between interactive and static plot modes
- Display a detailed comparison table below the plot
- Adjust plot aesthetics through the settings menu

**Note:** Significance is determined by the global FDR and LFC threshold settings.
