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
- Mousing over a point displays a gene's name and exact values
- Switch between comparing log2FoldChange values or adjusted p-values
- Toggle between interactive and static plot modes
- Display a detailed comparison table below the plot
- Adjust plot aesthetics through the settings menu
- Gene names entered in the `Gene scratchpad` are labeled on the static plot and circled on the interactive plot
  - Note that having many labels may overwhelm the static plot, so in that situation we recommend using the interactive plot

**Download Options:**
- Click the `Download` button to save a publication-ready PDF version of the static plot
- The downloaded plot preserves all current settings including:
  - Comparison selections
  - Plot type (LFC or p-value)
  - FDR and LFC thresholds
  - Gene highlights
  - Color scheme
  - Axis limits

**Note:** Significance is determined by the global FDR and LFC threshold settings.
