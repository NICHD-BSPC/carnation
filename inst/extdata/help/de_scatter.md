#### Scatter plot
-----------------

Interactive plot comparing two comparisons using either adjusted p-values (`padj`) or log2 fold changes (`LFC`).

- Mousing over a point displays a gene's name.
- Genes are defined as significant based on the FDR and LFC thresholds, which can be adjusted using the `Global settings`
  controls in the side panel.
- Genes may belong to one of 5 significance categories:
  1. Significant in both comparisons + same direction of change in both comparisons.
  2. Significant in both comparisons + opposite direction of change.
  3. Significant only in comparison 1 (x-axis) having either positive or negative LFC.
  4. Significant only in comparison 2 (y-axis) having either positive or negative LFC.
  5. Significant in neither comparison 1 or 2.
- The `Settings` button in the side panel allows you to configure the data that is plotted and to adjust a number
  of visual settings including axis limits and point color palettes.
- Gene names entered in the `Gene scratchpad` are labeled on the static plot and circled on the interactive plot.
  - Note that, having many labels may overwhelm the static plot, so in that situation we recommend using the interactive plot.
- A PDF image of the plot can be downloaded using the `Download` button.
