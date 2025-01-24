#### MA plot
------------

Interactive plot of average expression (`baseMean`) vs LFC (`log2FoldChange`)

- Mousing over a point displays the gene name.
- The FDR, LFC thresholds and the y-axis limits can be adjusted using the plotting controls in the sidebar.
- Genes that are statistically significant (`DE`) based on the FDR and LFC thresholds are colored red. If not, they
  are colored grey (`not DE`).
- Genes appearing beyond the y-axis limits are plotted at the y-axis limits using triangles
  - Triangles point up or down based on gene LFCs above or below y-axis limits.
  - Triangles are colored red or grey based on statistical significance.
- Gene names entered into the `Gene scratchpad` are labeled on the plot using dark circles.
- A PDF version of the plot can be downloaded using the
  `Download` button.
  - In the downloaded plot, labeled genes will have their names
    printed directly on the plot and colored red or grey based on
    statistical significance.

