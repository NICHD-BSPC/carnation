### PCA Controls
-----------------

#### Sample Group Selection

- `Sample group`: sample groups that exist in the analysis.
  - `all_samples` denotes a group containing all samples.
  - `subsets`: on selecting this option, a new menu `choose subset`
    will be shown where sample groups that exist in the analysis
    will be shown.

#### Plot Options

- `color by`: color samples by metadata columns
  - multiple selections result in samples being colored with all combinations
    of selected columns
- `x-axis`, `y-axis`, `z-axis`: choose different principal components (PCs) to
  show on the x-, y- or z-axis.
  - this is useful to explore the data when PC1 & PC2 don't account for
    a large fraction of the sample variation
  - Selecting a third principal component with `z-axis` (optional) creates a 3D plot
  - PC1-PC6 are currently supported
- `# of genes`: number of most-variable genes to use for PCA calculation
  - default is 500 genes
  - increasing this number may capture more subtle patterns
  - decreasing this number focuses on the most variable genes

#### Advanced Options

- `subset samples`: select samples to show on the PCA plot using a
  `grouping variable` selected from metadata
  - samples can be added or removed using levels of the selected grouping variable
    to observe the effect on the PCA plot.
- `gene loadings`: turn on gene loadings overlay and choose `# of top genes`
  to plot.
  - gene loadings show which genes contribute most to each principal component
  - arrows point in the direction of increasing expression
  - longer arrows indicate stronger contribution to the principal components
  - the `# of top genes` control lets you adjust how many gene loadings to display
