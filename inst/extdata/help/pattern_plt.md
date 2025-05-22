#### Pattern Analysis Plot
--------------------------

Visualization of gene expression patterns across conditions using hierarchical clustering.
This is particularly useful for experiments with a temporal component, helping identify major patterns of gene expression
changes over time or across treatment conditions.

**What it shows:**
- Groups of genes with similar expression trends
- Temporal or condition-dependent patterns
- Relative expression levels across your experimental conditions

**When to use it:**
- For time-series or dose-response experiments
- To identify co-regulated gene sets
- To discover predominant expression patterns in your data

**How to interpret:**
- Each line represents the mean expression pattern of a gene cluster
- The y-axis shows normalized expression values
- The x-axis represents your experimental conditions (e.g., time points, treatments)
- Shaded areas represent confidence intervals around the mean
- Genes found in the Gene scratchpad or selected Upset intersections are labeled on the plot

**Plot options:**
- Choose which metadata variable to use for the x-axis.
- Levels of the x-axis variable can be edited, e.g. to remove specific levels or to change the order in which they are shown.
- Filter out small clusters using the `min cluster size` control. Control which clusters are shown and/or their order
  using the `Clusters to show` menu.
- Adjust the visualization of trend lines and data points

**Interactive features:**
- Hover over lines to see gene names and exact values
- Download a publication-ready PDF version of the plot

**Note:** Changes to the FDR and LFC thresholds affect all DE analysis visualizations.
