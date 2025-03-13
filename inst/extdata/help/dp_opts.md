#### Plot options
-----------------

- `cluster column`
  - Use this menu to see results from *cuts* at different
    heights of the clustering dendrogram.
  - Cluster column names start with `cutoff` and are followed by a number denoting the height of the
    cut.
  - *Smaller numbers correspond to lower cuts and finer/smaller clusters and vice-versa.*
  - The final level of the cut is displayed by default (`cluster`).
- `x-axis variable`: metadata variable to show on x-axis
  - This denotes the temporal element of the analysis.
- `group by`: metadata variable to connect groups of samples
  - Different linetypes will highlight differing trends (default=`none`).
- `label`: which genes should be labeled?
  - By default, this is `gene_scratchpad`, so genes in the scratchpad will be
    highlighted in the pattern plot.
  - Alternatively, genes in selected `upset_intersections` can be labeled. If this
    option is selected, a dropdown menu to select UpSet intersections is shown below.
  - If any genes are labeled in the plot, a table is shown below the plot to indicate
    the cluster membership of labeled genes.
