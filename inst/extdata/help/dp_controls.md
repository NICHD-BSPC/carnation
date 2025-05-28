### Pattern Analysis Controls
-----------------------------

#### Analysis Selection

- `Analysis`: select pattern analysis to view from this menu
  - Each analysis may represent different clustering approaches or parameters

#### Plot Options

- `cluster column`: view results from cuts at different heights of the clustering dendrogram
  - Cluster column names start with `cutoff` and are followed by a number denoting the height of the cut
  - Smaller numbers correspond to lower cuts and finer/smaller clusters and vice-versa
  - The final level of the cut is displayed by default (`cluster`)

- `x axis variable`: choose metadata variable to show on x-axis
  - This denotes the temporal element of the analysis (e.g., time points,
   doses), typically representing the experimental condition that changes across samples

- `group by`: choose metadata variable to connect (color) groups of samples
  - Selecting `none` uses a single color for all pattern lines
  - Useful for distinguishing between experimental groups or conditions

- `label`: which genes should be labeled in the pattern plot
  - By default, this is `gene_scratchpad`, so genes in the scratchpad will be highlighted
  - Alternatively, genes in selected `upset_intersections` can be labeled
  - If any genes are labeled, a table is shown below the plot indicating cluster membership

#### Cluster Settings

- `Clusters to show`: select which clusters to display on the plot
  - All clusters are shown by default
  - Deselect clusters to focus on specific expression trends
  - Useful when there are many clusters and you want to focus on a few of interest
  - `Select all` and `Select none` buttons can be used to quickly reset/empty
    the list of clusters shown on the plot.
- `min cluster size`: minimum cluster size to show on the plot (default=10)
  - Increase this to filter out small clusters and make the plot less busy

#### X-axis Settings

- `variable levels`: drag-and-drop control to change the order of x-axis labels or to drop specific levels
  - Click and drag levels from `current` to `unused` to drop them from x-axis view and vice-versa
  - Drag levels up or down to change the order of levels
  - Use `Select all` or `Select none` to quickly add/remove all levels
  
- `rotate x labels`: rotation angle of x-axis labels
  - Change this to help declutter the display

#### Advanced Options

- `trendline`: type of trendline to connect groups
  - Can be `line` (default: straight lines), `smooth` (smoothed trend), or `none` (no connecting lines)

- `show boxes?`: toggle display of median-whisker boxplots for each group
  - Set to `TRUE` (default) to show boxplots

- `show points?`: toggle display of individual gene expression points
  - Set to `TRUE` to show points corresponding to each gene (default=`FALSE`)

- `show lines?`: toggle display of individual gene expression lines
  - Set to `TRUE` (default) to show individual gene trends


