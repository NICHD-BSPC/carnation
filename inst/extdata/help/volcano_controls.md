### Volcano Plot Controls
-------------------------

#### Main controls

- `Comparison`: select which differential expression comparison to display.
  - The same FDR and LFC thresholds are applied to each comparison
- `Interactive?`: choose between interactive and static plot modes.
  - `yes`: creates an interactive plot with hover functionality
  - `no`: creates a static plot suitable for publication

Changes in the main controls results in the plot being immediately redrawn.

#### Plot settings

- `Color by`: choose whether points are colored by significance or by average expression (base mean)
- `Opacity`: control the transparency of points
  - Lower values help visualize overlapping points
  - Higher values make individual points more distinct
- `Aspect ratio`: choose 'wide' or 'narrow' to switch between a wide or narrow view of the plot.
- `axis limits`
  - `x min`, `x max`: set the range of the x-axis (log2 fold change)
  - `y min`, `y max`: set the range of the y-axis (-log10 adjusted p-value)
  - Genes beyond these limits appear as triangles at the edges of the plot
  - Useful for zooming into the crowded base of the plot or isolating the most significant genes at the top
  - `Autoscale`: sets both axes to include all data points with a small margin
    - Use this to reset the view after zooming into a specific range

Note: You need to click `Refresh` for changes in "Plot Settings" to have an effect.
