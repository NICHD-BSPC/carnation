### MA Plot Controls
--------------------

#### Comparison Selection

- `Comparison`: select which differential expression comparison to display.
  - The same FDR and LFC thresholds are applied to each comparison

- `Interactive?`: choose between interactive and static plot modes.
  - `yes`: creates an interactive plot with hover functionality
  - `no`: creates a static plot suitable for publication

#### Plot Settings

- `Opacity`: control the transparency of points
  - Lower values help visualize overlapping points
  - Higher values make individual points more distinct
- `Aspect ratio`: choose 'wide' or 'narrow' to switch between a wide or narrow view of the plot.
- `y-axis limits`
  - `min`, `max`: adjust the range of the y-axis (log2 fold change).
    - Genes with fold changes beyond these limits appear as triangles at the top/bottom
    - Useful for focusing on specific fold change ranges of interest
  - `Autoscale`: clicking this button automatically sets the y-axis limits to include all data points with a small margin
    - Use this button when you want to reset the view after zooming in on a specific range

Note: You need to click `Refresh` for changes in "Plot Settings" to have an effect.
