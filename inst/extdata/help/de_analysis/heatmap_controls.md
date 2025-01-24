#### Heatmap controls
---------------------

Main controls to adjust the heatmap

- `sample group`: sample subset
  to be used for plotting heatmap (default: `all_samples`).
- `cluster rows`: should the rows be clustered?
  - If `TRUE` (default), the rows are clustered and a dendrogram is
    shown to the right of the heatmap.
- `cluster columns`: should the columns be clustered?
  - If `TRUE`, the columns are clustered and a dendrogram is shown above
    the plot.
  - By default, this is set to `FALSE` and column ordering is set via
    `column settings`
- `scale`: should the data be scaled?
  - By default, the data are scaled by `row`.
