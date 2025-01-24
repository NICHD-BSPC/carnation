#### Scatter plot comparisons
-----------------------------

- Select `Comparison 1` & `Comparison 2` to compare using either log2FoldChange (`LFC`) or
  adjusted p-values (`P-adj`). Switch the plotted values with the `Values to use` menu.
- An interactive or static plot can be selected using the `Interactive?` menu.
- If `Show table?` is `yes`, then a table is shown below the plot with `LFC` & `P-adj`
  values from the selected comparisons.
- Various elements of the plot can be further adjusted using the `Plot settings` menus
  - *Axes limits*

    Adjust x & y-axis limits or autoscale to fit the data.
  - *Point aesthetics*

    Adjust properties of plotted points (size & opacity) and the overall `Color palette`.
    Also, we can choose to show all points or skip those that fall outside axis limits.
    If `Show all points?` is *yes* then out-of-bound points are shown directly on the
    axes borders.
  - *Grid lines*

    Choose grid lines to show. There are four options that can be individually turned
    on or off (by default, all are plotted):
    - x = 0
    - y = 0
    - diagonal (x = y)
    - x-y grid
