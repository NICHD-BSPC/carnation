#### Upset Plot
---------------

Visualization of intersections between DE gene sets of interest

- The top panel shows a bar chart of the intersection size.
- The bottom panel shows the gene sets being compared on separate rows
  - Bars to the left of set names denote the size of the gene sets
  - Dots to the right of set names indicate the intersection.
     - A single dot on any row depicts genes present only in that gene set
     - Dots connected with lines depict genes shared by those gene sets
- `Intersections` selected in the sidebar are highlighted on the plot.
- If genes are selected in the `Gene scratchpad`, then
  a table is shown below the UpSet plot to show which
  intersections they are part of.
  - A subset of columns from the `Table` view are shown here.
- `FDR threshold` and `log2FC threshold` are used
  for filtering genes before comparison. These controls
  are shared with `Table`, `MA plot` & `Heatmap` tabs of `DE analysis`, so changes made here are also reflected
  there.
