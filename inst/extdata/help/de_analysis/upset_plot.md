#### Upset Plot
---------------

Visualization of intersections between DE gene sets of interest

- The view is split into two halves, the top half shows the
  Upset plot and the bottom shows a viewer of genes found
  in the intersections.
- `Upset plot`
  - The top panel shows a bar chart of the intersection size.
  - The bottom panel shows the gene sets being compared on separate rows
    - Bars to the left of set names denote the size of the gene sets
    - Dots to the right of set names indicate the intersection.
       - A single dot on any row depicts genes present only in that gene set
       - Dots connected with lines depict genes shared by those gene sets
- `Intersection viewer`
  - To view a specific intersection, select it from the `Intersections`
    drop-down menu.
     - Intersections shown in the Upset plot are labeled based on the
       order in which they appear. For example, the first set is
       called `set01`, the second set is `set02` and so on.
     - The label also includes the number of genes it contains.
  - Upon selection, the `comparisons` associated with the intersection
    are shown on the left, and the genes are shown in a
    table on the right.
- `Filters` shown in the sidebar are shared with `Table`, `MA plot` & `Heatmap` and used
  for filtering genes before comparison.
