#### Upset Intersection Table
-----------------------------

- To view a specific intersection, select it from the `Intersections`
  drop-down menu.
- Intersections shown in the Upset plot are labeled based on the
  order in which they appear. For example, the first set is
  called `set01`, the second set is `set02` and so on.
- The following columns are always shown:
  - `gene`: gene ID
  - `symbol`: gene symbol
  - `set`: intersection set numbered based on UpSet plot
  - `comparisons`: comparisons corresponding to selected intersection
  - `degree`: number of comparisons for which the gene is differentially expressed.
- The remaining columns show `log2FoldChange` and `padj`
  columns from the gene sets being compared (selected in
  `Comparisons`).
- By default, the table is sorted by the `set` column, and
  then by the first `log2FoldChange` column of the comparisons
  that are part of the set.
- `DE Filters` shown in the sidebar are used
  for filtering genes before comparison. These controls
  are shared with `Table`, `MA plot` & `Heatmap` tabs of `DE analysis`.
- `Selection options`
  - `Add to scratchpad` & `Remove from scratchpad` buttons will add/remove selected genes
    from the scratchpad.
  - `Reset selection` clears the selection.

