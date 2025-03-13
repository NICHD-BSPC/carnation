#### Upset options

- `FDR threshold`: false discovery rate (FDR) threshold (default: 0.1)
- `log2FC threshold`: log-2 fold-change (FC) threshold (default: 0)
  - Note that changes to these thresholds propagate to `Summary`, `Table`
    and `Heatmap` tabs.
- `Direction of change`:
  - `up & down`: all changed genes, both up- and down-regulated
  - `up`: upregulated genes
  - `down`: downregulated genes
- `Intersections`: intersections between comparisons shown in the UpSet plot.
  This menu is only shown for the `Table` tab.
  - Intersections are labeled based on the
    order in which they appear on the UpSet plot.
    - For example, the first set is called `set1`, the second set is `set2`, and so on.
    - Further, the name also contains within parentheses the comparisons sharing these genes and size of the intersection.
    - For instance, if `set1` consists of `25` genes shared by comparisons `comp_1` & `compn_2`, then
      the intersection is labeled `set 1 (comp_1 comp_2; n = 25)`.
  - By default, the first 10 intersections are selected on first load.
- `# of intersections`: The number of intersections shown on the UpSet plot (default: 40).
- `Min intersection size`: The minimum size of intersections shown on the plot (default: 5).
- `Text scale`: a numeric value used to scale text labels on the UpSet Plot.
  - Values > 1 will make the text larger and vice-versa.
