#### Summary overview
---------------------

Horizontal bar chart depicting the most significantly enriched terms in
the analysis.
  - *This plot can be used to easily highlight the most enriched terms*.
  - The x-axis depicts `-log10(p-value)` - higher values correspond to greater
    statistical significance.
    - The `p-value column` can be either `pvalue` or `p.adjust`.
  - The color of the circles can be used to depict either `z_score` (default) or `aggr_score`,
    - `z_score`: measure of the direction of change. If `u` is the number of upregulated genes
      and `d` the number of downregulated genes, then `z = (u - d)/sqrt(u + d)`.
    - `aggr_score` is a measure of the size of the effect.

- `More options`
  - Control the `# of terms` depicted in the plot (default=20).
