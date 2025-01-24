#### Summary overview
---------------------

Horizontal bar chart depicting the most significantly enriched terms in
the analysis.
  - *This plot can be used to easily highlight the most enriched terms*.
  - The x-axis depicts `-log10(p-value)` - higher values correspond to greater
    statistical significance.
  - The color of the circles can be used to depict either `z_score` (default) or `aggr_score`,
    - `z_score`: measure of the direction of change. if `u` is the number of upregulated genes associated with the term and `d` is
      the number of downregulated genes, then `z=(u-d)/sqrt(u + d)`.
    - `aggr_score`: measure of the size of the effect, defined as the mean LFC of genes associated with the term.

**Plot options**
- `# of terms` - Number of terms shown in the plot (default=20).
- `p-value column` - Column used to define p-value depicted on the x-axis.
  Can be either `pvalue` or `p.adjust`.
- `color by`: color circles to show a measure of direction (`z_score`) or magnitude of effect (`aggr_score`).
- `max name length`: Maximum number of characters for term names shown on y-axis (default=50).
