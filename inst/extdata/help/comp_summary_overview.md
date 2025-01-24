#### Summary overview
---------------------

Horizontal bar chart comparing two gene sets.
- **Note**: The circles with the solid black outline represents `Comparison 1`. The circles representing
  `Comparison 2` also have some transparency (`alpha = 0.7`) to help differentiate the two groups.
- The remaining details are shared with the single comparison, `summary_overview` plot.
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
