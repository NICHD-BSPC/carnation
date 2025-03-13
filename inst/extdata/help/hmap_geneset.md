#### Gene selection

- `Geneset`: the type of gene set being visualized. This can be `de_genes`,
  `upset_intersections` or `gene_scratchpad` (default: `de_genes`).
  - `de_genes`: top DE genes from selected comparison
  - `upset_intersections`: genes contained in specific upset intersection
  - `gene_scratchpad`: genes selected in `Gene scratchpad`
- `Comparison`: The comparison to choose DE genes from. This is shown only
  for `de_genes`.
- `Direction of change`
  - `up & down`: all changed genes, both up- and down-regulated
  - `up`: upregulated genes
  - `down`: downregulated genes
- When `Geneset` is `de_genes` or `upset_intersections`, the current FDR & log2FC thresholds
  are shown here (These can be changed from the `Summary`, `Table` or `Upset plot` tabs).
