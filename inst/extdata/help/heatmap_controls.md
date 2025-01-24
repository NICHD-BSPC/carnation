#### Heatmap controls
---------------------

Main controls to adjust the heatmap

- `sample group`: sample subset
  to be used for plotting heatmap (default: `all_samples`).
- `cluster by`: how should the heatmap be clustered?
  - If `row` (default), the rows are clustered and a dendrogram is
    shown to the right of the heatmap.
  - If `col`, the columns are clustered and a dendrogram is shown above
    the plot.
  - If `both`, both row and columns are clustered and dendrograms generated.
  - `none` turns off any clustering.
- `scale`: should the data be scaled?
  - By default, the data are scaled by `row`.
- `# genes to plot`: Number of top DE genes to show
- `ranking metric`: Column used to rank DE genes. Can be `padj` or `log2FoldChange`
  (default: `padj`).
  - This control is only shown when `Geneset` is `de_genes`.
