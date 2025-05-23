### Heatmap Controls
--------------------

#### Geneset Selection

- `Geneset`: the type of gene set being visualized. Can be:
  - `de_genes`: top DE genes from selected comparison
  - `upset_intersections`: genes contained in specific upset intersections
  - `gene_scratchpad`: genes selected in `Gene scratchpad`

- If `upset_intersections` is selected, the `Intersections` dropdown menu is shown.
  - The options shown here match the current intersections in the UpSet plot.
- If `de_genes` is selected, the `Comparison` and `Direction of change` menus are shown.
  - `Comparison`: The comparison to choose DE genes from.
  - `Direction of change`: filter genes based on expression direction, `up & down` (default), `up` or `down`.
- When `Geneset` is `de_genes` or `upset_intersections`, the current FDR & log2FC thresholds
  are shown here (these can be changed from the global settings menu in the sidebar).

#### Clustering Options

- `sample group`: sample subset to be used for plotting heatmap (default: `all_samples`).

- `cluster by`: how should the heatmap be clustered?
  - If `row` (default), the rows are clustered while `col` clusters columns. `both` clusters both rows and columns.
    Note that column dendrograms are not shown in the heatmap.
  - `none` turns off any clustering.

- `scale`: should the data be scaled?
  - `row` (default), `column`: scales each row (gene) or `column` (sample) to have mean=0 and standard deviation=1
  - `none`: no scaling, shows raw normalized expression values

#### Display Options

- `# genes to plot`: Number of top DE genes to show
  - If the selected gene set contains more than 150 genes, only the top 150 will be shown by default.
    Note that this limit can be changed using the `max # of genes` control in the `Advanced options` section.
  - For `de_genes`, genes are ranked by the selected ranking metric
  - For `upset_intersections`, genes are ranked by variability across samples

- `ranking metric`: Column used to rank DE genes. Can be:
  - `padj`: rank by statistical significance (default)
  - `log2FoldChange`: rank by magnitude of change
  - This control is only shown when `Geneset` is `de_genes`.

#### Column Settings

- `labels`: how should the columns be labeled?
  - Any column in the metadata (including custom columns) can be used for labeling columns
  - Multiple columns can be selected for more detailed labels
  - In case of duplicate labels, numeric prefixes (e.g. `1_`, `2_`) are added to make them distinct

- `group by`: controls for editing column order
  - Select any metadata column to group columns by that variable
  - Upon selection, columns are ordered by that variable and its levels are shown below
  - Add/remove levels to include/exclude corresponding samples from the heatmap
  - `Reset` returns the data to its initial state

#### More Options

- `row font` and `column font`: adjust the font size for row and column labels
