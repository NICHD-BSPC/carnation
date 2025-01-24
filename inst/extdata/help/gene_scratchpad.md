#### Gene scratchpad
--------------------

General-purpose scratch space to make note of genes of interest from plots or tables.

- This is accessible from most areas of the app. Genes selected here will be shown on the
  `Gene plot` and labeled on the `MA plot` & `Heatmap`.
- `Quick add`: controls to allow quick addition of genes to the scratchpad.
  - `Comparison`: which comparison to choose DE genes from.
  - `# of genes`: number of top genes to add (default: `6`).
  - `Add top genes by LFC` will add the top up- and down-regulated genes from the
    selected comparison (based on `log2FoldChange`).
  - `Add top genes by padj` adds the top genes with the lowest `padj` values.
- Genes can be added/removed in several other ways:
  - Select rows of the DE analysis `Table` and click `Add to scratchpad` or `Remove from scratchpad` button.
  - Select rows of the UpSet plot table and click `Add to scratchpad` or `Remove from scratchpad` button.
  - Copy & paste Comma-separated gene names here. Genes of interest can be selected from a functional enrichment table
    and added here in this manner.
  - Gene names can also be typed in directly. Autocomplete options will be shown as you type.

