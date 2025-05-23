#### Upset Intersection Table
-----------------------------

Interactive table displaying genes from selected intersections in the UpSet plot.

**What it shows:**
- Genes present in specific intersections between differential expression sets
- Log2 fold change and adjusted p-values for each gene across selected comparisons
- Intersection membership and degree information

**When to use it:**
- To explore genes shared across specific comparisons
- To examine detailed statistics for genes in intersections of interest
- To identify and select genes for further analysis

**How to interpret:**
- Each row represents a gene present in the selected intersection(s)
- The `gene` column shows gene IDs
- The `symbol` column shows gene symbols
- The `set` column indicates which intersection the gene belongs to
- The `comparisons` column lists which comparisons share this gene
- The `degree` column shows in how many comparisons the gene is differentially expressed
- Columns named for each comparison show 0 or 1, indicating whether the gene is DE in that comparison
- Additional columns show `log2FoldChange` and `padj` values from each comparison grouped under
  the comparison name

**Interactive features:**
- Select specific intersections using the `Intersections` dropdown menu
- Sort the table by any column by clicking the column header
- Search for specific genes using the search box
- Select genes to add to your Gene scratchpad
- Reset your selection using the `Reset selection` button

**Note:** Changes to the FDR and LFC thresholds affect all DE analysis visualizations.
