#### DE results table
---------------------

Table summarizing differential expression (DE) analysis results for
a selected comparison:

Column key:
- `gene`: gene ID
- `symbol`: gene symbol
- `baseMean`: average normalized count across all replicates in
  all conditions
- `log2FoldChange`: fold-difference observed in the comparison
  between treatment and control.
- `pvalue`: raw pvalue from Wald or LRT test performed by
  DESeq2.
  - Empty values indicate that the gene was flagged as an outlier and removed from the analysis.
- `padj`: adjusted p-value or false discovery rate (FDR).
  - Empty cells indicate that the gene was removed by independent
    filtering for having too low counts.

- Of the above, the most important columns are `padj` (which
  indicates statistical significance) and `log2FoldChange` (which
  shows the magnitude of the detected change).
  - `padj` is a stronger metric since low `padj` values always imply statistical
    significance, but a high `log2FoldChange` does not.

Other notes:
- Clicking a row on the table, adds the corresponding gene to the `Gene scratchpad`.
  To deselect, simply click on the row again.
- The table can be filtered by FDR (`FDR threshold`) and log2FoldChange (`log2FC threshold`).
- By default, genes are shown after applying filters. To see all genes,
  uncheck the `Only DE genes` checkbox.
- By default, the table is sorted by FDR (`padj`) with the lowest FDR (most significant genes) at the top.
  - To sort by a different column, simply click on the column name.
  - To sort by multiple columns, click the first column, then click the second column while holding down the Shift key.
- The table can also be filtered by gene names using the search box at the top right hand corner of the table.


