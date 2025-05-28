#### DE Analysis Table
---------------------

Interactive table displaying differential expression results for all genes in selected comparisons.

**What it shows:**
- Gene identifiers and descriptions
- Log2 fold changes and statistical significance values
- Base mean expression levels across selected comparisons

**When to use it:**
- To explore the complete differential expression results
- To search for specific genes of interest
- To export gene lists for downstream analysis

**How to interpret:**
- Lower `padj` values indicate higher statistical significance
- Positive `log2FoldChange` values indicate upregulation
- Negative `log2FoldChange` values indicate downregulation
- The `baseMean` column shows the average expression level across all samples in the selected comparisons

**Interactive features:**
- Click on a row to add/remove that gene to/from the Gene scratchpad
- Sort by any column by clicking the column header
- Multi-column sorting by holding Shift while clicking additional column headers
- Filter genes using the search box in the top-right corner of the table
- Toggle between showing only DE genes or all genes using the checkbox
- Adjust FDR and LFC thresholds using the sidebar controls

**Note:** Changes to the FDR and LFC thresholds affect all DE analysis visualizations.
