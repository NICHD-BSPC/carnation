### DE Analysis Filters
-----------------------

- `FDR threshold`: false discovery rate (FDR) threshold (default: 0.1)
  - Adjusts the statistical significance cutoff for calling differentially expressed genes
  - Lower values are more stringent (fewer genes pass the filter)
  - Higher values are more lenient (more genes pass the filter)
  - Commonly used values range from 0.01 (very stringent) to 0.1 (more inclusive)

- `log2FC threshold`: log-2 fold-change (FC) threshold (default: 0)
  - Sets the minimum absolute fold change required for a gene to be considered differentially expressed
  - Higher values focus on genes with larger expression changes
  - Setting to 0 includes all statistically significant genes regardless of fold change magnitude
  - Commonly used values range from 0.5 to 2 (corresponding to 1.4x to 4x fold change)

These filter controls are shared by the following `DE analysis` tabs:
- `Summary`: affects the counts of up/down regulated genes
- `Table`: filters the genes displayed in the table
- `MA plot`: highlights significant genes and affects the DE gene count display
- `Gene plot`: affects which genes are available for selection
- `Heatmap`: determines which genes are included in the heatmap

**Note:**
- Adjusting these thresholds immediately updates all DE analysis visualizations
- For exploratory analysis, consider using more lenient thresholds (higher FDR, lower log2FC)
- For candidate gene selection, use more stringent thresholds (lower FDR, higher log2FC)
- The total number of genes passing filters is displayed in various visualizations
- These global filters provide consistency across all DE analysis views
- Changes to these filters are not saved between sessions

