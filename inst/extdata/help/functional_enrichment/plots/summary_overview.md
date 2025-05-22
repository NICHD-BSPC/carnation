#### Summary Overview
---------------------

Horizontal bar chart depicting the most significantly enriched terms in your analysis.

**What it shows:**
- The most statistically significant functional terms from your enrichment analysis
- The direction of change (z-score) or magnitude of effect (aggregation score)
- Clear ranking of terms by significance

**When to use it:**
- To quickly identify the most important biological processes or pathways
- To create publication-ready visualizations of your top enrichment results
- To compare significance levels across multiple terms

**How to interpret:**
- The x-axis shows `-log10(p-value)` - higher values indicate greater statistical significance
- The color of circles represents either:
  - `z_score`: Direction of change, calculated as (u-d)/√(u+d), where u=upregulated genes, d=downregulated genes
  - `aggr_score`: Magnitude of effect, calculated as the mean log fold change of genes in the term

**Plot options:**
- `# of terms`: Control how many terms appear in the plot (default: 20)
- `p-value column`: Choose between raw p-values (`pvalue`) or adjusted p-values (`p.adjust`)
- `color by`: Select whether to color by direction (`z_score`) or magnitude (`aggr_score`)
- `max name length`: Limit the character length of term names for better display (default: 50)
