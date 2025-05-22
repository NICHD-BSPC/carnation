#### Summary Overview (Comparison view)
---------------------------------------

Horizontal bar chart comparing enriched terms between two different comparisons.
- This plot extends the standard summary overview to enable direct comparison between two different enrichment results.
- It's particularly useful for identifying biological processes that respond similarly or differently across experimental conditions.

**What it shows:**
- Side-by-side visualization of enriched terms from two different comparisons
- Relative significance of terms across conditions
- Shared and unique enrichment patterns

**When to use it:**
- To directly compare enrichment results between two conditions
- To identify terms that are consistently or differentially enriched
- To create publication-ready visualizations of comparative enrichment

**How to interpret:**
- Each term has two circles representing significance in each comparison
- Solid black outline circles represent `Comparison 1`
- Semi-transparent circles represent `Comparison 2`
- The x-axis shows `-log10(p-value)` - higher values indicate greater significance
- Terms are ordered by significance in `Comparison 1` by default

**Plot options:**
- `# of terms`: Control how many terms appear in the plot (default: 20)
- `p-value column`: Choose between raw p-values (`pvalue`) or adjusted p-values (`p.adjust`)
- `max name length`: Limit the character length of term names for better display (default: 50)
