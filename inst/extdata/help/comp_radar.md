#### Radar (Comparison view)
-----------------------------

Radar plot comparing enriched terms between two different comparisons.
- This plot provides an intuitive circular visualization for comparing enrichment patterns between two conditions.
- The overlaid patterns make it easy to identify both similarities and differences in functional enrichment.

**What it shows:**
- Circular visualization of enrichment significance from two comparisons
- Overlaid patterns revealing similarities and differences
- Relative importance of terms across conditions

**When to use it:**
- To compare enrichment patterns between two conditions
- To identify terms with similar or different significance
- To create a compact visualization of comparative enrichment

**How to interpret:**
- Each spoke represents a functional term or pathway
- The length of spokes represents `-log10(pvalue)` - longer spokes indicate greater significance
- Different colors or patterns represent each comparison
- Areas where patterns overlap indicate similar enrichment
- Areas where patterns diverge indicate differential enrichment

**Plot options:**
- `# of terms`: Control how many terms appear in the plot (default: 20)
- `p-value column`: Choose between raw p-values (`pvalue`) or adjusted p-values (`p.adjust`)
- `max name length`: Limit the character length of term names for better display
