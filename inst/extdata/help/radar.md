#### Radar
----------

Interactive plot depicting functional terms as spokes around a circle.
- The radar plot provides an alternative visualization to the summary overview, allowing you to compare the
  significance of multiple terms in a compact, circular layout.
- This can be particularly useful when you want to highlight the relative importance of different functional categories.

**What it shows:**
- Significance of enriched terms represented as distance from center
- Multiple terms arranged in a circular layout
- Relative importance of different functional categories

**When to use it:**
- To compare significance levels across multiple terms
- To create a compact visualization of top enrichment results
- To generate an alternative view of enrichment significance

**How to interpret:**
- Each spoke represents a functional term or pathway
- The length of spokes represents `-log10(pvalue)` - longer spokes indicate greater significance
- The circular arrangement allows comparison of many terms simultaneously

**Plot options:**
- `# of terms`: Control how many terms appear in the plot (default: 20)
- `p-value column`: Choose between raw p-values (`pvalue`) or adjusted p-values (`p.adjust`)
- `max name length`: Limit the character length of term names for better display