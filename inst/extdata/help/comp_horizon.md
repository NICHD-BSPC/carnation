#### Horizon
------------

Side-by-side visualization of enriched terms from two different comparisons with connecting lines.
- This plot provides a clear visualization of how enrichment patterns differ between two conditions.
- The connecting lines make it easy to track specific terms across comparisons, while the side-by-side layout facilitates direct comparison of significance levels.

**What it shows:**
- Direct comparison of term significance between two conditions
- Connected lines showing shared terms between comparisons
- Relative ranking and significance of terms in each comparison

**When to use it:**
- To directly compare enrichment results between two conditions
- To visualize differences in significance for shared terms
- To identify terms unique to each comparison

**How to interpret:**
- Left side represents enrichment results from `Comparison 1`
- Right side represents enrichment results from `Comparison 2`
- Terms are connected by solid lines if they appear in both comparisons
- The position along the x-axis represents significance (-log10 p-value)
- Higher positions indicate greater statistical significance
- The slope of connecting lines indicates difference in significance between comparisons

**Plot options:**
- `# of terms`: Control how many terms appear in the plot (default: 20)
- `max name length`: Limit the character length of term names for better display
- `sort by`: Choose between:
  - `clustered`: Orders terms according to which comparison has lower p-values
  - `first_set`: Gives precedence to `Comparison 1`
