#### Volcano Plot

------------
Interactive visualization that plots statistical significance against effect size, making it easy to spot genes that are both strongly *and* confidently changed.

**What it shows:**
- Log2 fold change (`log2FoldChange`) on the x-axis — the direction and magnitude of change
- Negative log10 adjusted p-value on the y-axis — genes higher up are more statistically significant
- Up-regulated genes sit to the right, down-regulated to the left, forming the characteristic two-winged "volcano" shape

**When to use it:**
- To find genes that combine a large fold change with strong statistical support (the upper corners)
- To weigh effect size against significance when prioritizing candidates for follow-up
- To check whether up- and down-regulation are roughly balanced, which can flag normalization or design issues

**How to interpret:**
- The most compelling hits sit in the upper-left and upper-right corners: large change *and* low p-value
- Highly significant genes near the center have small fold changes — often stably, highly expressed genes with low variance
- Triangles at the edges indicate genes with fold changes or p-values beyond the axis limits
- Genes from your "Gene scratchpad" are highlighted with dark circles
- **Significance mode:** red points pass your FDR and LFC thresholds; grey points do not
- **Base mean mode:** points are colored by average expression (see the color bar at right), useful for checking whether your top hits are being driven by highly expressed genes

**Interactive features:**
- Hover over points to see gene names
- Switch between significance and base mean coloring
- Adjust FDR and LFC thresholds using the sidebar controls
- Modify axis limits to zoom in on specific fold-change or significance ranges

**Download Options**
- Click the `Download` button to save a publication-ready PDF of the static plot
- The download preserves all current settings:
  - Axis limits
  - FDR and LFC thresholds
  - Gene highlights
  - Color mode and scheme

**Note:** Changes to the FDR and LFC thresholds affect all DE analysis visualizations.
