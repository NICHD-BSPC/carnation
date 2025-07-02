### Compare Results Controls
---------------------------

#### Main Controls

- `Comparison 1` & `Comparison 2`: differential expression comparisons to be compared
  - Results from these comparisons will be displayed side by side or overlaid

- `Direction`: direction of change in the DE genes
  - Can be `up & down` (default) for all changed genes, `up` for upregulated genes, or `down` for downregulated genes
  - Note that different directions within the same comparison can be compared

- `Database`: annotation database used for the analysis
  - Options may include Gene Ontology (GO), KEGG pathways, etc. depending on what was used during the analysis

- `Swap comparisons`: exchanges Comparison 1 and Comparison 2
  - Useful for quickly reversing the comparison order without manually reselecting

- `Swap directions`: exchanges the selected directions for both comparisons
  - If different directions are selected for each comparison, this swaps them

- The plot will immediately update to reflect any change in these controls

#### Available Visualizations

- `Type of plot`: select the visualization type for comparing enrichment results in the main panel
  Available options include:
  - `summary_overview`: side-by-side bar charts of top enriched terms
    - Allows direct comparison of term significance between conditions
    - Terms are ordered based on significance in `Comparison 1` and may be colored by z-score or aggregate score

  - `radar`: circular visualization comparing term significance
    - Displays terms around a circle with significance shown as distance from center
    - Allows easy identification of terms enriched in one comparison but not the other
    - Overlays results from both comparisons for direct visual comparison

  - `horizon`: horizontal bar chart with connected terms
    - Shows terms shared between comparisons with connecting lines
    - Terms are sorted by significance in either comparison
    - Particularly useful for identifying terms with different significance levels between comparisons

- Each plot type has specific options available in the `Plot options` section to control the appearance and behavior
  of the current plot

- Click `Refresh plot` to apply changes to plot options

**Note:**
- Any change in the main controls is immediately reflected in the plot
- Use the plot-specific options to customize the visualization after selecting the comparisons
