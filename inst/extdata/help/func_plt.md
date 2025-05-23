### Enrichment Plot Controls
-----------------------------

#### Main Controls

- `Comparison`: differential expression comparison that yielded the DE gene set being analyzed

- `Direction`: direction of change in the DE genes
  - Can be `up & down` (default) for all changed genes, `up` for upregulated genes, or `down` for downregulated genes

- `Database`: annotation database used for the analysis
  - Options may include Gene Ontology (GO), KEGG pathways, etc. depending on what was used during the analysis

#### Available Visualizations

- `Type of plot`: select the visualization type for enrichment results
  Available options include:

  - `summary_overview`: bar chart of top enriched terms
    - Provides clear ranking of terms by significance
    - Terms are colored by z-score or aggregation score

  - `enrichment_map`: network visualization of related terms
    - Shows functional terms as nodes and their relationships as edges
    - Terms with shared genes cluster together
    - Node size represents gene count, edge thickness represents gene overlap

  - `cnetplot`: network of genes and their associated terms
    - Shows connections between genes and functional categories
    - Identifies genes that link multiple biological processes
    - Visualizes the specific genes driving enrichment results

  - `radar`: circular visualization of enrichment significance
    - Displays terms around a circle with significance shown as distance from center
    - Useful for comparing significance across multiple terms
    - Provides a compact view of many terms simultaneously

  - `alluvial`: Sankey diagram showing gene-term relationships
    - Visualizes how genes are shared across multiple terms using a flow-based layout
    - Useful for identifying multi-functional genes

  - `dendrogram`: hierarchical clustering tree of related terms
    - Branches represent clusters of related biological processes
    - Useful for identifying major functional themes

  - `emap_distill`: network view of the Distill enrichment analysis
    - Shows clustered terms based on gene overlap similarity
    - Different colors represent different clusters
    - Helps reduce redundancy in enrichment results

  - `emap_fuzzy`: network view of the Fuzzy clustering analysis
    - Shows fuzzy clusters with terms that can belong to multiple clusters
    - Representative terms are outlined in black
    - Useful for understanding complex biological responses

#### Plot Options

- Each plot type has specific options available in the `Plot options` section
  - These options control the appearance and behavior of the current visualization
  - Click `Refresh plot` to apply changes to plot options

**Note:**
- Changes to the main controls (comparison, direction, database) are immediately reflected in both the plot and table views
- Use the plot-specific options to customize the visualization after selecting the main controls
