#### Enrichment Table
---------------------

Interactive table displaying functional enrichment results for differentially expressed genes.

**What it shows:**
- Enriched functional terms or pathways
- Associated genes for each term
- Statistical significance and enrichment metrics

**When to use it:**
- To explore the biological functions associated with your DE genes
- To search for specific pathways or biological processes
- To identify the genes driving specific functional enrichments

**How to interpret:**
- `ID`, `Description`: Identifier and description of the functional term
- `geneID`: Genes associated with the functional term
- `Count`: Number of DE genes associated with the term
- `pvalue`, `p.adjust`: Raw and adjusted p-values from hypergeometric test
- Lower p-values indicate stronger enrichment of the term

**Interactive features:**
- Use the "Search table" controls to filter results by genes and/or description using fuzzy matching
- Alternatively search for specific terms or genes using the search box
- Subset the table by genes in the Gene scratchpad or by UpSet plot intersections
- Select a row and click "Add to scratchpad" to add genes associated with that term to the Gene scratchpad
- Sort by any column to prioritize results
- Export results for further analysis or publication

**Note:** The table displays results based on the currently selected comparison, direction, and database, which can be changed using the controls in the sidebar.
