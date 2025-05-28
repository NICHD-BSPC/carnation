### Functional Enrichment Table Controls
-----------------------------------------

#### Main Controls

- `Comparison`: differential expression comparison that yielded the DE gene set being analyzed

- `Direction`: direction of change in the DE genes
  - Can be `up & down` (default) for all changed genes, `up` for upregulated genes, or `down` for downregulated genes

- `Database`: annotation database used for the analysis
  - Options may include Gene Ontology (GO), KEGG pathways, etc. depending on what was used during the analysis

#### Search Table

- `Search in`: select which fields to search within
  - `genes + description`: search in both gene lists and term descriptions (default)
  - `description`: search only in term descriptions
  - `genes`: search only in gene lists

- `Enter text to search`: type or select terms to search for
  - Multiple search terms can be entered. Supports fuzzy matching
  - Click `Refresh` to apply the search filter

- Matching terms are highlighted in **bold** in the table.

#### Subset Table

- `Subset by`: select a method to filter the table based on gene sets
  - `none`: no subsetting applied (default)
  - `gene_scratchpad`: show only terms containing genes from the Gene scratchpad
  - `upset_intersections`: show only terms containing genes from selected UpSet plot intersections

- `Upset intersection`: select which intersection to use for subsetting
  - Only shown when `Subset by` is set to `upset_intersections`
  - Dropdown shows available intersections from the current UpSet plot

- `Refresh`: apply the selected subsetting options to the table
  - Click after making changes to subset settings
  - Table will update to show only terms containing the selected genes

- Matching genes are highlighted in **bold** in the gene column

#### Table Display Options

- `Genes per line`: number of genes to be shown per line in the table (default: 6)
  - Higher values show more genes on a single line while lower values improve readability but make the table taller

- `# of terms`: number of terms to be used for clustering (default: 30)
  - Used for `Distill enrichment` or `Fuzzy clustering` tables
  - If this exceeds the number of rows in the enrichment table, it automatically defaults to the latter

#### Gene Selection

- Select one or more rows in the table and click `Add to scratchpad` to add associated genes to the scratchpad.
  - Useful for collecting genes of interest from enriched pathways
- `Reset selection` clears the current table selection, without affecting genes already added to the scratchpad

**Note:**
- Changes in the main controls are immediately reflected in both the table and plots
- The same controls also affect what is shown in the `Plots` tab
