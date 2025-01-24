## What's new?

Below, we summarize recent feature updates grouped by date. Major changes in the
user interface or major feature updates are tagged with <font color=red>NEW</font>.


### 2024-10-02
--------------

#### <font color=red>NEW</font> Custom gene set selection for UpSet plots

Fully customizable selection of gene sets is now possible for UpSet plots.
To use this, go to the UpSet plot tab, click on the settings button
and open up Upset options.

<img src="../www/upset_custom_menu.png" alt="upset custom menu" width="600"/>

Next, select *custom* from the 'Direction of change' menu, and click *Customize*.

<img src="../www/upset_custom_select.png" alt="upset custom select" width="300"/>
<img src="../www/upset_custom_click.png" alt="upset click customize" width="300"/>

Now click on the popup table to choose interesting gene sets and then click *Apply*.

<img src="../www/upset_tbl_select.png" alt="upset table select" width="500"/>

This will update the table to show selected gene sets with check marks.
Next, click *Plot*.

<img src="../www/upset_tbl_apply.png" alt="upset apply selections" width="500"/>

Now the plot will be refreshed to show an upset plot using only the selected gene sets.
To prevent confusion, the direction of change is also appended to the comparison names.

<img src="../www/upset_comp_dir_append.png" alt="upset direction added to comparison" width="600"/>



### 2024-07-18
--------------

#### <font color=red>NEW</font> Search and subset functional enrichment tables

New search and subset options are available for functional enrichment results.

<img src="../www/fe_search_settings.png" alt="search_settings" width="400"/>

**Use simple text to search in genes or description of functional enrichment tables**

<img src="../www/fe_text_search.png" alt="text_search" width="300"/>

Click refresh. <img src="../www/fe_text_search2.png" alt="text_search" width="300"/>

Search is case-insensitive and matches are highlighted in **bold** type. Also, you will now
only see rows of the table matching search terms.

<img src="../www/search_bold.png" alt="search_bold" width="800"/>

**Multiple search terms are combined with logical 'OR'**

For example, searching for "col11" and "cardiac", will filter the tables to show matches to *either*.

<img src="../www/search_multiple.png" alt="search_multiple" width="300"/>

<img src="../www/search_multiple2.png" alt="search_multiple2" width="800"/>

**Functional enrichment plots reflect filtered results**

Before filtering:

No search terms <img src="../www/no_filter.png" alt="no_filter" width="300"/>

Unfiltered plot <img src="../www/no_filter_plot.png" alt="no_filter_plot" width="600"/>

After filtering:

Multiple search terms <img src="../www/yes_filter.png" alt="yes_filter" width="300"/>

Now the plot is filtered <img src="../www/yes_filter_plot.png" alt="yes_filter_plot" width="600"/>

**Filter by gene scratchpad**

You can now filter the functional enrichment tables by genes of interest via the
*Gene scratchpad* or by UpSet plot intersections.

Add your genes of interest to the scratchpad <img src="../www/gene_scratch.png" alt="gene_scratch" width="300"/>

Select gene scratchpad from subset table options.
<img src="../www/subset_opts.png" alt="subset_opts" width="300"/>

Click refresh.
<img src="../www/subset_scratch.png" alt="subset_scratch" width="300"/>

Table will now show rows containing your genes of interest.

<img src="../www/subset_tbl_scratch.png" alt="subset_tbl_scratch" width="800"/>

**Filter by Upset intersection**

You can do a similar filtering using Upset intersections.

First, select an interesting gene set from the Upset plot, e.g. genes changed in both comparisons.

<img src="../www/upset_interesting.png" alt="upset_interesting" width="800"/>

Choose the interesting gene set from the dropdown menu and click 'Refresh'.

<img src="../www/subset_upset.png" alt="subset_upset" width="300"/>

Now the table is filtered by this upset intersection.

<img src="../www/subset_upset_tbl.png" alt="subset_upset_tbl" width="800"/>

**Easily add genes from functional enrichment tables to the gene scratchpad**

Select interesting functional terms from the table and click 'Add to scratchpad'.

<img src="../www/select_fe_rows.png" alt="select_fe_rows" width="800"/>

Now the genes from those rows will be shown in the gene scratchpad.

<img src="../www/select_fe_scratch.png" alt="select_fe_scratch" width="300"/>

Next you can view the expression of these genes using the gene plot.

<img src="../www/select_fe_gene_plot.png" alt="select_fe_gene_plot" width="800"/>

#### Other updates

- *Gene plot*
  - The gene plot now better handles many faceted plots. Instead of manually setting # of rows,
    plots are generated with a fixed number of genes per row: 8 by default
    and 6 if y-axes are 'free' to allow space for y-axis labels. The plot height is automatically
    adjusted to maintain overall viewability.
  - Added option to show boxes
  - Fixed the 'line' trendline option. Now it uses `stat_summary` to draw a line between groups
    instead of joining the individual samples.

- Gene patterns
  - Pattern plots now show number of labeled genes in the facet strip text.
  - The labeled gene table now shows labeled genes in a cluster in a compact single row-view.

### 2024-05-23
--------------

- Optimization and bug fixes

  - Carnation now supports native 'carnation-ready' objects which contain light-weight
    versions of functional enrichment results and gene patterns objects. This minimizes
    pre-processing and vastly improves load times, especially for complex projects.
  - Issues with using gene scratchpad in DE Analysis -> Heatmap are now fixed.
  - The heatmap now respects input gene order via the gene scratchpad if row clustering is turned off.

### 2024-04-17
--------------

#### <font color=red>NEW</font> Gene loadings in PCA

You can now overlay gene loadings onto the PCA plot using the new
`gene loadings` menu in the sidebar.

<img src="../www/pca_gene_loadings.png" alt="pca_loadings" width="800"/>

#### Other updates

- Now you can choose the number of top variable genes to run the PCA on (default: 500).
- Fixed a bug where choosing multiple comparisons in the `DE analysis` > `Table` was
  causing the app to crash due to mismatches in column names.

### 2024-04-04
-------------

#### <font color=red>NEW</font> Scatter Plot

This new plot provides a way of comparing differential gene expression results between
comparisons. In addition, a searchable and sortable table can be displayed below, showing
the plotted data that can be either log2 fold-changes or adjusted p-values.

- Genes are grouped and color-coded into five categories based on statistical significance,
  direction of change and whether they are shared between the comparisons, or unique to each.
- Genes of interest can be labeled on the plot via the gene scratchpad. The plot
  can also be downloaded to PDF.

The scatter plot can be accessed under the *DE Analysis* tab.

<img src="../www/scatter_tab.png" alt="scatter_tab" width="800"/>

#### Functional enrichment

- Functional enrichment plot selector is now moved from the settings dropdown to the main panel.
- Enrichment plots and tables now have a title showing the **comparison|direction|database**

<img src="../www/funenrich_plotsel.jpg" alt="plot_selector" width="800"/>

- Inputs on the *Table* tab of functional enrichment updates *Plots* & vice-versa.
  - This fixes an earlier issue where *Distill enrichment* & *Fuzzy clustering* plots
    were not updating after inputs on the *Table* tab and vice-versa.

#### Other updates

- The *Table* tab of *DE analysis* can now show multiple comparisons at once. This makes it easier to
  search for genes in multiple comparisons simultaneously.
- Recent feature updates will now be highlighted in the *What's new?* tab.
