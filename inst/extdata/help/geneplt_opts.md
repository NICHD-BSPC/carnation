#### Gene plot controls
-----------------------

Main options to control gene plot:

- `sample group`: which data subset should be used for the plot?
  - If a particular sample group was used for the selected
    comparison, then the gene plot uses it by default.
  - Alternatively, we can use data from `all_samples`.
- `x-axis variable`: which variable should be shown on the x-axis?
  - All variables available in the metadata (including custom
    columns) can be used here.
  - The variable `group` (if present) is selected by default.
- `facet by`: how should the plot be faceted (split)?
  - Again, all variables in the metadata can be selected here.
    In addition, the plot can be faceted by `gene`.
  - By default, this is empty and upto two faceting variables
    can be selected.
- `color by`: how should the samples be colored?
  - Again, all variables in the metadata plus `gene` (default)
    are options here.

