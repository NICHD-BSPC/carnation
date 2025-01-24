#### Column settings
--------------------


- `labels`: how should the columns be labeled?
  - Any column in the metadata (including custom columns) can be
    used for labeling columns. Multiple columns can be selected here.
  - In case of duplicate labels, numeric prefixes (e.g. `1_`, `2_`)
    are added to column labels to make them distinct.
- `group by`: controls for editing column order
  - Select any column in the metadata (including custom columns)
    from the dropdown menu to group columns.
  - Upon selecting a grouping variable, the columns are ordered
    by that variable and levels of the grouping variable are shown
    in the box below.
  - Add/remove levels from this box to remove the corresponding
    samples from the heatmap.
  - `Reset` returns the data to its initial state.
