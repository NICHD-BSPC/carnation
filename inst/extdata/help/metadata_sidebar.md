### Metadata Controls
---------------------

#### Sample Group Selection

- `Sample group`: sample groups that exist in the analysis.
  - `all_samples` denotes a group containing all samples.
  - If samples were divided into groups for different comparisons
    then these groups will be shown here.

#### Edit Metadata

- Metadata can be edited in three ways:
  - `Duplicate column`

    This duplicates an existing column. Users can select a custom column name
    or else the selected column name with a `_dup` suffix is used for the new
    column.
  - `Add empty column`

    Creates a new empty column with a default name (e.g., `user_col1`) or a custom name.
  - `Remove column`

    Allows deletion of one or more custom columns.

- Existing metadata is considered immutable and therefore, *only*
  custom columns are editable.
  - To edit, double click on a cell and enter the desired value.
  - When finished, hit `Enter` while holding down
    the `control` key.

    **Warning**: For multi-page tables, hit `Control+Enter` *before*
    switching pages in order not to lose edits.

- `Apply changes` will propagate edits through the app.
  - Now, custom columns will show up as an option in PCA plot, Gene plot & heatmap.
- To reject changes and reset to the original metadata hit the
  `Reset` button.

- Edits in metadata are only saved during the current
  session. To save changes to metadata, hit `Save changes` in the app header.
