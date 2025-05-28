### UpSet Plot/Table Controls
-----------------------------

#### Main Controls

- `Comparisons`: drag-and-drop control to select/deselect comparisons to compare
  - Click and drag gene sets from `current` to `unused` to drop them from the comparison
  - Click and drag in opposite direction to add them back
  - `Select all` moves all gene sets to `current`, resetting the plot to the starting condition
  - `Select none` moves all gene sets to `unused` and empties the plot
  - On first load, the first 5 comparisons are selected by default
  - Note that at least two comparisons must be selected to generate this plot
  - **Warning**: Moving comparisons to/from *Current* rapidly can cause the app to get stuck in an infinite loop
    - To avoid this behavior, please give a reasonable amount of time between moves
  - `Refresh`: update the plot/table with the current selection of comparisons

#### Shared Options

- `Direction of change`: select which genes to include based on their expression change
  - `up & down (changed)`: include all differentially expressed genes (default)
  - `up`, `down`: include only upregulated or downregulated genes
  - `custom`: create custom selections for each comparison
    - When selected, a `Customize` button appears that opens a custom selection matrix dialog
    - In the matrix, select which direction (up, down, or both) to use for each comparison
    - Click `Apply` to save your selections and `Plot` to update the visualization
    - `Reset` clears all custom selections and `Cancel` closes the dialog without saving any changes
    - The direction is appended to comparison names in the plot for clarity

#### Plot Options

The following options are only shown when viewing the Plot tab:

- `# of intersections`: control the number of intersections shown on the UpSet plot (default: 40)
  - Increasing this value shows more intersections but may make the plot more crowded
  - Decreasing this value focuses on the largest intersections

- `Min intersection size`: set the minimum size of intersections to display (default: 5)
  - Intersections smaller than this value will be hidden
  - Useful for filtering out small, potentially less meaningful intersections

- `Text scale`: adjust the size of text labels on the plot
  - Values > 1 or < 1 make the text larger or smaller, respectively
  - Useful for adjusting text visibility based on plot size

- `Refresh plot`: update the plot with the current display settings
  - Allows you to make multiple setting changes before refreshing the visualization

#### Table Options

The following options are only shown when viewing the Table tab:

- `Intersections`: select specific intersections to view in the Table tab
  - Intersections are labeled based on their order in the UpSet plot
    - For example, the first set is called `set1`, the second set is `set2`, and so on
    - Names include the comparisons sharing these genes and the size of the intersection
    - Example: `set1 (comp_1 comp_2; n = 25)` indicates 25 genes shared by comparisons *comp_1* and *comp_2*
  - By default, the first 10 intersections are selected
  - `Select all` and `Select none` buttons allow quick selection of all or no intersections

