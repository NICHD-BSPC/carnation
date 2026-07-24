### Gene Plot Controls
----------------------

#### Sample options

- `sample group`: which data subset should be used for the plot?
  - If a particular sample group was used for the selected comparison, then the gene plot uses it by default
  - Alternatively, we can use data from `all_samples`

- `normalization`: choose normalization method
  - `library size` (default): normalizes by sequencing depth
  - `vst` (variance stabilizing transformation): reduces dependence of variance on mean

#### Plot options

- `x-axis variable`: which variable should be shown on the x-axis?
  - All variables available in the metadata (including custom columns) can be used here
  - The variable `group` (if present) is selected by default

- `color by`: how should the samples be colored?
  - All variables in the metadata plus `gene` (default) are options here

- `facet by`: how should the plot be faceted (split)?
  - All variables in the metadata can be selected here. In addition, the plot can be faceted by `gene`
  - By default, this is empty and up to two faceting variables can be selected

- `free y axes`: should the y-axis of the faceted plots be scaled independently?
  - ignored if no faceting variables are set

#### X-axis Settings

- `variable levels`: drag-and-drop control to change the order of x-axis labels or to drop specific levels
  - Click and drag levels from `current` to `unused` to drop them from x-axis view
  - Click and drag in opposite direction to add them back
  - Drag levels up or down to change the order of levels
  - `Select all` moves all levels to `current`, resetting x-axis to the original state
  - `Select none` moves all levels to `unused` and empties the plot

- `rotate x labels`: rotate x-labels to prevent overlapping
  - Use up and down arrows to increment/decrement the rotate angle in steps of 15

#### Facet Settings

- `variable`: dropdown menu showing selected facet variable(s)
  - This menu will not show `gene` as an option, since this can be adjusted using the `Gene scratchpad`

- `levels`: menu to remove specific levels or to change their order
  - `Reset` returns the faceting variable to its original state

#### More Options

- `trendline`: choose line type to connect groups
  - `line` (default): connects group means with straight lines
  - `smooth` : shows smoothed trend

- `boxes`: show boxes on plot to show variability of groups
- `box position`: how should boxes be shown? can be 'identity' (default) where
  boxes can overlap or 'dodge' where boxes are shown side-by-side.
- `legend`: toggle the legend in the gene plot (default: on)
- `text scale`: scaling factor for text on plot (default: 1)
