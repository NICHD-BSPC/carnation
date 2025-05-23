### Scatter Plot Controls
-------------------------

#### Comparison Selection

- `Comparison 1 (x-axis)`: select which differential expression comparison to display on the x-axis

- `Comparison 2 (y-axis)`: select which differential expression comparison to display on the y-axis

- `Swap comparisons`: button that exchanges the x and y axis comparisons
  - Useful for quickly changing perspective without manually reselecting comparisons
- The plot will automatically update when any of these controls are changed

#### Plot Options

- `Values to use`: choose which values to compare between the two contrasts
  - `LFC`: compares log2 fold change values (default)
  - `P-adj`: compares -log10 transformed adjusted p-values

- `Interactive?`: toggle between interactive and static plot modes
  - `yes`: creates an interactive plot with hover functionality (default)
  - `no`: creates a static plot suitable for publication

- `Show table?`: control whether to display the comparison data table below the plot
  - `yes`: shows a sortable table with detailed statistics (default)
  - `no`: hides the table to focus on the visualization

#### Plot Settings

- `Point size`: adjust the size of points in the plot
  - Larger values make individual genes more visible
  - Smaller values reduce overlap in dense regions

- `Point opacity`: control the transparency of points
  - Lower values help visualize overlapping points
  - Higher values make individual points more distinct

- `Show grid?`: toggle grid lines on the plot
  - Grid lines can help with precise value estimation

- `Show quadrant lines?`: toggle the display of lines dividing the plot into quadrants
  - These lines help identify genes with consistent or opposite responses

- `Show significance lines?`: toggle the display of threshold lines
  - These lines indicate the FDR and LFC thresholds used for significance

- `Refresh plot`: button to update the plot with all current settings
  - Use this after making multiple setting changes to apply them all at once
