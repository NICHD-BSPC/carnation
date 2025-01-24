#### Comparisons

Drag-and-drop control to select/deselect comparisons to compare
- Click and drag gene sets from `current` to `unused` to drop them
  from the comparison. Click and drag in opposite direction to
  add them back.
- `Select all` moves all gene sets to `current`, i.e. resets the plot
  to the starting condition.
- `Select none` moves all gene sets to `unused` and empties the plot.
- On first load, the first 5 comparisons are selected by default.
  Note that, at least two comparisons must be selected to generate this plot.
- **Warning**: Moving comparisons to/from *Current* rapidly, can cause
  the app to get stuck in an infinite loop. This is a known issue and the only
  way to stop this is to refresh the app.
  - To avoid this behavior, please give a reasonable amount of time between moves.

