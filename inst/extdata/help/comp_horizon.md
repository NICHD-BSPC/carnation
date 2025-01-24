#### Horizon
------------

An alternative to the `summary_overview` plot where functional terms shared by two comparisons
are plotted side-by-side with terms corresponding to each comparison connected by solid lines.

**Plot options**
- `# of terms` - Number of terms shown in the plot (default=20).
- `sort by`: define sorting order for the plot. Options are:
  - `clustered`: order top shared terms according to comparisons for which they have
    lower p-values.
  - `first_set`: order shared terms according to `Comparison 1` (similar to summary overview).
- `max name length`: Maximum number of characters for term names shown on y-axis (default=50).
