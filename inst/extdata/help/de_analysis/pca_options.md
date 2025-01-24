#### PCA plot options
---------------------

- `color by`: color samples by metadata columns
  - multiple selections result in samples being colored with all combinations
    of selected columns
- `x-axis`, `y-axis`: choose different principal components (PCs) to
  show on the x- or y-axis.
  - this is useful to explore the data when PC1 & PC2 don't account for
    a large fraction of the sample variation
  - The app currently supports PC1 - PC6
- `sample selection`: select samples to show on the PCA plot using a
  `grouping variable` selected from metadata
  - samples can be added or removed using levels of the selected grouping variable
    to observe the effect on the PCA plot.
