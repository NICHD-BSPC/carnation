#### PCA Plot
-------------

Interactive Principal Component Analysis (PCA) visualization of sample relationships based on gene expression.
This is a useful tool for identifying potential outliers and assessing the overall quality of your samples.

**What it shows:**
- Sample clustering based on gene expression patterns
- Major sources of variation in your dataset
- Relationships between experimental conditions

**When to use it:**
- To assess sample quality and identify potential outliers
- To visualize separation between experimental groups
- To explore sources of variation beyond your experimental design

**How to interpret:**
- Samples that cluster together have similar overall expression profiles
- Distance between points represents degree of difference in expression
- The percentage on each axis shows how much variation that principal component explains in the data
- PC1 (first principal component) captures the largest source of variation in the data

**Interactive features:**
- Color points by any metadata variable
- Select different principal components to visualize (PC1-PC6)
- Selecting a third principal component creates a 3D plot
- Use the `subset samples` controls to filter samples based on metadata variables. This recomputes the PCA plot
  with only the selected samples.
- Use the `gene loadings` controls to overlay gene loadings on the plot, showing the contribution of each gene
  to the principal components.
- Hover over points to see sample details
- Download publication-ready versions of the plot

**Note:** PCA is performed on variance-stabilized or regularized log-transformed data, which helps normalize for sequencing depth and reduces the influence of highly variable genes.
