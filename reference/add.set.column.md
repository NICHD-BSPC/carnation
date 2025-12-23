# Add set column to UpSet plot matrix

This function adds a column denoting set number to a matrix generated
for an upset plot with fromList.with.names()

## Usage

``` r
add.set.column(df)
```

## Arguments

- df:

  binary matrix where row = genes & columns are gene sets, with 1
  indicating that a gene is present is that gene set and vice-versa

## Value

data.frame with added set column

## Examples

``` r
# list of genes
lst <- list(group1 = c(a = "gene1", b = "gene2", c = "gene3", d = "gene4"),
            group2 = c(c = "gene3", d = "gene4"))

# binarized matrix with group membership
df <- fromList.with.names(lst)

# matrix with added set column
ldf <- add.set.column(df)
```
