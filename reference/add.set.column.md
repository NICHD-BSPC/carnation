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
