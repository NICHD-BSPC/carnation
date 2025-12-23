# Prepare list for UpSet plots, but include rownames

Prepare list for UpSet plots, but include rownames

## Usage

``` r
fromList.with.names(lst)
```

## Arguments

- lst:

  List of sets to compare (same input as to UpSetR::fromList)

## Value

data.frame of 1 and 0 showing which genes are in which sets

## Examples

``` r
# list of genes
lst <- list(group1 = c(a = "gene1", b = "gene2", c = "gene3", d = "gene4"),
            group2 = c(c = "gene3", d = "gene4"))

# binarized matrix with group membership
df <- fromList.with.names(lst)
```
