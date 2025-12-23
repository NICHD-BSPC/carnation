# Generate upset plot table

Generate upset plot table

## Usage

``` r
get_upset_table(gene.lists, comp_split_pattern = ";")
```

## Arguments

- gene.lists:

  list with character vectors of gene names

- comp_split_pattern:

  character used to separate gene set names

## Value

list with upset table elements

## Examples

``` r
lst <- list(group1 = c(a = "gene1", b = "gene2", c = "gene3", d = "gene4"),
            group2 = c(b = "gene2", d = "gene4", e = "gene5"),
            group3 = c(d = "gene4", e = "gene5", f = "gene6"))

df <- get_upset_table(lst)
str(df)
#> List of 3
#>  $ tbl        :'data.frame': 6 obs. of  7 variables:
#>   ..$ symbol     : chr [1:6] "gene1" "gene3" "gene6" "gene5" ...
#>   ..$ set        : chr [1:6] "set1" "set1" "set2" "set3" ...
#>   ..$ group1     : int [1:6] 1 1 0 0 1 1
#>   ..$ group2     : int [1:6] 0 0 0 1 1 1
#>   ..$ group3     : int [1:6] 0 0 1 1 0 1
#>   ..$ comparisons: chr [1:6] "group1" "group1" "group3" "group2;group3" ...
#>   ..$ degree     : int [1:6] 1 1 1 2 2 3
#>  $ set_mapping:List of 5
#>   ..$ set1: chr "group1"
#>   ..$ set2: chr "group3"
#>   ..$ set3: chr [1:2] "group2" "group3"
#>   ..$ set4: chr [1:2] "group1" "group2"
#>   ..$ set5: chr [1:3] "group1" "group2" "group3"
#>  $ set_labels : Named chr [1:5] "set1" "set2" "set3" "set4" ...
#>   ..- attr(*, "names")= chr [1:5] "set1 (group1; n = 2)" "set2 (group3; n = 1)" "set3 (group2, group3; n = 1)" "set4 (group1, group2; n = 1)" ...
```
