# format gene names to look pretty in table output

This function works by grouping long lists of genes into groups of a
specified size. Each group is collapsed using commas, while groups are
separated by spaces so that datatable formatting is tricked into
separating space-separated groups and not comma-separated groups

## Usage

``` r
format_genes(g, sep = "\\/", genes.per.line = 6)
```

## Arguments

- g:

  vector of gene names

- sep:

  gene name separator

- genes.per.line:

  number of genes to show in a line
