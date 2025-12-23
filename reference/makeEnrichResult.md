# Make an enrichResult obj from a data frame

Most of the parameters are just placeholders and the dataframe must
contain the columns 'ID' and 'geneID'

## Usage

``` r
makeEnrichResult(
  df,
  split = "/",
  keytype = "UNKNOWN",
  ontology = "UNKNOWN",
  type = "enrichResult"
)
```

## Arguments

- df:

  data frame with functional enrichment results

- split:

  string, character used to split gene IDs

- keytype:

  type of gene ID

- ontology:

  ontology database being used

- type:

  string, can be 'enrichResult' or 'gseaResult'

## Value

enrichResult object

## Examples

``` r
# get enrichResult object
data(eres_dex, package='carnation')

# extract the results
df <- as.data.frame(eres_dex)

# convert to a stripped down enrichResult object
eres2 <- makeEnrichResult(df)
```
