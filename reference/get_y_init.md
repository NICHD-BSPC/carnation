# Get initial y-axis limits

Get initial y-axis limits

## Usage

``` r
get_y_init(df, y_delta, pseudocount)
```

## Arguments

- df:

  data.frame with counts. Must have column 'count'

- y_delta:

  y-axis padding for visualization, must be between 0 and 1

- pseudocount:

  pseudo-count to add to the data.frame
