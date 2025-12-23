# Create carnation python environment

This function installs 'plotly' and 'kaleido' python packages in an
environment to allow PDF downloads from plotly plots.

## Usage

``` r
install_carnation(envname, ...)
```

## Arguments

- envname:

  name of the python environment

- ...:

  parameters passed to reticulate::py_install

## Value

`NULL`, invisibly. The function is called for its side effects.

## Examples

``` r
if(interactive()){
    install_carnation()
}
```
