# Carnation

Interactive shiny dashboard for exploring RNA-Seq analysis.

## Usage

``` r
run_carnation(credentials = NULL, passphrase = NULL, enable_admin = TRUE, ...)
```

## Arguments

- credentials:

  path to encrypted sqlite db with user credentials.

- passphrase:

  passphrase for credentials db.

- enable_admin:

  if TRUE, admin view is shown. Note, this is only available if
  credentials have sqlite backend.

- ...:

  parameters passed to shinyApp() call

## Value

shinyApp object

## Examples

``` r
if(interactive()){
    shiny::runApp(
        run_carnation()
    )
}
```
