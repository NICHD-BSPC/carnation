# Carnation

Interactive shiny dashboard for exploring RNA-Seq analysis.

## Usage

``` r
run_carnation(
  credentials = NULL,
  passphrase = NULL,
  enable_admin = TRUE,
  config_path = NULL,
  ...
)
```

## Arguments

- credentials:

  path to encrypted sqlite db with user credentials.

- passphrase:

  passphrase for credentials db.

- enable_admin:

  if TRUE, admin view is shown. Note, this is only available if
  credentials have sqlite backend.

- config_path:

  optional path to a local config yaml override.

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
