# Get project name from path

This function takes in a path to an RDS file and returns a string to be
used as project name

## Usage

``` r
get_project_name_from_path(
  x,
  depth = 2,
  end_offset = 0,
  staging_dir = "dev",
  fsep = .Platform$file.sep
)
```

## Arguments

- x:

  character path to RDS file

- depth:

  integer how many levels below path to look?

- end_offset:

  integer how far from the end of path to end?

- staging_dir:

  name of staging directory

- fsep:

  file separator to split path with

## Value

project name parsed from path to object

## Examples

``` r
# path to carnation object
obj_path <- "/path/to/project/test/main.rnaseq.rds"

# parsed project name
get_project_name_from_path(obj_path, depth = 2, end_offset = 0)
#> [1] "project/test"
```
