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

project name

## Details

So if the path is: /path/to/project/test/main.*pattern*.rds and depth=2
& end_offset=0 this returns: project/test
