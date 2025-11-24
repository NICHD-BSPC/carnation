# Get data areas a user has access to

This function takes a username and returns a list with two elements:

## Usage

``` r
check_user_access(al, u, admin = "admin")
```

## Arguments

- al:

  list with access settings; should have two elements - user_group &
  data_area

- u:

  user name

- admin:

  Admin user group

## Details

user_group: one element vector data_area: vector of data areas
