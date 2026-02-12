# Pipe operator

See `magrittr::%>%` for details.

## Usage

``` r
lhs %>% rhs
```

## Arguments

- lhs:

  A value or the magrittr placeholder.

- rhs:

  A function call using the magrittr semantics.

## Value

`magrittr::%>%`

## Examples

``` r
x <- 1:100
x %>% head()
#> [1] 1 2 3 4 5 6
```
