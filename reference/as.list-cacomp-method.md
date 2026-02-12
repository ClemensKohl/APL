# Convert cacomp object to list.

Convert cacomp object to list.

## Usage

``` r
# S4 method for class 'cacomp'
as.list(x)
```

## Arguments

- x:

  A cacomp object.

## Value

A cacomp object.

## Examples

``` r
# Simulate counts
cnts <- mapply(function(x){rpois(n = 500, lambda = x)},
               x = sample(1:100, 50, replace = TRUE))
rownames(cnts) <- paste0("gene_", 1:nrow(cnts))
colnames(cnts) <- paste0("cell_", 1:ncol(cnts))

# Run correspondence analysis
ca <- cacomp(obj = cnts, princ_coords = 3)
#> Warning: 
#> Parameter top is >nrow(obj) and therefore ignored.
#> No dimensions specified. Setting dimensions to: 9
ca_list <- as.list(ca)
```
