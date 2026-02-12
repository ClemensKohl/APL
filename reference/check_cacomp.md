# Check if cacomp object was correctly created.

Checks if the slots in a cacomp object are of the correct size and
whether they are coherent.

## Usage

``` r
check_cacomp(object)
```

## Arguments

- object:

  A cacomp object.

## Value

TRUE if it is a valid cacomp object. FALSE otherwise.

## Examples

``` r
# Simulate scRNAseq data.
cnts <- data.frame(cell_1 = rpois(10, 5),
                   cell_2 = rpois(10, 10),
                   cell_3 = rpois(10, 20))
rownames(cnts) <- paste0("gene_", 1:10)
cnts <- as.matrix(cnts)

# Run correspondence analysis.
ca <- cacomp(obj = cnts, princ_coords = 3, top = 5)
#> No dimensions specified. Setting dimensions to: 2
#> Please consider setting the dimensions to a lower value to speed up the calculation.
#> Recommended dimensionality: << min(nrows, ncols) * 0.2

check_cacomp(ca)
#> [1] TRUE
```
