# Subset dimensions of a caobj

Subsets the dimensions according to user input.

## Usage

``` r
subset_dims(caobj, dims)
```

## Arguments

- caobj:

  A caobj.

- dims:

  Integer. Number of dimensions.

## Value

Returns caobj.

## Examples

``` r
# Simulate scRNAseq data.
cnts <- data.frame(cell_1 = rpois(10, 5),
                   cell_2 = rpois(10, 10),
                   cell_3 = rpois(10, 20))
rownames(cnts) <- paste0("gene_", 1:10)
cnts <- as.matrix(cnts)

# Run correspondence analysis.
ca <- cacomp(cnts)
#> Warning: 
#> Parameter top is >nrow(obj) and therefore ignored.
#> No dimensions specified. Setting dimensions to: 2
#> Please consider setting the dimensions to a lower value to speed up the calculation.
#> Recommended dimensionality: << min(nrows, ncols) * 0.2
ca <- subset_dims(ca, 2)
```
