# Access slots in a cacomp object

Access slots in a cacomp object

## Usage

``` r
cacomp_slot(caobj, slot)
```

## Arguments

- caobj:

  a cacomp object

- slot:

  slot to return

## Value

Chosen slot of the cacomp object

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

# access left singular vectors
cacomp_slot(ca, "U")
#>               Dim1       Dim2
#> gene_8  0.61493584 -0.5591409
#> gene_3  0.33693503  0.5993556
#> gene_2 -0.66907941 -0.1327659
#> gene_7 -0.06190695  0.4432526
#> gene_6 -0.23838828 -0.3376902
```
