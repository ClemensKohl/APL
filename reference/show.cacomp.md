# Prints cacomp object

Provides more user friendly printing of cacomp objects.

## Usage

``` r
show.cacomp(object)

# S4 method for class 'cacomp'
show(object)
```

## Arguments

- object:

  cacomp object to print

## Value

prints summary information about cacomp object.

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

ca
#> cacomp object with 3 columns, 5 rows and 2 dimensions.
#> Calc. standard coord.:  std_coords_rows, std_coords_cols
#> Calc. principal coord.: prin_coords_rows, prin_coords_cols
#> Calc. APL coord.:       
#> Explained inertia:      66.4% Dim1, 33.6% Dim2
```
