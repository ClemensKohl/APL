# Calculate Association Plot coordinates

Calculates the Association Plot coordinates for either the rows, columns
or both (default).

## Usage

``` r
apl_coords(caobj, group, calc_rows = TRUE, calc_cols = TRUE)
```

## Arguments

- caobj:

  A "cacomp" object with principal row coordinates and standardized
  column coordinates calculated.

- group:

  Numeric/Character. Vector of indices or column names of the columns to
  calculate centroid/x-axis direction.

- calc_rows:

  TRUE/FALSE. Whether apl row coordinates should be calculated. Default
  TRUE.

- calc_cols:

  TRUE/FALSE. Whether apl column coordinates should be calculated.
  Default TRUE.

## Value

Returns input "cacomp" object and adds components "apl_rows" and/or
"apl_cols" for row and column coordinates. In "group" the indices of the
columns used to calculate the centroid are saved.

## Details

Coordinates (x,y) of row vector \\\vec{r}\\ are defined as
\$\$x(\vec{r}) := \left\|\vec{r}\right\|\cos(\phi(\vec{r}))\$\$
\$\$y(\vec{r}) := \left\|\vec{r}\right\|\sin(\phi(\vec{r}))\$\$ The
x-direction is determined by calculating the centroid of the columns
selected with the indices in "group".

## References

Association Plots: Visualizing associations in high-dimensional
correspondence analysis biplots Elzbieta Gralinska, Martin Vingron
bioRxiv 2020.10.23.352096; doi:
https://doi.org/10.1101/2020.10.23.352096

## Examples

``` r
set.seed(1234)
# Simulate scRNAseq data
cnts <- data.frame(cell_1 = rpois(10, 5),
                   cell_2 = rpois(10, 10),
                   cell_3 = rpois(10, 20),
                   cell_4 = rpois(10, 20))
rownames(cnts) <- paste0("gene_", 1:10)
cnts <- as.matrix(cnts)

# Run correspondence analysis
ca <- cacomp(obj = cnts, princ_coords = 3, dims = 3)
#> Warning: 
#> Parameter top is >nrow(obj) and therefore ignored.
#> Please consider setting the dimensions to a lower value to speed up the calculation.
#> Recommended dimensionality: << min(nrows, ncols) * 0.2
# Calculate APL coordinates
ca <- apl_coords(ca, group = 3:4)
```
