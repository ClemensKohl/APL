# Calculate correspondence analysis row and column coordinates.

\`ca_coords\` calculates the standardized and principal coordinates of
the rows and columns in CA space.

## Usage

``` r
ca_coords(caobj, dims = NULL, princ_coords = 3, princ_only = FALSE)
```

## Arguments

- caobj:

  A "cacomp" object as outputted from \`cacomp()\`.

- dims:

  Integer indicating the number of dimensions to use for the calculation
  of coordinates. All elements of caobj (where applicable) will be
  reduced to the given number of dimensions. Default NULL (keeps all
  dimensions).

- princ_coords:

  Integer. Number indicating whether principal coordinates should be
  calculated for the rows (=1), columns (=2), both (=3) or none (=0).
  Default 3.

- princ_only:

  Logical, whether only principal coordinates should be calculated. Or,
  in other words, whether the standardized coordinates are already
  calculated and stored in \`caobj\`. Default \`FALSE\`.

## Value

Returns input object with coordinates added.
std_coords_rows/std_coords_cols: Standardized coordinates of
rows/columns. prin_coords_rows/prin_coords_cols: Principal coordinates
of rows/columns.

## Details

Takes a "cacomp" object and calculates standardized and principal
coordinates for the visualization of CA results in a biplot or to
subsequently calculate coordinates in an Association Plot.

## Examples

``` r
# Simulate scRNAseq data.
cnts <- data.frame(cell_1 = rpois(10, 5),
                   cell_2 = rpois(10, 10),
                   cell_3 = rpois(10, 20))
rownames(cnts) <- paste0("gene_", 1:10)
cnts <- as.matrix(cnts)

# Run correspondence analysis.
ca <- cacomp(obj = cnts, princ_coords = 1)
#> Warning: 
#> Parameter top is >nrow(obj) and therefore ignored.
#> No dimensions specified. Setting dimensions to: 2
#> Please consider setting the dimensions to a lower value to speed up the calculation.
#> Recommended dimensionality: << min(nrows, ncols) * 0.2
ca <- ca_coords(ca, princ_coords = 3)
```
