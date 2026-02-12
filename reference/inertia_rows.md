# Find most variable rows

Calculates the contributing inertia of each row which is defined as sum
of squares of pearson residuals and selects the rows with the largested
inertias, e.g. 5,000.

## Usage

``` r
inertia_rows(mat, top = 5000, ...)
```

## Arguments

- mat:

  A matrix with genes in rows and cells in columns.

- top:

  Number of genes to select.

- ...:

  Further arguments for \`comp_std_residuals\`

## Value

Returns a matrix, which consists of the top variable rows of mat.
