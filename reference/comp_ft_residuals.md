# Compute Freeman-Tukey residuals

Computes Freeman-Tukey residuals

## Usage

``` r
comp_ft_residuals(mat)
```

## Arguments

- mat:

  A numerical matrix or coercible to one by \`as.matrix()\`. Should have
  row and column names.

## Value

A named list. The elements are:

- "S": standardized residual matrix.

- "tot": grand total of the original matrix.

- "rowm": row masses.

- "colm": column masses.
