# Calculate residuals for Correspondence analysis

`calc_residuals` provides optional residuals as the basis for
Correspondence Analysis

## Usage

``` r
calc_residuals(mat, residuals = "pearson", clip = FALSE, cutoff = NULL)
```

## Arguments

- mat:

  A numerical matrix or coercible to one by \`as.matrix()\`. Should have
  row and column names.

- residuals:

  character string. Specifies which kind of residuals should be
  calculated. Can be "pearson" (default), "freemantukey" or "NB" for
  negative-binomial.

- clip:

  logical. Whether residuals should be clipped if they are higher/lower
  than a specified cutoff

- cutoff:

  numeric. Residuals that are larger than cutoff or lower than -cutoff
  are clipped to cutoff.

## Value

A named list. The elements are:

- "S": standardized residual matrix.

- "tot": grand total of the original matrix.

- "rowm": row masses.

- "colm": column masses.
