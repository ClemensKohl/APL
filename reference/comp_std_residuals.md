# Compute Standardized Residuals

\`comp_std_residuals\` computes the standardized residual matrix S based
on the Poisson model, which is the basis for correspondence analysis and
serves as input for singular value decomposition (SVD).

## Usage

``` r
comp_std_residuals(mat, clip = FALSE, cutoff = NULL)
```

## Arguments

- mat:

  A numerical matrix or coercible to one by \`as.matrix()\`. Should have
  row and column names.

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

## Details

Calculates standardized residual matrix S from the proportion matrix P
and the expected values E according to \\S = \frac{(P-E)}{sqrt(E)}\\.
