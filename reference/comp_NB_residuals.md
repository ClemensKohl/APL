# Compute Negative-Binomial residuals

Computes the residuals based on the negative binomial model. By default
a theta of 100 is used to capture technical variation.

## Usage

``` r
comp_NB_residuals(mat, theta = 100, clip = FALSE, cutoff = NULL, freq = TRUE)
```

## Arguments

- mat:

  A numerical matrix or coercible to one by \`as.matrix()\`. Should have
  row and column names.

- theta:

  Overdispersion parameter. By default set to 100 as described in Lause
  and Berens, 2021 (see references).

- clip:

  logical. Whether residuals should be clipped if they are higher/lower
  than a specified cutoff

- cutoff:

  numeric. Residuals that are larger than cutoff or lower than -cutoff
  are clipped to cutoff.

- freq:

  logical. Whether a table of frequencies (as used in CA) should be
  used.

## Value

A named list. The elements are:

- "S": standardized residual matrix.

- "tot": grand total of the original matrix.

- "rowm": row masses.

- "colm": column masses.

## References

Lause, J., Berens, P. & Kobak, D. Analytic Pearson residuals for
normalization of single-cell RNA-seq UMI data. Genome Biol 22, 258
(2021). https://doi.org/10.1186/s13059-021-02451-7
