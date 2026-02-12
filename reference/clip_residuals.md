# Perform clipping of residuals

Clips Pearson or negative-binomial residuals above or below a determined
value. For Pearson (Poisson) residuals it is set by default for 1, for
NB at sqrt(n).

## Usage

``` r
clip_residuals(S, cutoff = sqrt(ncol(S)))
```

## Arguments

- S:

  Matrix of residuals.

- cutoff:

  Value above/below which clipping should happen.

## Value

Matrix of clipped residuals.

## References

Lause, J., Berens, P. & Kobak, D. Analytic Pearson residuals for
normalization of single-cell RNA-seq UMI data. Genome Biol 22, 258
(2021). https://doi.org/10.1186/s13059-021-02451-7
