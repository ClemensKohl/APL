# Find most variable rows

Calculates the variance of the chi-square component matrix and selects
the rows with the highest variance, e.g. 5,000.

## Usage

``` r
var_rows(mat, residuals = "pearson", top = 5000, ...)
```

## Arguments

- mat:

  A numeric matrix. For sequencing a count matrix, gene expression
  values with genes in rows and samples/cells in columns. Should contain
  row and column names.

- residuals:

  character string. Specifies which kind of residuals should be
  calculated. Can be "pearson" (default), "freemantukey" or "NB" for
  negative-binomial.

- top:

  Integer. Number of most variable rows to retain. Default 5000.

- ...:

  Further arguments for \`calc_residuals\`.

## Value

Returns a matrix, which consists of the top variable rows of mat.

## Examples

``` r
set.seed(1234)

# Simulate counts
cnts <- mapply(function(x){rpois(n = 500, lambda = x)},
              x = sample(1:20, 50, replace = TRUE))
rownames(cnts) <- paste0("gene_", 1:nrow(cnts))
colnames(cnts) <- paste0("cell_", 1:ncol(cnts))

# Choose top 5000 most variable genes
cnts <- var_rows(mat = cnts, top = 5000)
#> Warning: Top is larger than the number of rows in matrix. Top was set to nrow(mat).

```
