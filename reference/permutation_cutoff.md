# Calculates permuted association plot coordinates

Calculates matrix of apl coordinates when permuting the original data.

## Usage

``` r
permutation_cutoff(
  caobj,
  mat,
  group = caobj@group,
  dims = caobj@dims,
  reps = 10,
  store_perm = FALSE,
  python = TRUE
)
```

## Arguments

- caobj:

  A "cacomp" object with principal row coordinates and standardized
  column coordinates calculated.

- mat:

  A numeric matrix. For sequencing a count matrix, gene expression
  values with genes in rows and samples/cells in columns. Should contain
  row and column names.

- group:

  Vector of indices of the columns to calculate centroid/x-axis
  direction.

- dims:

  Integer. Number of CA dimensions to retain. Needs to be the same as in
  caobj!

- reps:

  Integer. Number of permutations to perform.

- store_perm:

  Logical. Whether permuted data should be stored in the CA object. This
  implementation dramatically speeds up computation compared to
  \`svd()\` in R.

- python:

  DEPRACTED. A logical value indicating whether to use singular-value
  decomposition from the python package torch.

## Value

List with permuted apl coordinates ("apl_perm") and, a list of saved ca
components ("saved_ca") that allow for quick recomputation of the CA
results. For random_direction_cutoff this saved_ca is empty.
