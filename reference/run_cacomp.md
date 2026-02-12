# Internal function for \`cacomp\`

\`run_cacomp\` performs correspondence analysis on a matrix and returns
the transformed data.

## Usage

``` r
run_cacomp(
  obj,
  coords = TRUE,
  princ_coords = 3,
  python = FALSE,
  dims = 100,
  top = 5000,
  inertia = TRUE,
  rm_zeros = TRUE,
  residuals = "pearson",
  cutoff = NULL,
  clip = FALSE,
  ...
)
```

## Arguments

- obj:

  A numeric matrix or Seurat/SingleCellExperiment object. For sequencing
  a count matrix, gene expression values with genes in rows and
  samples/cells in columns. Should contain row and column names.

- coords:

  Logical. Indicates whether CA standard coordinates should be
  calculated.

- princ_coords:

  Integer. Number indicating whether principal coordinates should be
  calculated for the rows (=1), columns (=2), both (=3) or none (=0).

- python:

  DEPRACTED. A logical value indicating whether to use singular-value
  decomposition from the python package torch. This implementation
  dramatically speeds up computation compared to \`svd()\` in R when
  calculating the full SVD. This parameter only works when dims==NULL or
  dims==rank(mat), where caculating a full SVD is demanded.

- dims:

  Integer. Number of CA dimensions to retain. If NULL: (0.2 \*
  min(nrow(A), ncol(A)) - 1 ).

- top:

  Integer. Number of most variable rows to retain. Set NULL to keep all.

- inertia:

  Logical. Whether total, row and column inertias should be calculated
  and returned.

- rm_zeros:

  Logical. Whether rows & cols containing only 0s should be removed.
  Keeping zero only rows/cols might lead to unexpected results.

- residuals:

  character string. Specifies which kind of residuals should be
  calculated. Can be "pearson" (default), "freemantukey" or "NB" for
  negative-binomial.

- cutoff:

  numeric. Residuals that are larger than cutoff or lower than -cutoff
  are clipped to cutoff.

- clip:

  logical. Whether residuals should be clipped if they are higher/lower
  than a specified cutoff

- ...:

  Arguments forwarded to methods.

## Value

Returns a named list of class "cacomp" with components U, V and D: The
results from the SVD. row_masses and col_masses: Row and columns masses.
top_rows: How many of the most variable rows/genes were retained for the
analysis. tot_inertia, row_inertia and col_inertia: Only if inertia =
TRUE. Total, row and column inertia respectively.

## Details

The calculation is performed according to the work of Michael Greenacre.
When working with large matrices, CA coordinates and principal
coordinates should only be computed when needed to save computational
time.

## References

Greenacre, M. Correspondence Analysis in Practice, Third Edition, 2017.
