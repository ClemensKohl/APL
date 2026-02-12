# Recompute missing values of cacomp object.

The caobj needs to have the std_coords_cols, the prin_coords_rows and D
calculated. From this the remainder will be calculated. Future updates
might extend this functionality.

## Usage

``` r
recompute(calist, mat, ...)
```

## Arguments

- calist:

  A list with std_coords_cols, the prin_coords_rows and D.

- mat:

  A matrix from which the cacomp object is derived from.

- ...:

  Further arguments forwarded to cacomp.

## Value

A cacomp object with additional calculated row_masses, col_masses,
std_coords_rows, U and V.
