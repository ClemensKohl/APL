# Random direction association plot coordinates

Calculates matrix of apl coordinates for random directions

## Usage

``` r
random_direction_cutoff(caobj, dims = caobj@dims, reps = 100)
```

## Arguments

- caobj:

  A "cacomp" object with principal row coordinates and standardized
  column coordinates calculated.

- dims:

  Integer. Number of CA dimensions to retain. Needs to be the same as in
  caobj!

- reps:

  Integer. Number of permutations to perform.

## Value

List with permuted apl coordinates ("apl_perm") and, a list of saved ca
components ("saved_ca") that allow for quick recomputation of the CA
results. For random_direction_cutoff this saved_ca is empty.
