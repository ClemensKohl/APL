# An S4 class that contains all elements needed for CA.

This class contains elements necessary to computer CA coordinates or
Association Plot coordinates, as well as other informative data such as
row/column inertia, gene-wise APL-scores, etc. ...

Creates new cacomp object.

## Usage

``` r
new_cacomp(...)
```

## Arguments

- ...:

  slot names and objects for new cacomp object.

## Value

cacomp object

## Slots

- `U`:

  class "matrix". Left singular vectors of the original input matrix.

- `V`:

  class "matrix". Right singular vectors of the original input matrix.

- `D`:

  class "numeric". Singular values of the original inpt matrix.

- `std_coords_rows`:

  class "matrix". Standardized CA coordinates of the rows.

- `std_coords_cols`:

  class "matrix". Standardized CA coordinates of the columns.

- `prin_coords_rows`:

  class "matrix". Principal CA coordinates of the rows.

- `prin_coords_cols`:

  class "matrix". Principal CA coordinates of the columns.

- `apl_rows`:

  class "matrix". Association Plot coordinates of the rows for the
  direction defined in slot "group"

- `apl_cols`:

  class "matrix". Association Plot coordinates of the columns for the
  direction defined in slot "group"

- `APL_score`:

  class "data.frame". Contains rows sorted by the APL score. Columns:
  Rowname (gene name in the case of gene expression data), APL score
  calculated for the direction defined in slot "group", the original row
  number and the rank of the row as determined by the score.

- `dims`:

  class "numeric". Number of dimensions in CA space.

- `group`:

  class "numeric". Indices of the chosen columns for APL calculations.

- `row_masses`:

  class "numeric". Row masses of the frequency table.

- `col_masses`:

  class "numeric". Column masses of the frequency table.

- `top_rows`:

  class "numeric". Number of most variable rows chosen.

- `tot_inertia`:

  class "numeric". Total inertia in CA space.

- `row_inertia`:

  class "numeric". Row-wise inertia in CA space.

- `col_inertia`:

  class "numeric". Column-wise inertia in CA space.

- `permuted_data`:

  class "list". Storage slot for permuted data.

- `params`:

  class "list". List of parameters.

## Examples

``` r
set.seed(1234)

# Simulate counts
cnts <- mapply(function(x){rpois(n = 500, lambda = x)}, 
               x = sample(1:20, 50, replace = TRUE))
rownames(cnts) <- paste0("gene_", 1:nrow(cnts))
colnames(cnts) <- paste0("cell_", 1:ncol(cnts))

res <-  APL:::comp_std_residuals(mat=cnts)
SVD <- svd(res$S)
names(SVD) <- c("D", "U", "V")
SVD <- SVD[c(2, 1, 3)]

ca <- new_cacomp(U = SVD$U,
                 V = SVD$V,
                 D = SVD$D,
                 row_masses = res$rowm,
                 col_masses = res$colm)
```
