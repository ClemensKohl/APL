# Correspondance Analysis

\`cacomp\` performs correspondence analysis on a matrix or
Seurat/SingleCellExperiment object and returns the transformed data.

\`cacomp.seurat\` performs correspondence analysis on an assay from a
Seurat container and stores the standardized coordinates of the columns
(= cells) and the principal coordinates of the rows (= genes) as a
DimReduc Object in the Seurat container.

\`cacomp.SingleCellExperiment\` performs correspondence analysis on an
assay from a SingleCellExperiment and stores the standardized
coordinates of the columns (= cells) and the principal coordinates of
the rows (= genes) as a matrix in the SingleCellExperiment container.

## Usage

``` r
cacomp(
  obj,
  coords = TRUE,
  princ_coords = 3,
  python = FALSE,
  dims = NULL,
  top = 5000,
  inertia = TRUE,
  rm_zeros = TRUE,
  residuals = "pearson",
  cutoff = NULL,
  clip = FALSE,
  ...
)

# S4 method for class 'matrix'
cacomp(
  obj,
  coords = TRUE,
  princ_coords = 3,
  python = FALSE,
  dims = NULL,
  top = 5000,
  inertia = TRUE,
  rm_zeros = TRUE,
  residuals = "pearson",
  cutoff = NULL,
  clip = FALSE,
  ...
)

# S4 method for class 'dgCMatrix'
cacomp(
  obj,
  coords = TRUE,
  princ_coords = 3,
  python = FALSE,
  dims = NULL,
  top = 5000,
  inertia = TRUE,
  rm_zeros = TRUE,
  residuals = "pearson",
  cutoff = NULL,
  clip = FALSE,
  ...
)

# S4 method for class 'Seurat'
cacomp(
  obj,
  coords = TRUE,
  princ_coords = 3,
  python = FALSE,
  dims = NULL,
  top = 5000,
  inertia = TRUE,
  rm_zeros = TRUE,
  residuals = "pearson",
  cutoff = NULL,
  clip = FALSE,
  ...,
  assay = SeuratObject::DefaultAssay(obj),
  slot = "counts",
  return_input = FALSE
)

# S4 method for class 'SingleCellExperiment'
cacomp(
  obj,
  coords = TRUE,
  princ_coords = 3,
  python = FALSE,
  dims = NULL,
  top = 5000,
  inertia = TRUE,
  rm_zeros = TRUE,
  residuals = "pearson",
  cutoff = NULL,
  clip = FALSE,
  ...,
  assay = "counts",
  return_input = FALSE
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

  Other parameters

- assay:

  Character. The assay from which extract the count matrix for SVD, e.g.
  "RNA" for Seurat objects or "counts"/"logcounts" for
  SingleCellExperiments.

- slot:

  character. The slot of the Seurat assay. Default "counts".

- return_input:

  Logical. If TRUE returns the input (SingleCellExperiment/Seurat
  object) with the CA results saved in the reducedDim/DimReduc slot
  "CA". Otherwise returns a "cacomp". Default FALSE.

## Value

Returns a named list of class "cacomp" with components U, V and D: The
results from the SVD. row_masses and col_masses: Row and columns masses.
top_rows: How many of the most variable rows were retained for the
analysis. tot_inertia, row_inertia and col_inertia: Only if inertia =
TRUE. Total, row and column inertia respectively.

If return_imput = TRUE with Seurat container: Returns input obj of class
"Seurat" with a new Dimensional Reduction Object named "CA". Standard
coordinates of the cells are saved as embeddings, the principal
coordinates of the genes as loadings and the singular values (= square
root of principal intertias/eigenvalues) are stored as stdev. To
recompute a regular "cacomp" object without rerunning cacomp use
\`as.cacomp()\`.

If return_input =TRUE for SingleCellExperiment input returns a
SingleCellExperiment object with a matrix of standardized coordinates of
the columns in reducedDim(obj, "CA"). Additionally, the matrix contains
the following attributes: "prin_coords_rows": Principal coordinates of
the rows. "singval": Singular values. For the explained inertia of each
principal axis calculate singval^2. "percInertia": Percent explained
inertia of each principal axis. To recompute a regular "cacomp" object
from a SingleCellExperiment without rerunning cacomp use
\`as.cacomp()\`.

## Details

The calculation is performed according to the work of Michael Greenacre.
Singular value decomposition can be performed either with the base R
function \`svd\` or preferably by the faster pytorch implementation
(python = TRUE). When working with large matrices, CA coordinates and
principal coordinates should only be computed when needed to save
computational time.

## References

Greenacre, M. Correspondence Analysis in Practice, Third Edition, 2017.

## Examples

``` r
# Simulate scRNAseq data.
cnts <- data.frame(cell_1 = rpois(10, 5),
                   cell_2 = rpois(10, 10),
                   cell_3 = rpois(10, 20))
rownames(cnts) <- paste0("gene_", 1:10)
cnts <- as.matrix(cnts)

# Run correspondence analysis.
ca <- cacomp(obj = cnts, princ_coords = 3, top = 5)
#> No dimensions specified. Setting dimensions to: 2
#> Please consider setting the dimensions to a lower value to speed up the calculation.
#> Recommended dimensionality: << min(nrows, ncols) * 0.2

###########
# Seurat  #
###########
library(SeuratObject)
set.seed(1234)

# Simulate counts
cnts <- mapply(function(x){rpois(n = 500, lambda = x)},
                     x = sample(1:20, 50, replace = TRUE))
rownames(cnts) <- paste0("gene_", 1:nrow(cnts))
colnames(cnts) <- paste0("cell_", 1:ncol(cnts))

# Create Seurat object
seu <- CreateSeuratObject(counts = cnts)
#> Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
#> Warning: Data is of class matrix. Coercing to dgCMatrix.

# Run CA and save in dim. reduction slot
seu <- cacomp(seu, return_input = TRUE, assay = "RNA", slot = "counts")
#> Warning: 
#> Parameter top is >nrow(obj) and therefore ignored.
#> No dimensions specified. Setting dimensions to: 9

# Run CA and return cacomp object
ca <- cacomp(seu, return_input = FALSE, assay = "RNA", slot = "counts")
#> Warning: 
#> Parameter top is >nrow(obj) and therefore ignored.
#> No dimensions specified. Setting dimensions to: 9

########################
# SingleCellExperiment #
########################
library(SingleCellExperiment)
set.seed(1234)

# Simulate counts
cnts <- mapply(function(x){rpois(n = 500, lambda = x)},
               x = sample(1:20, 50, replace = TRUE))
rownames(cnts) <- paste0("gene_", 1:nrow(cnts))
colnames(cnts) <- paste0("cell_", 1:ncol(cnts))
logcnts <- log2(cnts + 1)

# Create SingleCellExperiment object
sce <- SingleCellExperiment(assays=list(counts=cnts, logcounts=logcnts))

# run CA and save in dim. reduction slot.
sce <- cacomp(sce, return_input = TRUE, assay = "counts") # on counts
#> No dimensions specified. Setting dimensions to: 9
sce <- cacomp(sce, return_input = TRUE, assay = "logcounts") # on logcounts
#> No dimensions specified. Setting dimensions to: 9

# run CA and return cacomp object.
ca <- cacomp(sce, return_input = FALSE, assay = "counts")
#> No dimensions specified. Setting dimensions to: 9
```
