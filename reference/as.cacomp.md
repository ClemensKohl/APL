# Create cacomp object from Seurat/SingleCellExperiment container

Converts the values stored in the Seurat/SingleCellExperiment
dimensional reduction slot "CA" to a cacomp object. If recompute = TRUE
additional parameters are recomputed from the saved values without
rerunning SVD (need to specify assay to work).

as.cacomp.cacomp returns input without any calculations.

Recomputes missing values and returns cacomp object from a list. If you
have a \*complete\* cacomp object in list form, use do.call(new_cacomp,
obj).

as.cacomp.Seurat: Converts the values stored in the Seurat DimReduc slot
"CA" to an cacomp object.

as.cacomp.SingleCellExperiment: Converts the values stored in the
SingleCellExperiment reducedDim slot "CA" to a cacomp object.

## Usage

``` r
as.cacomp(obj, ...)

# S4 method for class 'cacomp'
as.cacomp(obj, ...)

# S4 method for class 'list'
as.cacomp(obj, ..., mat = NULL)

# S4 method for class 'Seurat'
as.cacomp(obj, ..., assay = "RNA", slot = "counts")

# S4 method for class 'SingleCellExperiment'
as.cacomp(obj, ..., assay = "counts")
```

## Arguments

- obj:

  An object of class "Seurat" or "SingleCellExperiment" with a dim.
  reduction named "CA" saved. For obj "cacomp" input is returned.

- ...:

  Further arguments.

- mat:

  Original input matrix.

- assay:

  Character. The assay from which extract the count matrix, e.g. "RNA"
  for Seurat objects or "counts"/"logcounts" for SingleCellExperiments.

- slot:

  character. Slot of the Seurat assay to use. Default "counts".

## Value

A cacomp object.

## Details

By default extracts std_coords_cols, D, prin_coords_rows, top_rows and
dims from obj and outputs a cacomp object. If recompute = TRUE the
following are additionally recalculated (doesn't run SVD): U, V,
std_coords_rows, row_masses, col_masses.

## Examples

``` r
#########
# lists #
#########

# Simulate counts
cnts <- mapply(function(x){rpois(n = 500, lambda = x)},
               x = sample(1:100, 50, replace = TRUE))
rownames(cnts) <- paste0("gene_", 1:nrow(cnts))
colnames(cnts) <- paste0("cell_", 1:ncol(cnts))

# Run correspondence analysis
ca <- cacomp(obj = cnts, princ_coords = 3)
#> Warning: 
#> Parameter top is >nrow(obj) and therefore ignored.
#> No dimensions specified. Setting dimensions to: 9
ca_list <- as.list(ca)

# Only keep subset of elements for demonstration
ca_list <- ca_list[c("U", "std_coords_rows", "std_coords_cols", "params")]

# convert (incomplete) list to cacomp object.
ca <- as.cacomp(ca_list, mat = cnts)
#> Calling cacomp to recompute from matrix.
#> Please consider setting the dimensions to a lower value to speed up the calculation.
#> Recommended dimensionality: << min(nrows, ncols) * 0.2

##########
# Seurat #
##########
library(SeuratObject)
set.seed(1234)

# Simulate counts
cnts <- mapply(function(x){rpois(n = 500, lambda = x)},
               x = sample(1:100, 50, replace = TRUE))
rownames(cnts) <- paste0("gene_", 1:nrow(cnts))
colnames(cnts) <- paste0("cell_", 1:ncol(cnts))

seu <- CreateSeuratObject(counts = cnts)
#> Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
#> Warning: Data is of class matrix. Coercing to dgCMatrix.
seu <- cacomp(seu, return_input = TRUE)
#> Warning: 
#> Parameter top is >nrow(obj) and therefore ignored.
#> No dimensions specified. Setting dimensions to: 9

ca <- as.cacomp(seu, assay = "RNA", slot = "counts")

########################
# SingleCellExperiment #
########################
library(SingleCellExperiment)
#> Loading required package: SummarizedExperiment
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: ‘matrixStats’
#> The following objects are masked from ‘package:Biobase’:
#> 
#>     anyMissing, rowMedians
#> 
#> Attaching package: ‘MatrixGenerics’
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> The following object is masked from ‘package:Biobase’:
#> 
#>     rowMedians
#> Loading required package: GenomicRanges
#> Loading required package: Seqinfo
#> 
#> Attaching package: ‘SummarizedExperiment’
#> The following object is masked from ‘package:SeuratObject’:
#> 
#>     Assays
set.seed(1234)

# Simulate counts
cnts <- mapply(function(x){rpois(n = 500, lambda = x)},
               x = sample(1:100, 50, replace = TRUE))
rownames(cnts) <- paste0("gene_", 1:nrow(cnts))
colnames(cnts) <- paste0("cell_", 1:ncol(cnts))

sce <- SingleCellExperiment(assays=list(counts=cnts))
sce <- cacomp(sce, return_input = TRUE)
#> No dimensions specified. Setting dimensions to: 9

ca <- as.cacomp(sce, assay = "counts")
```
