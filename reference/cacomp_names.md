# Prints slot names of cacomp object

Prints slot names of cacomp object

## Usage

``` r
cacomp_names(caobj)
```

## Arguments

- caobj:

  a cacomp object

## Value

Prints slot names of cacomp object

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

# show slot names:
cacomp_names(ca)
#>  [1] "U"                "V"                "D"                "std_coords_rows" 
#>  [5] "std_coords_cols"  "prin_coords_rows" "prin_coords_cols" "apl_rows"        
#>  [9] "apl_cols"         "APL_score"        "params"           "dims"            
#> [13] "group"            "row_masses"       "col_masses"       "top_rows"        
#> [17] "tot_inertia"      "row_inertia"      "col_inertia"      "permuted_data"   
```
