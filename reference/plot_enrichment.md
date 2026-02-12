# Generates plot for results from apl_topGO

Plots the results from the data frame generated via apl_topGO.

## Usage

``` r
plot_enrichment(genenr, ntop = 10)
```

## Arguments

- genenr:

  data.frame. gene enrichment results table.

- ntop:

  numeric. Number of elements to plot.

## Value

Returns a ggplot plot.

## Examples

``` r
library(SeuratObject)
set.seed(1234)
cnts <- SeuratObject::LayerData(pbmc_small, assay = "RNA", layer = "counts")
cnts <- as.matrix(cnts)

# Run CA on example from Seurat

ca <- cacomp(pbmc_small,
             princ_coords = 3,
             return_input = FALSE,
             assay = "RNA",
             slot = "counts")
#> Warning: 
#> Parameter top is >nrow(obj) and therefore ignored.
#> No dimensions specified. Setting dimensions to: 15

grp <- which(Idents(pbmc_small) == 2)
ca <- apl_coords(ca, group = grp)
ca <- apl_score(ca,
                mat = cnts)
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======                                                               |  10%  |                                                                              |==============                                                        |  20%  |                                                                              |=====================                                                 |  30%  |                                                                              |============================                                          |  40%  |                                                                              |===================================                                   |  50%  |                                                                              |==========================================                            |  60%  |                                                                              |=================================================                     |  70%  |                                                                              |========================================================              |  80%  |                                                                              |===============================================================       |  90%  |                                                                              |======================================================================| 100%

enr <- apl_topGO(ca,
                 ontology = "BP",
                 organism = "hs")
#> 
#> groupGOTerms:    GOBPTerm, GOMFTerm, GOCCTerm environments built.
#> 
#> Building most specific GOs .....
#>  ( 1453 GO terms found. )
#> 
#> Build GO DAG topology ..........
#>  ( 3595 GO terms and 7787 relations. )
#> 
#> Annotating nodes ...............
#>  ( 205 genes annotated to the GO terms. )
#> 
#>           -- Elim Algorithm -- 
#> 
#>       the algorithm is scoring 492 nontrivial nodes
#>       parameters: 
#>           test statistic: fisher
#>           cutOff: 0.01
#> 
#>   Level 12:  2 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 11:  6 nodes to be scored    (0 eliminated genes)
#> 
#>   Level 10:  17 nodes to be scored   (8 eliminated genes)
#> 
#>   Level 9:   27 nodes to be scored   (11 eliminated genes)
#> 
#>   Level 8:   57 nodes to be scored   (11 eliminated genes)
#> 
#>   Level 7:   77 nodes to be scored   (11 eliminated genes)
#> 
#>   Level 6:   98 nodes to be scored   (23 eliminated genes)
#> 
#>   Level 5:   90 nodes to be scored   (23 eliminated genes)
#> 
#>   Level 4:   67 nodes to be scored   (23 eliminated genes)
#> 
#>   Level 3:   37 nodes to be scored   (23 eliminated genes)
#> 
#>   Level 2:   13 nodes to be scored   (23 eliminated genes)
#> 
#>   Level 1:   1 nodes to be scored    (23 eliminated genes)

plot_enrichment(enr)
```
