

#' Compute Standard Residuals
#'
#' @description
#' `comp_std_residuals` computes the standard Residuals matrix S, which isthe basis for correspondence analysis and serves as input for singular value decomposition.
#'
#' @details
#' Calculates standardized residual matrix S from the proportions matrix P and
#' the expected values E according to \eqn{S = \frac{(P-E)}{sqrt(E)}}.
#'
#' @param mat A numerical matrix or coercible to one by `as.matrix()` Should have row and column names.
#' @return
#' A named list with standard residual matrix "S", grand total of the original matrix "tot"
#'  as well as row and column masses "rowm" and "colm" respectively.
#'
comp_std_residuals <- function(mat){


  if (!is(mat, "matrix")){
    mat <- as.matrix(mat)
  }
  stopifnot("Input matrix does not have any rownames!" = !is.null(rownames(mat)))
  stopifnot("Input matrix does not have any colnames!" = !is.null(colnames(mat)))

  tot <- sum(mat)
  P <- mat/tot               # proportions matrix
  rowm <- rowSums(P)          # row masses
  colm <- colSums(P)          # column masses

  E <- rowm %o% colm      # expected proportions
  S <-  (P - E) / sqrt(E)         # standardized residuals
  S[is.nan(S)] <- 0

  out <- list("S"=S, "tot"=tot, "rowm"=rowm, "colm"=colm)
  return(out)
}

#' Find most variable rows
#'
#' @description
#' Calculates the variance of the chi-square components matrix and selects the top e.g 5000 rows.
#'
#' @return
#' Returns a matrix, which consists of the top variable rows of mat.
#'
#' @param mat A numeric matrix. For sequencing a count matrix, gene expression values with genes in rows and samples/cells in columns.
#' Should contain row and column names.
#' @param top Integer. Number of most variable rows to retain. Default 5000.
#'
var_rows <- function(mat, top = 5000){

  res <-  comp_std_residuals(mat=mat)

  if(top>nrow(mat)) warning("Top is larger than the number of rows in matrix. Top was set to nrow(mat).")
  top <- min(nrow(mat), top)

  chisquare <- res$tot * (res$S^2)		# chi-square components matrix
  variances <- apply(chisquare,1,var) #row-wise variances
  ix_var <- order(-variances)
  mat <- mat[ix_var[1:top],] # choose top rows
  return(mat)

}

#' Correspondance Analysis
#'
#' @description
#' `cacomp` performs correspondence analysis on a matrix or Seurat/SingleCellExperiment object and returns the transformed data.
#'
#' @details
#' The calculation is performed according to Greenacre. Singular value decomposition
#' can be performed either with the base R function `svd` or preferably by the much faster
#' pytorch implementation (python = TRUE). When working on large matrices, CA coordinates and
#' principal coordinates should only computed when needed to save computational time.
#'
#' @return
#' Returns a named list of class "cacomp" with components
#' U, V and D: The results from the SVD.
#' row_masses and col_masses: Row and columns masses.
#' top_rows: How many of the most variable rows were retained for the analysis.
#' tot_inertia, row_inertia and col_inertia: Only if inertia = TRUE. Total, row and column inertia respectively.
#' @references
#' Greenacre, M. Correspondence Analysis in Practice, Third Edition, 2017.

#' @param obj A numeric matrix or Seurat/SingleCellExperiment object. For sequencing a count matrix, gene expression values with genes in rows and samples/cells in columns.
#' Should contain row and column names.
#' @param coords Logical. Indicates whether CA standard coordinates should be calculated. Default TRUE
#' @param python A logical value indicating whether to use singular-value decomposition from the python package torch.
#' This implementation dramatically speeds up computation compared to `svd()` in R.
#' @param princ_coords Integer. Number indicating whether principal coordinates should be calculated for the rows (=1), columns (=2), both (=3) or none (=0).
#' Default 1.
#' @param dims Integer. Number of CA dimensions to retain. Default NULL (keeps all dimensions).
#' @param top Integer. Number of most variable rows to retain. Default NULL.
#' @param inertia Logical.. Whether total, row and column inertias should be calculated and returned. Default TRUE.
#' @param rm_zeros Logical. Whether rows & cols containing only 0s should be removed. Keeping zero only rows/cols might lead to unexpected results. Default TRUE.
#' @param ... Arguments forwarded to methods.
#' @export
cacomp <- function(obj, coords=TRUE, princ_coords = 1, python = TRUE, dims = NULL, top = NULL, inertia = TRUE, rm_zeros = TRUE, ...){
  UseMethod("cacomp")
}

#' @rdname cacomp
#' @export
cacomp.default <- function(obj, coords=TRUE, princ_coords = 1, python = TRUE, dims = NULL, top = NULL, inertia = TRUE, rm_zeros = TRUE, ...){
  stop(paste0("cacomp does not know how to handle objects of class ",
              class(obj),
              ". Currently only objects of class 'matrix' or objects coercible to one, 'Seurat' or 'SingleCellExperiment' are supported."))
}


#' @rdname cacomp
#' @export
cacomp.matrix <- function(obj, coords=TRUE, princ_coords = 1, python = TRUE, dims = NULL, top = NULL, inertia = TRUE, rm_zeros = TRUE, ...){

  # chkDots(...)

  stopifnot("Input matrix does not have any rownames!" = !is.null(rownames(obj)))
  stopifnot("Input matrix does not have any colnames!" = !is.null(colnames(obj)))

  if (rm_zeros == TRUE){
    no_zeros_rows <- rowSums(obj) > 0
    no_zeros_cols <- colSums(obj) > 0
    if (sum(!no_zeros_rows) != 0){
      ## Delete genes with only with only zero values across all conditions
      warning("Matrix contains rows with only 0s. These rows were removed. If undesired set rm_zeros = FALSE.")
      obj <- obj[no_zeros_cols,]
    }
    if (sum(!no_zeros_cols) != 0){
      ## Delete genes with only with only zero values across all conditions
      warning("Matrix contains columns with only 0s. These columns were removed. If undesired set rm_zeros = FALSE.")
      obj <- obj[,no_zeros_cols]
    }
  }


  # Choose only top # of variable  genes
  if (is.null(top) || top == nrow(obj)) {
    res <-  comp_std_residuals(mat=obj)
    toptmp <- nrow(obj)
  } else if (!is.null(top) && top < nrow(obj)){
    # message("Subsetting matrix by the top ", top," most variable rows (chisquare variances) ...")
    obj <- var_rows(mat = obj, top = top)
    res <-  comp_std_residuals(mat=obj)
    toptmp <- top
  } else if (top > nrow(obj)) {
    warning("Parameter top is >nrow(obj) and therefore ignored.")
    res <-  comp_std_residuals(mat=obj)
    toptmp <- nrow(obj)
  } else {
    warning("Unusual input for top, argument ignored.")
    res <-  comp_std_residuals(mat=obj)
    toptmp <- nrow(obj)
  }

  S <- res$S
  tot <- res$tot
  rowm <- res$rowm
  colm <- res$colm
  rm(res)

  n <- nrow(S)
  p <- ncol(S)
  k <- min(n,p)

  # S <- (diag(1/sqrt(r)))%*%(P-r%*%t(c))%*%(diag(1/sqrt(c)))
  message("Running singular value decomposition ...")

  if (python == TRUE){
    # require(reticulate)
    # source_python('./python_svd.py')
    reticulate::source_python(system.file("python/python_svd.py", package = "APL"))
    SVD <- svd_torch(S)
    names(SVD) <- c("U", "D", "V")
    SVD$D <- as.vector(SVD$D)

  } else {
    SVD <- svd(S, nu = k, nv = k)
    names(SVD) <- c("D", "U", "V")
    SVD <- SVD[c(2, 1, 3)]
  }

  j <- seq_len(k)
  dimnames(SVD$V) <- list(colnames(S), paste0("Dim", j))
  dimnames(SVD$U) <- list(rownames(S), paste0("Dim", j))


  if(inertia == TRUE){
    #calculate inertia
    SVD$tot_inertia <- sum(SVD$D^2)
    SVD$row_inertia <- rowSums(S^2)
    SVD$col_inertia <- colSums(S^2)
  }

  SVD$row_masses <- rowm
  SVD$col_masses <- colm
  SVD$top_rows <- toptmp

  SVD <- new_cacomp(SVD)
  # class(SVD) <- "cacomp"

  if (coords == TRUE){
    message("Calculating coordinates ...")

    SVD <- ca_coords(caobj = SVD,
                    dims = dims,
                    princ_coords = princ_coords,
                    princ_only = FALSE)
  } else {
    if(!is.null(dims)){
      if (dims >= length(SVD$D)){
        if (dims > length(SVD$D)){
          warning("Chosen dimensions are larger than the number of dimensions obtained from the singular value decomposition. Argument ignored.")
        }
        SVD$dims <- length(SVD$D)
      } else {
        dims <- min(dims, length(SVD$D))
        SVD$dims <- dims

        dims <- seq(dims)

        # subset to number of dimensions
        SVD$U <- SVD$U[,dims]
        SVD$V <- SVD$V[,dims]
        SVD$D <- SVD$D[dims]
      }
    } else {
      SVD$dims <- length(SVD$D)
    }
  }


  return(SVD)

}


#' Correspondance Analysis for Seurat objects
#'
#' @description
#' `cacomp.seurat` performs correspondence analysis on a assay from a Seurat container and stores the standard coordinates
#'  of the columns (= cells) and the principal coordinates of the rows (= genes) as a DimReduc Object in the Seurat container.
#'
#' @return
#' If return_imput = TRUE with Seurat container: Returns input obj of class "Seurat" with a new Dimensional Reduction Object named "CA".
#' Standard coordinates of the cells are saved as embeddings,
#' the principal coordinates of the genes as loadings and
#' the singular values (= square root of principal intertias/eigenvalues)
#' are stored as stdev.
#' To recompute a regular "cacomp" object without rerunning cacomp use `as.cacomp()`.
#' @param assay Character. The assay from which extract the count matrix for SVD, e.g. "RNA" for Seurat objects or "counts"/"logcounts" for SingleCellExperiments.
#' @param return_input Logical. If TRUE returns the input (SingleCellExperiment/Seurat object) with the CA results saved in the reducedDim/DimReduc slot "CA".
#'  Otherwise returns a "cacomp". Default FALSE.
#' @param ... Other parameters
#' @rdname cacomp
#' @export
cacomp.Seurat <- function(obj, coords=TRUE, princ_coords = 1, python = TRUE, dims = NULL, top = NULL, inertia = TRUE, rm_zeros = TRUE, ...,assay = DefaultAssay(obj), return_input = FALSE){

  stopifnot("obj doesnt belong to class 'Seurat'" = is(obj, "Seurat"))

  stopifnot("Set coords = TRUE when inputting a Seurat object and return_input = TRUE." = coords == TRUE)


  seu <- Seurat::GetAssayData(object = obj, assay = assay, slot = "data")
  seu <- as.matrix(seu)

  caobj <- cacomp(obj = seu,
                  coords = coords,
                  top = top,
                  princ_coords = princ_coords,
                  dims = dims,
                  python = python,
                  rm_zeros = rm_zeros,
                  inertia = inertia)

  if (return_input == TRUE){
    colnames(caobj$V) <- paste0("DIM_", seq(ncol(caobj$V)))
    colnames(caobj$U) <- paste0("DIM_", seq(ncol(caobj$U)))

    obj[["CA"]] <- Seurat::CreateDimReducObject(embeddings = caobj$std_coords_cols,
                                               loadings = caobj$prin_coords_rows,
                                               stdev = caobj$D,
                                               key = "DIM_",
                                               assay = assay)

    return(obj)
  } else {
    return(caobj)
  }

}



#' @description
#' `cacomp.SingleCellExperiment` performs correspondence analysis on an assay from a SingleCellExperiment and stores the standard coordinates
#'  of the columns (= cells) and the principal coordinates of the rows (= genes) as a matrix in the SingleCellExperiment container.
#'
#' @return
#' If return_input =TRUE for SingleCellExperiment input returns a SingleCellExperiment object with a matrix of standard coordinates of the columns in
#' reducedDim(obj, "CA"). Additionally, the matrix contains the following attributes:
#' "prin_coords_rows": Principal coordinates of the rows.
#' "singval": Singular values. For the explained inertia of each principal axis calculate singval^2.
#' "percInertia": Percent explained inertia of each principal axis.
#' To recompute a regular "cacomp" object from a SingleCellExperiment without rerunning cacomp use `as.cacomp()`.
#' @param assay Character. The assay from which extract the count matrix for SVD, e.g. "RNA" for Seurat objects or "counts"/"logcounts" for SingleCellExperiments.
#' @param return_input Logical. If TRUE returns the input (SingleCellExperiment/Seurat object) with the CA results saved in the reducedDim/DimReduc slot "CA".
#'  Otherwise returns a "cacomp". Default FALSE.
#' @rdname cacomp
#' @export
cacomp.SingleCellExperiment <- function(obj, coords=TRUE, princ_coords = 1, python = TRUE, dims = NULL, top = NULL, inertia = TRUE, rm_zeros = TRUE, ..., assay = "counts", return_input = FALSE){

  stopifnot("obj doesnt belong to class 'SingleCellExperiment'" = is(obj, "SingleCellExperiment"))
  stopifnot("Set coords = TRUE when inputting a SingleCellExperiment object and return_input = TRUE." = coords == TRUE)

  mat <- SummarizedExperiment::assay(obj, assay)
  mat <- as.matrix(mat)

  top <- min(nrow(mat), top)

  caobj <- cacomp.matrix(obj = mat,
                         coords = coords,
                         top = top,
                         princ_coords = princ_coords,
                         dims = dims,
                         python = python,
                         rm_zeros = rm_zeros,
                         inertia = inertia)

  if (return_input == TRUE){
    prinInertia <- caobj$D^2
    percentInertia <- prinInertia / sum(prinInertia) * 100 #TODO IS THIS CORRECT ??

    # Saving the results
    ca <- caobj$std_coords_cols
    attr(ca, "prin_coords_rows") <- caobj$prin_coords_rows
    attr(ca, "singval") <- caobj$D
    attr(ca, "percInertia") <- percentInertia

    SingleCellExperiment::reducedDim(obj, "CA") <- ca

    return(obj)

  } else {
    return(caobj)
  }

  # old code
  # lme <- LinearEmbeddingMatrix(sampleFactors = SVD$std_coords_cols,
  #                             featureLoadings =SVD$prin_coords_rows,
  #                             factorData = DataFrame("D" = SVD$D, row.names = paste0("Dim", seq(length(SVD$D)))),
  #                             metadata = list())


}


#' Subset dimensions of a caobj
#'
#' @description Subsets the dimensions according to user input.
#'
#' @return Returns caobj.
#'
#' @param caobj A caobj.
#' @param dims Integer. Number of dimensions.
subset_dims <- function(caobj, dims){

  stopifnot(is(caobj, "cacomp"))

  if(dims > length(caobj$D)){
    warning("dims is larger than the number of available dimensions. Argument ignored")
  } else if (dims == length(caobj$D)){
    caobj$dims <- dims
    return(caobj)
  }

  dims <- min(dims, length(caobj$D))
  caobj$dims <- dims
  dims <- seq(dims)
  caobj$U <- caobj$U[,dims]
  caobj$V <- caobj$V[,dims]
  caobj$D <- caobj$D[dims]

  if (is.null(caobj$std_coords_cols)){
    caobj$std_coords_cols <- caobj$std_coords_cols[,dims]
    if (is.null(caobj$prin_coords_cols)){
      caobj$prin_coords_cols <- caobj$prin_coords_cols[,dims]

    }
  }

  if (is.null(caobj$std_coords_rows)){
    caobj$std_coords_rows <- caobj$std_coords_rows[,dims]
    if (is.null(caobj$prin_coords_rows)){
      caobj$prin_coords_rows <- caobj$prin_coords_rows[,dims]
    }
  }

  return(caobj)
}


#' Calculate correspondence analysis row and column coordinates.
#'
#' @description `ca_coords` calculates the standard and principal coordinates of the rows and columns in CA space.
#'
#' @details
#' Takes a "cacomp" object and calculates standard and principal coordinates for the visualisation of CA results in a biplot or
#' to subsequently calculate coordinates in an association plot.
#'
#' @return
#' Returns input object with coordinates added.
#' std_coords_rows/std_coords_cols: Standard coordinates of rows/columns.
#' prin_coords_rows/prin_coords_cols: Principal coordinates of rows/columns.
#'
#' @param caobj A "cacomp" object as outputted from `cacomp()`.
#' @param dims Integer indicating the number of dimensions to use for the calculation of coordinates.
#' All elements of caobj (where applicable) will be reduced to the given number of dimensions. Default NULL (keeps all dimensions).
#' @param princ_only Logical, whether only principal coordinates should be calculated.
#' Or, in other words, whether the standard coordinates are already calculated and stored in `caobj`. Default `FALSE`.
#' @param princ_coords Integer. Number indicating whether principal coordinates should be calculated for the rows (=1), columns (=2), both (=3) or none (=0).
#' Default 3.
#' @export
ca_coords <- function(caobj, dims=NULL, princ_coords = 3, princ_only = FALSE){

  stopifnot(is(caobj, "cacomp"))
  stopifnot(dims <= length(caobj$D))

  if(!is.null(dims)){
    if (dims > length(caobj$D)){
      warning("Chosen dimensions are larger than the number of dimensions obtained from the singular value decomposition. Argument ignored.")
     } #else {
    #   if ("dims" %in% names(caobj)) {
    #     warning("The caobj was previously already subsetted to ", caobj$dims, " dimensions. Subsetting again!")
    #   }
      caobj <- subset_dims(caobj = caobj, dims = dims)
      # dims <- min(dims, length(caobj$D))
      # caobj$dims <- dims
      # dims <- seq(dims)
      # # subset to number of dimensions
      # caobj$U <- caobj$U[,dims]
      # caobj$V <- caobj$V[,dims]
      # caobj$D <- caobj$D[dims]
      #
      # if (princ_only == TRUE){
      #   stopifnot(!is.null(caobj$std_coords_rows))
      #   stopifnot(!is.null(caobj$std_coords_cols))
      #
      #   caobj$std_coords_rows <- caobj$std_coords_rows[,dims]
      #   caobj$std_coords_cols <- caobj$std_coords_cols[,dims]
      # }
    }




  if(princ_only == FALSE){

    #standard coordinates
    caobj$std_coords_rows <- sweep(caobj$U, 1, sqrt(caobj$row_masses), "/")
    caobj$std_coords_cols <- sweep(caobj$V, 1, sqrt(caobj$col_masses), "/")

    # Ensure no NA/Inf after dividing by 0.
    caobj$std_coords_rows[is.na(caobj$std_coords_rows)] <- 0
    caobj$std_coords_cols[is.na(caobj$std_coords_cols)] <- 0
    caobj$std_coords_rows[is.infinite(caobj$std_coords_rows)] <- 0
    caobj$std_coords_cols[is.infinite(caobj$std_coords_cols)] <- 0

  }


  stopifnot("princ_coords must be either 0, 1, 2 or 3" = (princ_coords == 0 || princ_coords == 1 || princ_coords == 2 || princ_coords == 3))

  if(princ_coords != 0){
    stopifnot(!is.null(caobj$std_coords_rows))
    stopifnot(!is.null(caobj$std_coords_cols))

      if (princ_coords == 1){
        #principal coordinates for rows
        caobj$prin_coords_rows <- sweep(caobj$std_coords_rows, 2, caobj$D, "*")
      } else if (princ_coords == 2) {
        #principal coordinates for columns
        caobj$prin_coords_cols <- sweep(caobj$std_coords_cols, 2, caobj$D, "*")
      } else if (princ_coords  == 3) {
        #principal coordinates for rows
        caobj$prin_coords_rows <- sweep(caobj$std_coords_rows, 2, caobj$D, "*")
        #principal coordinates for columns
        caobj$prin_coords_cols <- sweep(caobj$std_coords_cols, 2, caobj$D, "*")
      }

  }

  return(caobj)
}


#' Scree Plot
#'
#'@description Plots a scree plot.
#'
#'@return
#'Returns a ggplot object.
#'
#'@param df A data frame with columns "dims" and "inertia".
scree_plot <- function(df){

  stopifnot(c("dims", "inertia") %in% colnames(df))

  avg_inertia <- 100/nrow(df)
  max_num_dims <- nrow(df)

  screeplot <- ggplot2::ggplot(df, ggplot2::aes(x=dims, y=inertia)) +
    ggplot2::geom_col(fill="#4169E1") +
    # geom_point(color="#B22222")+
    ggplot2::geom_line(color="#B22222", size=1) +
    ggplot2::geom_abline(slope = 0, intercept = avg_inertia, linetype=3, alpha=0.8, color ="#606060")+
    ggplot2::labs(title = "Scree plot of explained inertia per dimensions and the average inertia",
         y="Explained inertia [%]",
         x="Dimension") +
    ggplot2::annotate(
      "text",
      x = max_num_dims*0.9,
      y = avg_inertia,
      label = "avg. inertia")+
    # scale_color_identity(name = "Explained Inertia of:",
    #                      breaks = c("#B22222"),
    #                      labels = c("Data"),
    #                      guide = "legend")+
    ggplot2::theme_bw()
  return(screeplot)
}


#' Compute statistics to help choose the number of dimensions
#'
#' @description
#' Allow the user to choose from 4 different methods ("avg_inertia", "maj_inertia", "scree_plot" and "elbow_rule")
#' to estimate the number of dimensions that best represent the data.
#'
#' @details
#' "avg_inertia" calculates the number of dimensions in which the inertia is above the average inertia.
#' "maj_inertia" calculates the number of dimensions in which cumulatively explain up to 80\% of the total inertia.
#' "scree_plot" plots a scree plot.
#' "elbow_rule" Formalization of the commonly used elbow rule. Permutes the rows for each column and reruns `cacomp()` for a total of `reps` times.
#'  The number of relevant dimensions is obtained from the point where the line for the explained inertia of the permuted data intersects with the actual data.
#'
#' @return
#' For `avg_inertia`, `maj_inertia` and `elbow_rule` (when `return_plot=FALSE`) returns an integer, indicating the suggested number of dimensions to use.
#' `scree_plot` returns a ggplot object.
#' `elbow_rule` (for `return_plot=TRUE`) returns a list with two elements: "dims" contains the number of dimensions and "plot" a ggplot.
#'
#' @param obj A "cacomp" object as outputted from `cacomp()`,
#' a "Seurat" object with a "CA" DimReduc object stored,
#' or a "SingleCellExperiment" object with a "CA" dim. reduction stored.
#' @param mat A numeric matrix. For sequencing a count matrix, gene expression values with genes in rows and samples/cells in columns.
#' Should contain row and column names.
#' @param method String. Either "scree_plot", "avg_inertia", "maj_inertia" or "elbow_rule" (see Details section). Default "scree_plot".
#' @param reps Integer. Number of permutations to perform when choosing "elbow_rule". Default 3.
#' @param return_plot TRUE/FALSE. Whether a plot should be returned when choosing "elbow_rule". Default FALSE.
#' @param python A logical value indicating whether to use singular-value decomposition from the python package torch.
#' This implementation dramatically speeds up computation compared to `svd()` in R.
#' @param ... Arguments forwarded to methods.
#' @export
pick_dims <- function(obj, mat = NULL, method="scree_plot", reps=3, python = TRUE, return_plot = FALSE, ...){
  UseMethod("pick_dims")
}


#' @rdname pick_dims
#' @export
pick_dims.default <- function(obj, mat = NULL, method="scree_plot", reps=2, python = TRUE, return_plot = FALSE, ...){
  stop(paste0("pick_dims does not know how to handle objects of class ",
              class(obj),
              ". Currently only objects of class 'cacomp', 'Seurat' or 'SingleCellExperiment' are supported."))
}




#' @rdname pick_dims
#' @export
pick_dims.cacomp <- function(obj, mat = NULL, method="scree_plot", reps=3, python = TRUE, return_plot = FALSE, ...){

  if (!is(obj,"cacomp")){
    stop("Not a CA object. Please run cacomp() first!")
  }

  ev <- obj$D^2
  expl_inertia <- (ev/sum(ev)) *100
  max_num_dims <- length(obj$D)

  if (method == "avg_inertia"){
    # Method 1: Dim's > average inertia
    avg_inertia <- 100/max_num_dims		# percentage of inertia explained by 1 dimension (on average)
    dim_num <- sum(expl_inertia > avg_inertia)  # result: number of dimensions, all of which explain more than avg_inertia
    return(dim_num)

  } else if (method == "maj_inertia"){
    # Method 2: Sum of dim's > 80% of the total inertia
    dim_num <- min(which(cumsum(expl_inertia)>80)) # the first dimension for which the cumulative sum of inertias (from dim1 up to given dimension) is higher than 80%
    return(dim_num)

  } else if (method == "scree_plot"){
    # Method 3: Graphical representation of explained inertias (scree plot)
    # the user can set the threshold based on the scree plot

    df <- data.frame(dims = 1:max_num_dims,
                     inertia = expl_inertia)

    screeplot <- scree_plot(df)

    return(screeplot)

  } else if (method == "elbow_rule") {

    if(is.null(mat)){
      cat("When running method=\"elbow_rule\", please provide the original data matrix (paramater mat) which was earlier submitted to cacomp()!")
      stop()
    }
    # Method 4: Formalization of the elbow rule

    # 1.Generate an artificial data matrix by randomly permuting the responses (rows)
    #   of each sample (column)  of the original expression matrix
    # 2.Perform the CA calculations (I) for the generated data
    # 3.Repeat the steps 1. to 3. 5 times and save the calculated
    #   explained inertia vector (expl_inertia_perm) for each artificial matrix
    #   as a row in an overall inertia matrix (matrix_expl_inertia_perm)

    matrix_expl_inertia_perm <- matrix(0, nrow = max_num_dims , ncol = reps)

    for (k in seq(reps)) {
      message("Running permutation ", k, " out of ", reps, " for elbow rule ...")
      mat <- as.matrix(mat)
      mat_perm <- apply(mat, 2, FUN=sample)
      colnames(mat_perm) <- colnames(mat)
      rownames(mat_perm) <- 1:nrow(mat_perm)

      obj_perm <- cacomp(obj=mat_perm, top = obj$top_rows, dims = obj$dims, coords = FALSE, python = python)

      ev_perm <- obj_perm$D^2
      expl_inertia_perm <- (ev_perm/sum(ev_perm))*100

      matrix_expl_inertia_perm[,k] <- expl_inertia_perm
      colnames(matrix_expl_inertia_perm) <- paste0("perm",1:reps)
    }

    if (return_plot == TRUE){
      df <- data.frame(dims = 1:max_num_dims,
                       inertia = expl_inertia)

      df <- cbind(df, matrix_expl_inertia_perm)

      screeplot <- scree_plot(df)

      for (k in 1:reps) {

        colnm <- paste0("perm",k)
        screeplot <- screeplot +
          ggplot2::geom_line(data = df, ggplot2::aes(x=dims, y=.data[[colnm]]), color="black", alpha=0.8, linetype=2)

      }
    }
    # 5.Identify the largest number of dimension based on intersection of the real data's scree plot
    #   with the average of simulated scree plots

    avg_inertia_perm <- rowMeans(matrix_expl_inertia_perm)		# average artificial inertia vector

    tmp <- as.integer(expl_inertia>avg_inertia_perm)						# check for each position (dimension) if expl_inertia is higher than the avg_artificial_inertia. If yes =1, if no =0
    if (sum(tmp)==0 || sum(tmp)==max_num_dims){
      dim_number <- max_num_dims									# result: if the lines do not intersect, choose max_num_dims
    } else if (tmp[1] == 0){
      stop("Average inertia of the permutated data is above the explained inertia of the data in the first dimension. Please either try more permutations or a different method.")
    }else{
      dim_number <- length(tmp[cumsum(tmp == 0)<1 & tmp!=0])		# result: if the lines intersect at at 1 or more than 1 positions, choose the intersection with the lowest (BEFORE: it was with highest) x-coordinate (:= the first intersection!)
    }
    #TODO if permutations have higher average inertia in the beginning this results in 0.
    if (return_plot == FALSE){
      return(dim_number)
    } else {
      return(list("dims" = dim_number, "plot" = screeplot))
    }

  } else {
    cat("Please pick a valid method!")
    stop()
  }
}



#' @param assay Character. The assay from which extract the count matrix for SVD, e.g. "RNA" for Seurat objects or "counts"/"logcounts" for SingleCellExperiments.
#'
#' @rdname pick_dims
#' @export
pick_dims.Seurat <- function(obj, mat = NULL, method="scree_plot", reps=3, python = TRUE, return_plot = FALSE, ..., assay){

  stopifnot("obj doesn't belong to class 'Seurat'" = is(obj, "Seurat"))

  if (method == "elbow_rule") {
    seu <- Seurat::GetAssayData(object = obj, assay = assay, slot = "data")
    seu <- as.matrix(seu)
  } else {
    seu <- NULL
  }

  if ("CA" %in% Seurat::Reductions(obj)){
    caobj <- as.cacomp(obj, assay = assay, recompute = TRUE)
  } else {
    stop("No 'CA' dim. reduction object found. Please run cacomp(seurat_obj, top, coords = FALSE, return_input=TRUE) first.")
  }


  pick_dims.cacomp(obj = caobj,
                    mat = seu,
                    method = method,
                    reps = reps,
                    return_plot = return_plot,
                   python = python)
}


#' @param assay Character. The assay from which extract the count matrix for SVD, e.g. "RNA" for Seurat objects or "counts"/"logcounts" for SingleCellExperiments.
#'
#' @rdname pick_dims
#' @export
pick_dims.SingleCellExperiment <- function(obj, mat = NULL, method="scree_plot", reps=3, python = TRUE, return_plot = FALSE, ..., assay){

  stopifnot("obj doesn't belong to class 'SingleCellExperiment'" = is(obj, "SingleCellExperiment"))

  if (method == "elbow_rule") {
    mat <- SummarizedExperiment::assay(obj, assay)
  } else {
    mat <- NULL
  }

  if ("CA" %in% SingleCellExperiment::reducedDimNames(obj)){
    caobj <- as.cacomp(obj, assay = assay, recompute = TRUE)
  } else {
    stop("No 'CA' dim. reduction object found. Please run cacomp(sce, top, coords = FALSE, return_input=TRUE) first.")
  }


  pick_dims.cacomp(obj = caobj,
                   mat = mat,
                   method = method,
                   reps = reps,
                   return_plot = return_plot,
                   python = python)

}







