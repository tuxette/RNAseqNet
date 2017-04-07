# manage donors and imputation
#
################################################################################
# donor pools, selection of sigma, imputed datasets
################################################################################
#' @import hot.deck
#' @importFrom hot.deck hot.deck
#'
#' @title Impute missing row datasets with multiple hot deck.
#' @export
#'
#' @description
#' \code{imputeHD} performs multiple hot-deck imputation on an input data frame with missing rows. Each missing row
#' is imputed with a unique donor. This method requires an auxiliary dataset to compute similaritities between 
#' individuals and create the pool of donors.
#'
#' @param X n x p numeric matrix containing RNA-seq expression with missing rows
#'  (numeric matrix or data frame)
#' @param Y auxiliary dataset (numeric matrix or data frame)
#' @param sigma threshold for hot-deck imputation, obtained with the 
#'  function\code{\link{chooseSigma}} (numeric, positive)
#' @param m number of replicat in multiple imputation (integer), default to 50.
#' @param seed a single value, interpreted as an in integer, in
#' order to control the hazard
#'
#' @author {Alyssa Imbert, \email{alyssa.imbert@inra.fr}
#'
#' Nathalie Villa-Vialaneix, \email{nathalie.villa-vialaneix@inra.fr}}
#'
#' @references {Picheny, V., Servien, R. and Villa-Vialaneix, N. (2016)
#' Interpretable sparse SIR for digitized functional data. \emph{Preprint}.}
#'
#' @seealso \code{\link{chooseSigma}}, \code{\link{ImputedGLMnetwork}}
#'
#' @examples
#' data(lung)
#' X_full <- lung
#' data(thyroid)
#' #' # create missing observations in lung dataset
#' miss <- 0.2
#' nobs <- dim(X_full)[[1]]
#' miss_ind <- sample(1:nobs, round(miss * nobs), replace = FALSE)
#' X_miss <- X_full
#' X_miss[miss_ind, ] <- NA
#' lung <- na.omit(X_miss)
#' # Multiple hot-deck imputation
#' sigma <- 2
#' don <- imputeHD(lung,thyroid,sigma)
#'
#' @return S3 object of class \code{HDImputed}: a list consisting of
#' \itemize{
#'    \item{\code{donors}}{a list , each element of this list 
#'    contains the donors for every missing observations}
#'    \item{\code{draws}}{ a data frame which indicates which donor was
#'    chosen for each missing samples}
#'    \item{\code{data}}{ a list of m imputed datasets}
#'  }
#'

imputeHD <- function(X, Y, sigma, m = 50, seed=NULL) {
  if (!is.null(seed)) set.seed(seed)

  # prepare data
  if (is.null(rownames(X)) | is.null(rownames(Y))) {
    warning("rownames are not provided: the first ", nrow(X),
            "rows of Y are supposed to correspond to rows of X")
    sample_miss <- c(1:nrow(X), rep(NA, nrow(Y) - nrow(X)))
    merged_data <- data.frame("samples" = sample_miss, Y)
  } else {
    s_size <- length(union(rownames(X), rownames(Y)))
    sample_miss <- c(intersect(rownames(X), rownames(Y)),
                     setdiff(rownames(Y), rownames(X)))
    merged_data <- data.frame("samples" = sample_miss,
                              Y[match(sample_miss, rownames(Y)), ])
    merged_data$samples[is.na(match(merged_data$samples, rownames(X)))] <- NA
  }

  # hot deck imputation
  suppressWarnings(hd_res <- hot.deck(merged_data, sdCutoff = sigma, m = m))
  imputation <- lapply(hd_res$draws, function(alist) match(alist, rownames(X)))
  imputation <- matrix(unlist(imputation), byrow = TRUE,
                       nrow = sum(is.na(merged_data$samples)))
  imputation <- data.frame(imputation)

  # imputed datasets
  X_MI <- sapply(imputation, function(acol) {
    res <- rbind(X, X[acol, ])
    return(res)
  }, simplify = FALSE)

  # donors
  donors <- lapply(hd_res$donors, function(alist) match(alist, rownames(X)))

  res <- list("donors" = donors, "draws" = imputation, "data" = X_MI)
  class(res) <- "HDImputed"
  return(res)
}

## Methods for objects of class HDImputed
#' @title Print 'HDImputed' object.
#' @name HDImputed
#' @exportClass HDImputed
#' @export
#' @aliases summary.HDImputed
#' @aliases print.HDImputed
#' @aliases HDImputed-class
#' @description Print a summary of the result of \code{\link{imputeHD}} (
#' \code{HDImputed} object)
#' @param object a \code{HDImputed} object
#' @param x a \code{HDImputed} object
#' @param ... not used
#' @author {Alyssa Imbert, \email{alyssa.imbert@inra.fr}
#'
#' Nathalie Villa-Vialaneix, \email{nathalie.villa-vialaneix@inra.fr}}
#' @seealso \code{\link{imputeHD}}
#'
summary.HDimputed <- function(object, ...) {
  print(object)
  pool_len <- sapply(object$donors, length)
  cat("\nThe number of donors per missing sample is:\n")
  print(summary(pool_len))
  unique_draws <- apply(object$draws, 2, function(acol) length(unique(acol)))
  cat("The number of unique imputations per sample is:\n")
  print(summary(unique_draws))
}

#' @export
#' @rdname HDImputed
print.HDImputed <- function(x, ...) {
  cat(ncol(x$draws), "imputed datasets for", nrow(x$draws), "missing samples:\n",
      "donors are in '$donors' (list)\n",
      "chosen individuals are in '$draws' (matrix)\n",
      "imputed datasets are in '$data' (list)")
}

createDonors <- function(X, Y, sigma, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # prepare data
  if (is.null(rownames(X)) | is.null(rownames(Y))) {
    warning("rownames are not provided: the first ", nrow(X),
            "rows of Y are supposed to correspond to rows of X")
    sample_miss <- c(1:nrow(X), rep(NA, nrow(Y) - nrow(X)))
    merged_data <- data.frame("samples" = sample_miss, Y)
  } else {
    s_size <- length(union(rownames(X), rownames(Y)))
    sample_miss <- c(intersect(rownames(X), rownames(Y)),
                     setdiff(rownames(Y), rownames(X)))
    merged_data <- data.frame("samples" = sample_miss,
                              Y[match(sample_miss, rownames(Y)), ])
    merged_data$samples[is.na(match(merged_data$samples, rownames(X)))] <- NA
  }

  # hot deck imputation
  suppressWarnings(hd_res <- hot.deck(merged_data, sdCutoff = sigma, m = 1))
  donors <- lapply(hd_res$donors, function(alist) match(alist, rownames(X)))
  return(donors)
}

## TODO: faire documentation => export
#
################################################################################
#  
################################################################################
#' @import stats
#' @title Average variance intra-donors
#' @export
#'
#' @description
#' \code{varIntra} does compute the average variance intra-D(i) (Vintra) for a list of donors associated with each
#' missing individuals.
#'
#' @param X n x p numeric matrix containing RNA-seq expression with missing rows
#'  (numeric matrix or data frame)
#' @param Y auxiliary dataset (numeric matrix or data frame)
#' @param donors the list of vectos which contains the donors for each missing samples in X
#'
#' @author {Alyssa Imbert, \email{alyssa.imbert@inra.fr}
#'
#' Nathalie Villa-Vialaneix, \email{nathalie.villa-vialaneix@inra.fr}}
#'
#' @references {}
#'
#' @seealso \code{\link{chooseSigma}}
#' 
#' @return \code{varIntra} returns a numeric value which corresponds to the average variance for the group of donors 

varIntra <- function(X, Y, donors) {
  # prepare data
  if (is.null(rownames(X)) | is.null(rownames(Y))) {
    warning("rownames are not provided: the first ", nrow(X),
            "rows of Y are supposed to correspond to rows of X")
    indX <- 1:nrow(X)
    indY <- (nrow(X)+1):(nrow(Y))
  } else {
    indX <- intersect(rownames(X), rownames(Y))
    indY <- setdiff(rownames(Y), rownames(X))
  }

  all_varintra <- mapply(function(alist, ind) {
    all_dist <- apply(Y[indX, ][alist, , drop = FALSE], 1, function(arow)
      dist(rbind(arow, Y[ind, ]))^2
    )
    sum(all_dist) / length(all_dist)
  }, donors, indY)
  return(mean(all_varintra))
}

## TODO: faire documentation => export
################################################################################
# choose the threshold sigma for multiple hot-deck imputation
################################################################################
#'
#' @title Select the threshold sigma for hd-MI
#' @export
#'
#' @description
#' \code{chooseSigma} does computes the average variance intra-D(i) for different values of sigma in order to 
#' select a sigma which makes a good trade-off between similarity within the pool of donors and variety 
#' (large enough number of donors in every pool).
#'
#' @param X n x p numeric matrix containing RNA-seq expression with missing rows
#'  (numeric matrix or data frame)
#' @param Y auxiliary dataset (numeric matrix or data frame)
#' @param sigma_list a sequence of increasing positive number of sigma (numeric vector)
#' @param seed a singlue value, interpreted as an integern in order to control the hazard
#'
#' @author {Alyssa Imbert, \email{alyssa.imbert@inra.fr}
#'
#' Nathalie Villa-Vialaneix, \email{nathalie.villa-vialaneix@inra.fr}}
#'
#' @references {}
#'
#' @seealso \code{\link{varIntra}}
#' 
#' @examples 
#' data(lung)
#' X_full <- lung
#' data(thyroid)
#' miss <- 0.2
#' nobs <- dim(X_full)[[1]]
#' miss_ind <- sample(1:nobs, round(miss * nobs), replace = FALSE)
#' X_miss <- X_full
#' X_miss[miss_ind, ] <- NA
#' lung <- na.omit(X_miss)
#' sig <- chooseSigma(lung, thyroid,1:5)
#' plot(sig, type="b")
#' 
#' @return a data frame associating a sigma of the \code{sigma_list} with the corresponding \code{varIntra}
#'
chooseSigma <- function(X, Y, sigma_list, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  donors <- lapply(sigma_list, function(sigma)
    createDonors(X, Y, sigma = sigma)
  )

  all_vars <- lapply(donors, function(alist) varIntra(X, Y, alist)) %>% unlist()

  res <- data.frame("sigma" = sigma_list, "varintra" = all_vars)
  return(res)
}

