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
#' \code{imputeHD} performs multiple hot-deck imputation on an input data frame
#' with missing rows. Each missing row is imputed with a unique donor. This
#' method requires an auxiliary dataset to compute similaritities between
#' individuals and create the pool of donors.
#'
#' @param X n x p numeric matrix containing RNA-seq expression with missing rows
#' (numeric matrix or data frame)
#' @param Y auxiliary dataset (n' x q numeric matrix or data frame)
#' @param sigma threshold for hot-deck imputation (numeric, positive)
#' @param m number of replicates in multiple imputation (integer). Default to 50
#' @param seed single value, interpreted as an in integer, used to initialize
#' the random number generation state. Default to \code{NULL} (not used in this
#' case)
#'
#' @details Missing values are identified by matching rownames in \code{X} and
#' \code{Y}. If rownames are not provided the missing rows in \code{X} are
#' supposed to correspond to the last rows of \code{Y}.
#'
#' @author {Alyssa Imbert, \email{alyssa.imbert@inra.fr}
#'
#' Nathalie Villa-Vialaneix, \email{nathalie.villa-vialaneix@inra.fr}}
#'
#' @references {Imbert, A., Valsesia, A., Le Gall, C., Armenise, C., Lefebvre,
#' G. Gourraud, P.A., Viguerie, N. and Villa-Vialaneix, N. (2018) Multiple
#' hot-deck imputation for network inference from RNA sequencing data.
#' \emph{Bioinformatics}. \doi{10.1093/bioinformatics/btx819}.}
#'
#' @seealso \code{\link{chooseSigma}}, \code{\link{imputedGLMnetwork}}
#'
#' @examples
#' data(lung)
#' data(thyroid)
#' nobs <- nrow(lung)
#' miss_ind <- sample(1:nobs, round(0.2 * nobs), replace = FALSE)
#' lung[miss_ind, ] <- NA
#' lung <- na.omit(lung)
#' imputed_lung <- imputeHD(lung, thyroid, sigma = 2)
#'
#' @return S3 object of class \code{HDImputed}: a list consisting of
#' \itemize{
#'   \item{\code{donors}}{ a list. Each element of this list contains the donor
#'   pool for every missing observations}
#'   \item{\code{draws}}{ a data frame which indicates which donor was chosen
#'   for each missing samples}
#'   \item{\code{data}}{ a list of \code{m} imputed datasets}
#' }

imputeHD <- function(X, Y, sigma, m = 50, seed = NULL) {
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
#' @title Methods for 'HDImputed' objects.
#' @name HDImputed
#' @exportClass HDImputed
#' @export
#'
#' @aliases summary.HDImputed
#' @aliases print.HDImputed
#' @aliases HDImputed-class
#'
#' @description Methods for the result of \code{\link{imputeHD}}
#' (\code{HDImputed} object)
#' @param object \code{HDImputed} object
#' @param x \code{HDImputed} object
#' @param ... not used
#' @author {Alyssa Imbert, \email{alyssa.imbert@inra.fr}
#'
#' Nathalie Villa-Vialaneix, \email{nathalie.villa-vialaneix@inra.fr}}
#'
#' @seealso \code{\link{imputeHD}}

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

################################################################################
# variability of a donor pool
################################################################################
#' @importFrom stats dist
#' @title Average intra-donor pool variance.
#' @export
#'
#' @description
#' \code{varIntra} computes the average intra-donor pool variance.
#'
#' @param X n x p numeric matrix containing RNA-seq expression with missing rows
#' (numeric matrix or data frame)
#' @param Y auxiliary dataset (n' x q numeric matrix or data frame)
#' @param donors donor pool (a list, as given \code{$donors} obtained from the
#' function \code{\link{imputeHD}})
#'
#' @author {Alyssa Imbert, \email{alyssa.imbert@inra.fr}
#'
#' Nathalie Villa-Vialaneix, \email{nathalie.villa-vialaneix@inra.fr}}
#'
#' @references {Imbert, A., Valsesia, A., Le Gall, C., Armenise, C., Lefebvre,
#' G. Gourraud, P.A., Viguerie, N. and Villa-Vialaneix, N. (2018) Multiple
#' hot-deck imputation for network inference from RNA sequencing data.
#' \emph{Bioinformatics}. \doi{10.1093/bioinformatics/btx819}.}
#'
#' @seealso \code{\link{imputeHD}}, \code{\link{chooseSigma}}
#'
#' @return \code{varIntra} returns a numeric value which is the average
#' intra-donor pool variance, as described in (Imbert \emph{et al.}, 2018).

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

################################################################################
# diagnostic tool for choosing sigma
################################################################################
#'
#' @title Select the threshold sigma for hd-MI.
#' @export
#'
#' @description
#' \code{chooseSigma} computes the average intra-donor pool variance for
#' different values of sigma. It helps choosing a sigma that makes a good
#' trade-off between homogeneity within the pool of donors and variety (large
#' enough number of donors in every pool).
#'
#' @param X n x p numeric matrix containing RNA-seq expression with missing rows
#' (numeric matrix or data frame)
#' @param Y auxiliary dataset (n' x q numeric matrix or data frame)
#' @param sigma_list a sequence of increasing positive values for sigma (numeric
#' vector)
#' @param seed single value, interpreted as an in integer, used to initialize
#' the random number generation state
#'
#' @author {Alyssa Imbert, \email{alyssa.imbert@inra.fr}
#'
#' Nathalie Villa-Vialaneix, \email{nathalie.villa-vialaneix@inra.fr}}
#'
#' @references {Imbert, A., Valsesia, A., Le Gall, C., Armenise, C., Lefebvre,
#' G. Gourraud, P.A., Viguerie, N. and Villa-Vialaneix, N. (2018) Multiple
#' hot-deck imputation for network inference from RNA sequencing data.
#' \emph{Bioinformatics}. \doi{10.1093/bioinformatics/btx819}.}
#'
#' @details The average intra-donor pool variance is described in (Imbert
#' \emph{et al.}, 2018).
#'
#' @seealso \code{\link{varIntra}}
#'
#' @examples
#' data(lung)
#' data(thyroid)
#' nobs <- nrow(lung)
#' miss_ind <- sample(1:nobs, round(0.2 * nobs), replace = FALSE)
#' lung[miss_ind, ] <- NA
#' lung <- na.omit(lung)
#' sigma_stats <- chooseSigma(lung, thyroid, 1:5)
#' \dontrun{plot(sigma_stats, type = "b")}
#'
#' @return a data frame with the values of sigma and the corresponding
#' intra-donor pool variances

chooseSigma <- function(X, Y, sigma_list, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  donors <- lapply(sigma_list, function(sigma)
    createDonors(X, Y, sigma = sigma)
  )

  all_vars <- lapply(donors, function(alist) varIntra(X, Y, alist)) %>% unlist()

  res <- data.frame("sigma" = sigma_list, "varintra" = all_vars)
  return(res)
}
