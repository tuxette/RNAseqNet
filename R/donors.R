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
#' \code{imputeHD} does blabla
#'
#' @param X tructruc
#'
#' @author {Alyssa Imbert, \email{alyssa.imbert@inra.fr}
#'
#' Nathalie Villa-Vialaneix, \email{nathalie.villa-vialaneix@inra.fr}}
#'
#' @references {Picheny, V., Servien, R. and Villa-Vialaneix, N. (2016)
#' Interpretable sparse SIR for digitized functional data. \emph{Preprint}.}
#'
#' @seealso \code{\link{chooseSigma}}, \code{\link{netImputeHD}}
#'
#' @examples
#' set.seed(1140)
#' tsteps <- seq(0, 1, length = 50)
#' \dontrun{print(res_ridge)}
#'
#' @return S3 object of class \code{HDImputed}: a list consisting of
#' \itemize{
#'    \item{\code{donors}}{ blabla}
#'    \item{\code{draws}}{ blabla}
#'    \item{\code{data}}{ blabla}
#'  }
#'

imputeHD <- function(X, Y, sigma, m = 50) {
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
chooseSigma <- function(X, Y, sigma_list, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  donors <- lapply(sigma_list, function(sigma)
    createDonors(X, Y, sigma = sigma)
  )

  all_vars <- lapply(donors, function(alist) varIntra(X, Y, alist)) %>% unlist()

  res <- data.frame("sigma" = sigma_list, "varintra" = all_vars)
  return(res)
}

