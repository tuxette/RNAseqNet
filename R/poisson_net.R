# network inference from RNAseq datasets (LLGM)
#
################################################################################
# network inference
################################################################################
#' @import glmnet
#' @import PoiClaClu
#' @importFrom glmnet glmnet
#' @importFrom PoiClaClu FindBestTransform
#'
#' @title Infer a network from RNA-seq expression.
#' @export
#'
#' @description
#' \code{GLMnetwork} infers a network from RNA-seq expression with the
#' log-linear Poisson graphical model of (Allen and Liu, 2012).
#'
#' @param counts a n x p matrix of RNA-seq expression (numeric matrix or data
#' frame)
#' @param lambdas a sequence of decreasing positive numbers to control the
#' regularization (numeric vector). Default to \code{NULL}
#' @param normalize logical value to normalize predictors in the log-linear
#' Poisson graphical model. If \code{TRUE}, log normalization and scaling are
#' performed prior the model is fit. Default to \code{TRUE}
#'
#' @author {Alyssa Imbert, \email{alyssa.imbert@gmail.com}
#'
#' Nathalie Vialaneix, \email{nathalie.vialaneix@inrae.fr}}
#'
#' @details When input \code{lambdas} are null the default sequence of
#' \code{\link[glmnet]{glmnet}} for the first model (the one with the first
#' column of \code{count} as the target) is used.
#'
#' @seealso \code{\link{stabilitySelection}}
#'
#' @references {Allen, G. and Liu, Z. (2012) A log-linear model for inferring
#' genetic networks from high-throughput sequencing data. In \emph{Proceedings
#' of IEEE International Conference on Bioinformatics and Biomedecine (BIBM)}.}
#'
#' @examples
#' data(lung)
#' lambdas <- 4 * 10^(seq(0, -2, length = 10))
#' ref_lung <- GLMnetwork(lung, lambdas = lambdas)
#'
#' @return S3 object of class \code{GLMnetwork}: a list consisting of
#' \itemize{
#'   \item{\code{lambda}}{ regularization parameters used for LLGM path(vector)}
#'   \item{\code{path}}{ a list having the same length than \code{lambda}. It
#'   contains the estimated coefficients (in a matrix) along the path}
#'  }

GLMnetwork <- function(counts, lambdas = NULL, normalize = TRUE) {
  nvar <- ncol(counts)
  nobs <- nrow(counts)
  # transformation of data (into Poisson distribution)
  alpha <- FindBestTransform(counts)
  counts <- counts^alpha

  # glm fitting for every variable
  all_glms <- lapply(1:ncol(counts), function(avar) {
    target <- counts[ ,avar,drop = FALSE]
    predictors <- counts[ ,-avar]
    if (normalize) {
      predictors <- log(predictors + 1)
      predictors <- scale(predictors, center = FALSE)
    }
    cur_glm <- glmnet(x = predictors, y = target, family = "poisson",
                      lambda = lambdas)

    if ((avar == 1) & is.null(lambdas)) lambdas <<- cur_glm$lambda

    return(cur_glm$beta)
  })

  # re-order result by lambda and fill the diagonal with 1s
  lambda_path <- lapply(seq_along(lambdas), function(indl) {
    all_betas <- lapply(all_glms, function(alist) alist[ ,indl]) %>% unlist()
    all_betas <- matrix(all_betas, byrow = TRUE, nrow = ncol(counts))
    res <- matrix(0, ncol = ncol(counts), nrow = ncol(counts))
    res[upper.tri(res)] <- all_betas[upper.tri(all_betas, diag = TRUE)]
    res[lower.tri(res)] <- all_betas[lower.tri(all_betas)]
    diag(res) <- 1
    colnames(res) <- rownames(res) <- colnames(counts)
    return(res)
  })

  lambda_path <- list("lambda" = lambdas, "path" = lambda_path)

  class(lambda_path) <- "GLMpath"
  return(lambda_path)
}

## Methods for objects of class GLMpath
#' @title Methods for 'GLMpath' objects.
#' @name GLMpath
#' @exportClass GLMpath
#' @export
#'
#' @aliases summary.GLMpath
#' @aliases print.GLMpath
#' @aliases GLMpath-class
#'
#' @description Methods for the result of \code{\link{GLMnetwork}}
#' (\code{GLMpath} object)
#' @param object \code{GLMpath} object
#' @param x \code{GLMpath} object
#' @param ... not used
#' @author {Alyssa Imbert, \email{alyssa.imbert@gmail.com}
#'
#' Nathalie Vialaneix, \email{nathalie.vialaneix@inrae.fr}}
#'
#' @seealso \code{\link{GLMnetwork}}

summary.GLMpath <- function(object, ...) {
  print(object)
}

#' @export
#' @rdname GLMpath
print.GLMpath <- function(x, ...) {
  cat("Object of class 'GLMpath'...\n",
      length(x$lambda), "predicted sets of coefficients from Poisson GLM (LLGM).\n",
      "Regularization parameters are in '$lambda' and coefficient in '$path'.")
}

################################################################################
# StARS for network selection
################################################################################
#' @title Selection of the regularization parameter by StARS (Liu et al., 2010).
#' @export
#'
#' @description
#' \code{stabilitySelection} implements the regularization parameter selection
#' of (Liu et al., 2010) called 'Stability Approach to Regularization Selection'
#' (StARS).
#'
#' @param counts a n x p matrix of RNA-seq expression (numeric matrix or data
#' frame)
#' @param lambdas a sequence of decreasing positive numbers to control the
#' regularization (numeric vector). Default to \code{NULL}
#' @param B number of iterations for stability selection. Default to 20
#'
#' @references {Liu, H., Roeber, K. and Wasserman, L. (2010) Stability approach
#' to regularization selection (StARS) for high dimensional graphical models. In
#' \emph{Proceedings of Neural Information Processing Systems (NIPS 2010)},
#' \strong{23}, 1432-1440, Vancouver, Canada.}
#'
#' @details When input \code{lambdas} are null the default sequence of
#' \code{\link[glmnet]{glmnet}} (see \code{\link{GLMnetwork}} for details).
#'
#' @author {Alyssa Imbert, \email{alyssa.imbert@gmail.com}
#' Nathalie Vialaneix, \email{nathalie.vialaneix@inrae.fr}}
#'
#' @seealso \code{\link{GLMnetwork}}
#'
#' @examples
#' data(lung)
#' lambdas <- 4 * 10^(seq(0, -2, length = 5))
#' stability_lung <- stabilitySelection(lung, lambdas = lambdas, B = 4)
#' \dontrun{plot(stability_lung)}
#'
#' @return S3 object of class \code{stabilitySelection} : a list consisting of
#' \itemize{
#'   \item{\code{lambdas}}{ numeric regularization parameters used for
#'   regularization path}
#'   \item{\code{B}}{ number of iterations for stability selection}
#'   \item{\code{best}}{ index of the regularization parameter selected by StARS
#'   in \code{lambdas}}
#'   \item{\code{variabilities}}{ numeric vector having same length than lambdas
#'   and providing the variability value as defined by StARS along the path}
#' }

stabilitySelection <- function(counts, lambdas = NULL, B = 20) {
  nvar <- ncol(counts)
  nobs <- nrow(counts)

  # collect all subsampling results
  all_res <- lapply(1:B, function(aboot) {
    selected <- sample(1:nobs, round(0.8 * nobs), replace = FALSE)
    # network inference
    res_boot <- GLMnetwork(counts[selected, ], lambdas = lambdas)
    if ((aboot == 1) & is.null(lambdas)) lambdas <<- res_boot$lambda
    return(res_boot$path)
  })
  # order by lambda and create adjacency
  all_res <- lapply(seq_along(lambdas), function(indl) {
    lapply(all_res, function(alist) alist[[indl]] != 0)
  })

  # compute variabilities
  all_variabilities <- lapply(all_res, function(alist) {
    variability <- Reduce("+", alist) / B
    diag(variability) <- 0
    variability <- 4 * sum(variability * (1 - variability))
    variability <- variability / (nvar * (nvar - 1))
  }) %>% unlist()
  all_variabilities <- sapply(seq_along(all_variabilities), function(ind) {
    max(all_variabilities[1:ind])
  })

  selected <- all_variabilities < 0.05
  sel_ind <- which.max(all_variabilities[selected])

  res <- list("lambdas" = lambdas, "B" = B, "best" = sel_ind,
              "variabilities" = all_variabilities)
  class(res) <- "stars"
  return(res)
}

## Methods for objects of class stars
#' @title Methods for 'stars' objects.
#' @name stars
#' @exportClass stars
#' @export
#'
#' @aliases summary.stars
#' @aliases print.stars
#' @aliases plot.stars
#' @aliases stars-class
#'
#' @description Methods for the result of \code{\link{stabilitySelection}}
#' (\code{stars} object)
#' @param object \code{stars} object
#' @param x \code{stars} object
#' @param ... not used
#' @author {Alyssa Imbert, \email{alyssa.imbert@gmail.com}
#'
#' Nathalie Vialaneix, \email{nathalie.vialaneix@inrae.fr}}
#'
#' @seealso \code{\link{stabilitySelection}}

summary.stars <- function(object, ...) {
  print(object)
}

#' @export
#' @rdname stars
print.stars <- function(x, ...) {
  cat("Best lambda found with", x$B, "bootstrap simulations is",
      x$lambdas[x$best], "\n",
      "(variability:", x$variabilities[x$best], ")")
}

#' @import graphics
#' @export
#' @rdname stars
plot.stars <- function(x, ...) {
  args <- list(...)
  args$x <- x$lambdas
  args$y <- x$variabilities

  if (is.null(args$xlab)) args$xlab <- expression(lambda)
  if (is.null(args$ylab)) args$ylab <- "StARS variability"

  do.call(plot, args)
  points(x$lambdas[x$best], x$variabilities[x$best], pch = "+", cex = 1.5,
         col = "red")
}
