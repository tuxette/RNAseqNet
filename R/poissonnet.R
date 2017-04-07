## TODO : documentation avec import de glmnet

# donor pools, selection of sigma, imputed datasets
################################################################################
#' @import glmnet
#' @import PoiClaClu
#' @importFrom glmnet glmnet
#' @importFrom PoiClaClu FindBestTransform
#'
#' @title Infer a network from RNA-seq expression
#' @export
#'
#' @description
#' \code{GLMnetwork} infers a network from RNA-seq expression by using a log-linear Poisson
#' graphical model 
#'
#' @param counts a n x p matrix of RNA-seq expression (numeric matrix or data frame)
#' @param lambdas a sequence of decresing positive numbers to control
#'  the regularization (numeric vector)
#' @param normalize logical value to normalize predictors in the log-linear
#' Poisson graphical model bu using a log transformation, default to TRUE
#'
#' @author {Alyssa Imbert, \email{alyssa.imbert@inra.fr}
#'
#' Nathalie Villa-Vialaneix, \email{nathalie.villa-vialaneix@inra.fr}}
#'
#' @seealso \code{\link{stabilitySelection}}
#' 
#'
#' @references{}
#' @return S3 object of class \code{GLMnetwork} ....
#'    

GLMnetwork <- function(counts, lambdas, normalize = TRUE) {
  nvar <- ncol(counts)
  nobs <- nrow(counts)
  # transformation of data (into Poisson distribution)
  alpha <- FindBestTransform(counts)
  tcounts <- counts^alpha

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
    return(res)
  })

  class(lambda_path) <- "GLMpath"
  return(lambda_path)
}

summary.GLMpath <- function(object, ...) {
  print(object)
}

#' @export
#' @rdname GLMpath
print.GLMpath <- function(x, ...) {
  cat("Object of class 'GLMpath'...\n",
      length(x), "predicted sets of coefficients from Poisson GLM.")
}

############################################################
## Selection of lambda
############################################################
#'
#' @import glmnet
#' @import PoiClaClu
#' @importFrom glmnet glmnet
#' @importFrom PoiClaClu FindBestTransform
#'
#' @title Selection of the regularization parameter for high dimensional undirected graph
#' @export
#'
#' @description
#' \code{stabilitySelection} implements the regularization parameter selection by using 
#' Stability Approach to Regularization Selection (StARS).
#'
#' @param counts a n x p matrix of RNA-seq expression (numeric matrix or data frame)
#' @param lambdas a sequence of decresing positive numbers to control
#'  the regularization (numeric vector)
#' @param B number of iterations for stability selection, default to 20
#'
#' @author {Alyssa Imbert, \email{alyssa.imbert@inra.fr}
#'
#' Nathalie Villa-Vialaneix, \email{nathalie.villa-vialaneix@inra.fr}}
#' 
#' @seealso \code{\link{GLMnetwork}}
#' 
#' @examples 
#' data(lung)
#' lambdas <- seq(5, 0.1, -0.1)
#' stab <- stabilitySelection(lung, lambdas)
#' plot(stab)
#'
#' @return S3 object of class \code{stabilitySelection} : a list consisiting of
#'  \itemize{
#'        \item{\code{lambdas}}{ numeric vector used for regularization path}
#'        \item{\code{B}}{number of iterations for stability selection}
#'        \item{\code{best}}{lambda value that gives the optimal network}
#'        \item{\code{variabilities}}{numeric vector ...}
#'    }
#'

stabilitySelection <- function(counts, lambdas, B = 20) {
  nvar <- ncol(counts)
  nobs <- nrow(counts)
  # transformation of data (into Poisson distribution)
  alpha <- FindBestTransform(counts)
  tcounts <- counts^alpha

  # collect all subsampling results
  all_res <- lapply(1:B, function(aboot) {
    selected <- sample(1:nobs, round(0.8 * nobs), replace = FALSE)
    # network inference
    res_boot <- GLMnetwork(tcounts[selected, ], lambdas = lambdas)
    return(res_boot)
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
  if (is.null(args$ylab)) args$ylab <- "stars variability"

  do.call(plot, args)
  points(x$lambdas[x$best], x$variabilities[x$best], pch = "+", cex = 1.5,
         col = "red")
}
