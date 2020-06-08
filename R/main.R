# network inference from RNAseq datasets with hd-MI
#
################################################################################
# hd-MI
################################################################################
#' @title Multiple hot-deck imputation and network inference from RNA-seq data.
#'
#' @description
#' \code{imputedGLMnetwork} performs a multiple hot-deck imputation and infers a
#' network for each imputed dataset with a log-linear Poisson graphical model
#' (LLGM).
#'
#' @param X n x p numeric matrix containing RNA-seq expression with missing rows
#' (numeric matrix or data frame)
#' @param Y auxiliary dataset (n' x q numeric matrix or data frame)
#' @param sigma affinity threshold for donor pool
#' @param m number of replicates in multiple imputation (integer). Default to 50
#' @param lambdas a sequence of decreasing positive numbers to control the
#' regularization (numeric vector). Default to \code{NULL}
#' @param B number of iterations for stability selection. Default to 20
#'
#' @export
#'
#' @author {Alyssa Imbert, \email{alyssa.imbert@gmail.com}
#'
#' Nathalie Vialaneix, \email{nathalie.vialaneix@inrae.fr}}
#'
#' @references {Imbert, A., Valsesia, A., Le Gall, C., Armenise, C., Lefebvre,
#' G. Gourraud, P.A., Viguerie, N. and Villa-Vialaneix, N. (2018) Multiple
#' hot-deck imputation for network inference from RNA sequencing data.
#' \emph{Bioinformatics}. \doi{10.1093/bioinformatics/btx819}.}
#'
#' @details When input \code{lambdas} are null the default sequence of
#' \code{\link[glmnet]{glmnet}} for the first model (the one with the first
#' column of \code{count} as the target) is used. A common default sequence is
#' generated for all imputed datasets using this method.
#'
#' @examples
#' data(lung)
#' data(thyroid)
#' nobs <- nrow(lung)
#' miss_ind <- sample(1:nobs, round(0.2 * nobs), replace = FALSE)
#' lung[miss_ind, ] <- NA
#' lung <- na.omit(lung)
#' lambdas <- 4 * 10^(seq(0, -2, length = 10))
#' \dontrun{
#' lung_hdmi <- imputedGLMnetwork(lung, thyroid, sigma = 2, lambdas = lambdas,
#'                                m = 10, B = 5)
#' }
#'
#' @return S3 object of class \code{HDpath}: a list consisting of
#' \itemize{
#'   \item{\code{path}}{ a list of \code{m} data frames, each containing the
#'   adjacency matrix of the inferred network obtained from the corresonding
#'   imputed dataset. The regularization parameter is selected by StARS}
#'   \item{\code{efreq}}{ a numeric matrix of size p x p,  which indicates the
#'   number of times an edge has been predicted among the \code{m} inferred
#'   networks}
#' }

imputedGLMnetwork <- function(X, Y, sigma, m = 50, lambdas = NULL, B = 20) {
  # hot deck imputation
  imputed <- imputeHD(X, Y, sigma, m)

  # lambda selection
  selected_lambda <- lapply(imputed$data, function(counts) {
    test <- stabilitySelection(counts, lambdas, B)
    sel_ind <- test$best
    if (is.null(lambdas)) lambdas <<- test$lambdas
    return(sel_ind)
  }) %>% unlist()
  lambdas <- lambdas[1:max(selected_lambda)]

  inferred_adj <- mapply(function(counts, indl) {
    # network inference
    res_glmnet <- GLMnetwork(counts, lambdas)
    adjacency <- res_glmnet$path[[indl]]
    return(adjacency)
  }, imputed$data, selected_lambda, SIMPLIFY = FALSE)

  edge_distribution <- lapply(inferred_adj, function(ares) {
    edges <- (ares !=0) * t(ares != 0)
  })
  edge_distribution <- Reduce("+", edge_distribution)

  res <- list("path" = inferred_adj, "efreq" = edge_distribution)
  class(res) <- "HDpath"
  return(res)
}

## Methods for objects of class HDpath
#' @title Methods for 'HDpath' objects.
#' @name HDpath
#' @export
#'
#' @aliases summary.HDpath
#' @aliases print.HDpath
#' @aliases HDpath-class
#'
#' @description Methods for the result of \code{\link{imputedGLMnetwork}}
#' (\code{HDpath} object)
#' @param object \code{HDpath} object
#' @param x \code{HDpath} object
#' @param ... not used
#' @author {Alyssa Imbert, \email{alyssa.imbert@gmail.com}
#'
#' Nathalie Vialaneix, \email{nathalie.vialaneix@inrae.fr}}
#'
#' @examples
#' data(lung)
#' data(thyroid)
#' nobs <- nrow(lung)
#' miss_ind <- sample(1:nobs, round(0.2 * nobs), replace = FALSE)
#' lung[miss_ind, ] <- NA
#' lung <- na.omit(lung)
#' lambdas <- 4 * 10^(seq(0, -2, length = 10))
#' \dontrun{
#' lung_hdmi <- imputedGLMnetwork(lung, thyroid, sigma = 2, lambdas = lambdas,
#'                                m = 10, B = 5)
#' plot(lung_hdmi)
#' }
#'
#' @seealso \code{\link{imputedGLMnetwork}}

summary.HDpath <- function(object, ...) {
  print(object)
}

#' @export
#' @rdname HDpath

print.HDpath <- function(x, ...) {
  cat("Object of class 'HDpath'...\n",
      length(x$path), "predicted sets of coefficients from Poisson GLM after HD imputation.")
}

#' @import ggplot2
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_histogram
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 element_text
#' @export
#' @rdname HDpath

plot.HDpath <- function(x, ...) {
  df <- data.frame("Frequency" = x$efreq[upper.tri(x$efreq)])
  p <- ggplot(data = df, aes_string(x = "Frequency")) + geom_histogram() +
    theme_bw() + ggtitle("Distribution of edge frequency") +
    ylab("count") +
    theme(title = element_text(size = 10))

  return(p)
}

################################################################################
# convert the result of imputedGLMnetwork into an 'igraph' object
################################################################################
#' @import igraph
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph simplify
#' @importFrom methods is
#'
#' @title Convert the result of imputedGLMnetwork or a matrix into a network.
#' @export
#'
#' @description
#' \code{GLMnetToGraph} combines the m inferred networks, obtained from m
#' imputed datasets, into a single stable network or convert a matrix of
#' coefficients of a GLM model into a network (non zero coefficients are
#' converted to edges)
#'
#' @param object an object of class \code{HDpath} as obtained from the function
#' \code{\link{imputedGLMnetwork}} or a squared matrix with zero and non zero
#' values
#' @param threshold the percentage of times, among the m imputed networks, that
#' an edge has to be predicted to be in the final network. Used only for objects
#' of class \code{HDpath}. Default to 0.9
#'
#' @author {Alyssa Imbert, \email{alyssa.imbert@gmail.com}
#'
#' Nathalie Vialaneix, \email{nathalie.vialaneix@inrae.fr}}
#'
#' @references {Imbert, A., Valsesia, A., Le Gall, C., Armenise, C., Lefebvre,
#' G. Gourraud, P.A., Viguerie, N. and Villa-Vialaneix, N. (2018) Multiple
#' hot-deck imputation for network inference from RNA sequencing data.
#' \emph{Bioinformatics}. \doi{10.1093/bioinformatics/btx819}.}
#'
#' @seealso \code{\link{imputedGLMnetwork}}, \code{\link[igraph]{igraph}}
#'
#' @examples
#' data(lung)
#' data(thyroid)
#' nobs <- nrow(lung)
#' miss_ind <- sample(1:nobs, round(0.2 * nobs), replace = FALSE)
#' lung[miss_ind, ] <- NA
#' lung <- na.omit(lung)
#' lambdas <- 4 * 10^(seq(0, -2, length = 10))
#' \dontrun{
#' lung_hdmi <- imputedGLMnetwork(lung, thyroid, sigma = 2, lambdas = lambdas,
#'                                m = 10, B = 5)
#' lung_net <- GLMnetToGraph(lung_hdmi, 0.75)
#' lung_net
#' plot(lung_net)
#' }
#'
#' @return an 'igraph' object. See \code{\link[igraph]{igraph}}

GLMnetToGraph <- function(object, threshold = 0.9) {
  if (is(object, "HDpath")) {
    if (threshold < 0 | threshold > 1) {
      stop("'threshold' must be a ratio (between 0 and 1)")
    }
    adj_sel <- object$efreq > round(length(object$path) * threshold)
    net <- graph_from_adjacency_matrix(adj_sel, mode = "undirected")
    net <- simplify(net)
  } else {
    adj_sel <- object != 0
    net <- graph_from_adjacency_matrix(adj_sel, mode = "undirected")
    net <- simplify(net)
  }
  return(net)
}
