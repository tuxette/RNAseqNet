
#' @title Multiple hot-deck imputation and network inference from RNA-seq data
#' 
#' @description
#' \code{ImputedGLMnetwork}performs a multiple hot-deck imputation with an affinity 
#' score and infer a network for each imputed datasets with a log-linear Poisson graphical model
#' 
#'  @param X RNA-seq datasets with missing rows (numeric matrix or data frame)
#'  @param Y  auxiliary data (numeric matrix or data frame)
#'  @param sigma threshold for hot-deck imputation, obtained with the function \code{\link{chooseSigma}} (numeric, positive)
#'  @param m number of replicat in multiple imputation (integer). Default to 50.
#'  @param lambdas a sequence of decresing positive numbers to control
#'  the regularization (numeric vector)
#'  
#' @export  
#'  
#' @author {Alyssa Imbert, \email{alyssa.imbert@inra.fr}
#'
#' Nathalie Villa-Vialaneix, \email{nathalie.villa-vialaneix@inra.fr}}
#' 
#' @return S3 object of class \code{ImputedGLMnetwork} : a list consisiting of
#'  \itemize{
#'        \item{\code{path}}{a list of m data frame, 
#'        each data frame corresponds to a adjacency matrix (one per imputed dataset)}
#'        \item{\code{efreq}}{a numeric matrix of size p x p  which 
#'        is the sum of the m adjacency matrix}
#'    } 

ImputedGLMnetwork <- function(X, Y, sigma, m = 50, lambdas) {
  # hot deck imputation
  imputed <- imputeHD(X, Y, sigma, m)

  # lambda selection TODO: parallel computing... with snow?
  selected_lambda <- lapply(imputed$data, function(counts) {
    test <- stabilitySelection(counts, lambdas)
    #selected <- test < 0.05
    #sel_ind <- which.max(test[selected])
    sel_ind <- test$best
    return(sel_ind)
  }) %>% unlist()

  inferred_adj <- mapply(function(counts, indl) {
    # network inference
    res_glmnet <- GLMnetwork(counts, lambdas)
    adjacency <- res_glmnet[[indl]]
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

summary.HDpath <- function(object, ...) {
  print(object)
}

#' @export
#' @rdname HDpath
print.HDpath <- function(x, ...) {
  cat("Object of class 'HDpath'...\n",
      length(x), "predicted sets of coefficients from Poisson GLM after HD imputation.")
}

#' @import ggplot2
#' @export
#' @rdname HDpath
plot.HDpath <- function(x, ...) {
  df <- data.frame("Frequency" = x$efreq[upper.tri(x$efreq)])
  p <- ggplot(data = df, aes(x = Frequency)) + geom_histogram() +
    theme_bw() + ggtitle("Distribution of edge frequency") +
    ylab("count") +
    theme(title = element_text(size=10))
  print(p)

  return(p)
}

##################################################################################
##
##################################################################################
#' @import igraph
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph simplify
#' @title Combine m network into a single and final network
#' @export
#' 
#' @description
#' \code{GLMnetToGraph} is the last part in multiple imputation : this function combines
#' the m networks obtained from m imputed dataset into a single compromise network.
#' 
#' @param object an object of class \code{ImputedGLMnetwork} as obtained from 
#' the function \code{\link{ImputedGLMnetwork}}
#' @param threshold the percent of stable edges which will be kept in the final network 
#'  
#' @author {Alyssa Imbert, \email{alyssa.imbert@inra.fr}
#'
#' Nathalie Villa-Vialaneix, \email{nathalie.villa-vialaneix@inra.fr}}
#' 
#' @references {}
#'
#' @seealso \code{\link{ImputedGLMnetwork}
#' 
#' @return \code{GLMnetToGraph} returns the compromise network 

GLMnetToGraph <- function(object, threshold) {
  adj_sel <- object$efreq > round(length(object$path) * threshold)
  net <- graph_from_adjacency_matrix(adj_sel, mode = "undirected")
  net <- simplify(net)
  return(net)
}
