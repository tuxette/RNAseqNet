ImputedGLMnetwork <- function(X, Y, sigma, m = 50, lambdas) {
  # hot deck imputation
  imputed <- imputeHD(X, Y, sigma, m = 50)

  # lambda selection TODO: parallel computing... with snow?
  selected_lambda <- lapply(imputed$data, function(counts) {
    test <- stabilitySelection(counts, lambdas)
    selected <- test < 0.05
    sel_ind <- which.max(test[selected])
    return(sel_ind)
  }) %>% unlist()

  inferred_adj <- mapply(function(counts, indl) {
    # network inference
    res_glmnet <- GLMnetwork(tX_MI, lambdas)
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

#' @export
#' @rdname HDpath
plot.HDpath <- function(x, ...) {
  df <- data.frame("frequency" = x$efreq)
  p <- ggplot(data = df, aes(x = frequency)) + geom_histogram() +
    theme_bw() + ggtitle("Distribution of edge frequency") +
    ylab("count") +
    theme(title = element_text(size=10))
  print(p)

  return(p)
}

GLMnetToGraph <- function(object, threshold) {
  adj_sel <- object$efreq > round(length(object$path) * threshold)
  net <- graph_from_adjacency_matrix(adj_sel, mode = "undirected")
  net <- simplify(net)
  return(net)
}
