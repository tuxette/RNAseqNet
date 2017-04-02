## TODO : documentation avec import de glmnet
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
