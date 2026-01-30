#' @export
boltzmannTest <- function(x, ...) UseMethod("boltzmannTest")

## what do we need
## 1. G the coefficient matrix
## 2. H = eta == vector of hypothesis conditions
## 3. A = mu == vector of observed means
## 4. A' = subset of observed means, confounders?

## boltzmannTest <- function(G, eta, mu, subs)
## 1. project empirical to hypothesis linear family H: G p = eta -> h
## 2. if subs are present project h to to A_k: G_k p = mu_k -> a_k
## 3. if subs are not present a_k = h (in a sense only normalization)
## 4. project a_k to A: G * p = mu -> f

#' @export
boltzmannTest.entities_tibble <-function(entities, G, eta, nested = NULL, strict = FALSE, maxit = 10000L, convTolerance = .Machine$double.eps){
  if (NCOL(G) != length(eta)){
    stop("the number of constraints ")
  }
  if (!is.null(nested)){
    nested <- as.vector(nested)
    if (any(nested) < 1 || any(nested) > NROW(G)){
      stop("the indices of the nested conditions are out of range")
    }
  }
  ## 1. project empirical distribution to hypothesis linear family
  h <- iProjector(G = G, eta = eta, v = empirical(entities), maxit = maxit, convTolerance = convTolerance)
  ## check for convergence
  if(h$converged != 0){
    ## "update" hypothesis distribution
    if (!is.null(nested)){

    }
    notNull <- which(!h$p < convTolerance)
    a <- if(any(h$p < convTolerance)){
      warning("hypothesis conditions lead to zero probabilities")
      iProjector(G[, notNull], G %*% empirical(entities), v = h$p[notNull], maxit = maxit, convTolerance = convTolerance)
    } else{
      res <- list(
        p = empirical(entities),
        converged = 1,
        iter = 0,
        error = ""
      )
      class(res) = "iprojection"
      res
    }

    idiv <- iDivergence(a$p, h$p)

    statistic <- 2 * sampleSize(entities) * idiv
    df = 1
    p.value = pchisq(statistic, df = df, lower.tail = FALSE)



  }




}

