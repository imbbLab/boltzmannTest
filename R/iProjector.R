#' Compute the I-projection
#'
#' Computes the I-projection of v onto the linear family G p = eta. Throws a warning if any condition cannot be met and returns the
#' reference distribution `v`.
#'
#' @param G matrix with `nrow(G)` constraints (including normalization) for `ncol(G)` entities
#' @param eta vector with one value for each row of `G`
#' @param v reference distribution with one value for each column of `G`
#' @param maxit single integer value (default = 10000L) for the maximimal number of iterations
#' @param convTolerance double value (default = `.Machine$double.eps` for the convergence tolerance
#'
#' @export
iProjector <- function(G, eta, v, maxit = 10000L, convTolerance = .Machine$double.eps){
  ## coerce G to a matrix
  G <- as.matrix(G)
  if (!is.numeric(G) || !is.atomic(G)){
    stop("`G must be a numeric matrix")
  }

  ## coerce eta to a vector
  eta <- as.vector(eta)
  if(!is.numeric(eta) || !is.atomic(eta)){
    stop("`eta` must be a numeric vector")
  }

  ## coerce v to a vector
  v <- as.vector(v)
  if(!is.numeric(v) || !is.atomic(v)){
    stop("`v` must be a numeric vector")
  }

  ## check dimensions
  if (NROW(G) != length(eta)){
    stop("`eta` must have one value per row of `G`")
  }
  if (NCOL(G) != length(v)){
    stop("`v` must have one value per column of `G`")
  }

  ## check for the normalization condition
  if (! any(apply(G, 1, function(x) all(x == 1))) || ! any(eta == 1)){
    stop("normalization condition is missing")
  }
  ## check whether G has full row rank
  if (! matrix_hasFullRowRank(G)){
    stop("`G` has not full row rank")
  }
  ## check whether eta is a feasible solution
  if (eta_isFeasible(G, eta)){
    res <- iProjector_cpp(G, eta, v, maxit, convTolerance)
    class(res) <- "iprojection"
    res
  } else{
    warning("`eta` is not a feasible")
    res <- list(
      p = v,
      converged = 0,
      iter = NA,
      error = ""
    )
    class(res) = "iprojection"
    res
  }
}

#' @export
iDivergence <- function(p, q, tolerance = .Machine$double.eps){
  p <- as.vector(p)
  q <- as.vector(q)
  if (length(p) != length(q)){
    stop("`p` and `q` must have the same length")
  }
  if(!is.numeric(p) || !is.atomic(p)){
    stop("`p` must be a numeric vector")
  }
  if (!is.numeric(q) || !is.atomic(q)){
    stop("`q` must be a numeric vector")
  }
  if (any(p < 0)){
    stop("`p` must be non-negative")
  }
  if (any(q < tolerance)){
    stop("`q` must be positive")
  }
  if(abs(sum(p) - 1) > sqrt(tolerance)){
    stop("`p` must sum to 1")
  }
  if(abs(sum(q) - 1) > sqrt(tolerance)){
    stop("`q` must sum to 1")
  }

  notNull <- which(! p < tolerance)
  sum(p[notNull] * (log(p[notNull]) - log(q[notNull])))
}
