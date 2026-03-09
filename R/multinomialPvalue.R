#' Compute exact multinomial p-values
#'
#' Computes the exact multinomial p-value using
#'
#' @param G matrix with `nrow(G)` constraints (does not require normalization)
#' for `ncol(G)` entities
#' @param eta vector with one value for each row of `G`
#' @param mu vector with the empirical values
#' @param v reference distribution with one value for each column of `G`
#' @param N samples size
#'
#' @export

multinomialPvalue <- function(G, eta, mu, v, N){
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

  ## coerce mu to a vector
  mu <- as.vector(mu)
  if(!is.numeric(mu) || !is.atomic(mu)){
    stop("`mu` must be a numeric vector")
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
  if (length(eta) != length(mu)){
    stop("`eta` and `mu` must be the same length")
  }
  if (NCOL(G) != length(v)){
    stop("`v` must have one value per column of `G`")
  }
  multinomialPvalue_incremental_cpp(G, eta, mu, v, N)

}
