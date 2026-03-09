#' Compute the I-projection
#'
#' Computes the I-projection of v onto the linear family G p = eta. If the
#' conditions cannot be met, the function throws a warning and returns the
#' reference distribution v
#'
#' @param G matrix with `nrow(G)` constraints (including normalization) for
#' `ncol(G)` entities
#' @param eta vector with one value for each row of `G`
#' @param v reference distribution with one value for each column of `G`
#' @param maxit single integer value (default = 10000L) for the maximimal
#' number of iterations
#' @param convTolerance double value (default = `.Machine$double.eps` for the
#' convergence tolerance
#'
#' @returns
#' A list with the following entries
#' \describe{
#'  \item{p}{the I-projection of v onto the linear family described by G p = eta}
#'  \item{converged}{a integer value, with 0 no convergence, 1 converged & conditions are met and 2 converged & conditions are not met}
#'  \item{iter}{the number of iterations used}
#'  \item{error}{an error message if any}
#' }
#' @examples
#'
#' ## projecting the kidney stone data to a study, where the
#' ## treatments A and B received the same proportion of small and large stones
#' data(kidneyStones)
#' outcomes <- outcomes_tibble(kidneyStones)
#' G <- with(
#'   outcomes,
#'   rbind(
#'     norm = 1,
#'     # structural expectation fraction treatment A
#'     treatment_A = treatment == "A",
#'     # ambient expectation fraction small stone size given treatment B
#'     stoneSize_small.treatment_B =
#'       (stoneSize == "small" & treatment == "B") / 0.5,
#'     # ambient expectation fraction small stone size given treatment A
#'     stoneSize_small.treatment_A =
#'       (stoneSize == "small" & treatment == "A") / 0.5
#'   )
#' )
#' eta <- c(
#'   norm = 1,
#'   treatment_A = 0.5,
#'   stoneSize_small.treatment_B = 0.51,
#'   stoneSize_small.treatment_A = 0.51
#'
#' )
#' ## perform I-projection
#' iProj <- iProjector(G, eta = eta, v = empirical(outcomes))
#'
#' ## the old difference in the treatment success rate was
#' with(
#'   outcomes,
#'   sum(
#'     empirical(outcomes) * (treatment == "B") * (success == "yes")
#'   ) / sum(empirical(outcomes) * (treatment == "B")) -
#'   sum(
#'     empirical(outcomes) * (treatment == "A") * (success == "yes")
#'   ) / sum(empirical(outcomes) * (treatment == "A"))
#' )
#'
#' ## the "corrected" difference in the treatment success rate is
#' with(
#'   outcomes,
#'   sum(
#'     iProj$p * (treatment == "B") * (success == "yes")
#'   ) / sum(iProj$p * (treatment == "B")) -
#'   sum(
#'     iProj$p * (treatment == "A") * (success == "yes")
#'   ) / sum(iProj$p * (treatment == "A"))
#' )
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
    res
  } else{
    warning("`eta` is not a feasible")
    res <- list(
      p = v,
      converged = 0,
      iter = NA,
      error = "`eta` is not a feasible"
    )
    res
  }
}

#' Computes the I-divergence
#'
#' Given two probability vectors p and q compute the I-divergence
#' also known as Kullback-Leibler divergence D(p||q) of p from q
#'
#' @param p a numeric vector that should sum up to one (candidate distribution)
#' @param q a numeric vector that should sum up to one (reference distribution)
#'
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
  if(abs(sum(p) - 1) >= sqrt(tolerance)){
    warning("`p` should sum to 1. Normalizing to 1")
    p <- p / sum(p)

  }
  if(abs(sum(q) - 1) >= sqrt(tolerance)){
    warning("`q` must sum to 1. Normalizing to 1")
    q <- q / sum(q)
  }

  notNull <- which(! p < tolerance)
  sum(p[notNull] * (log(p[notNull]) - log(q[notNull])))
}
