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

  if (!is.numeric(G)){
    stop("`G` must be a numeric matrix")
  }
  if (any(!is.finite(G))){
    stop("all entries in `G` must be finite")
  }

  ## coerce eta to a vector
  eta <- as.vector(eta)
  if(!is.numeric(eta)){
    stop("`eta` must be a numeric vector")
  }
  if (any(!is.finite(eta))){
    stop("all entries in `eta` must be finite")
  }

  ## coerce v to a vector
  v <- as.vector(v)
  if(!is.numeric(v)){
    stop("`v` must be a numeric vector")
  }
  if (any(!is.finite(v))){
    stop("all entries in `v` must be finite")
  }

  if(any(v < 0)){
    stop("`v` must be non-negative")
  }


  if (abs(1 - sum(v)) >= sqrt(convTolerance)){
    warning("`v` should sum to 1. Normalizing to 1")
    v <- v / sum(v)
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
    warning("`eta` is not feasible")
    res <- list(
      p = v,
      converged = 0,
      iter = NA,
      error = "`eta` is not feasible"
    )
    res
  }
}
