#' Computes the I-divergence
#'
#' Given two probability vectors p and q compute the I-divergence
#' also known as Kullback-Leibler divergence D(p||q) of p from q
#'
#' @param p a numeric vector that should sum up to one (candidate distribution)
#' @param q a numeric vector that should sum up to one (reference distribution)
#' @param tolerance the numerical tolerance (default .Machine$double.eps)
#'
#' @returns the value of the I-divergence (a.k.a. Kullback-Leibler divergence)
#' of q to p
#'
#' @examples
#'
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
#' iDivergence(iProj$p, empirical(outcomes))
#'
#'
#' @export
iDivergence <- function(p, q, tolerance = .Machine$double.eps){
  p <- as.vector(p)
  q <- as.vector(q)
  if (length(p) != length(q)){
    stop("`p` and `q` must have the same length")
  }
  if(!is.numeric(p)){
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
    warning("`q` should sum to 1. Normalizing to 1")
    q <- q / sum(q)
  }

  notNull <- which(! p < tolerance)
  sum(p[notNull] * (log(p[notNull]) - log(q[notNull])))
}
