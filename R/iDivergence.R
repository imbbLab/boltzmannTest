#' Computes the I-divergence
#'
#' Given two probability vectors p and q compute the I-divergence
#' also known as Kullback-Leibler divergence D(p||q) of p from q
#'
#' @param p a numeric vector that should sum up to one (candidate distribution)
#' @param q a numeric vector that should sum up to one (reference distribution)
#'
#' @returns the value of the I-divergence (a.k.a. Kullback-Leibler divergence)
#' of q to p
#'
#' @examples
#'
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
