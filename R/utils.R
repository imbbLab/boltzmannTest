#' Check whether a matrix has full row rank
#'
#' Uses the QR decomposition to check whether a matrix `G` has full row rank
#'
#' @param G a numeric matrix with #rows <= #columns
#' @param tol double (default `.Machine$double.eps`) tolerance
#'
#' @returns bool with `TRUE` if G has full row rank and `FALSE` if not
#'
#' @keywords internal
matrix_hasFullRowRank <- function(G, tol = .Machine$double.eps){
  ## coerce G to a matrix
  G <- as.matrix(G)
  if (!is.numeric(G) || !is.atomic(G)){
    stop("`G must be a numeric matrix")
  }
  matrix_hasFullRowRank_cpp(G, tol)
}

#' Check whether a vector is a feasible linear combination on the simplex
#'
#' Checks whether a vector `eta` = `G` `p` can be obtained by a vector `p` with p_i >= 0
#' and sum p_i = 1
#'
#' @param G a numeric matrix
#' @param eta a numeric vector with one value per row of G
#' @param tol double (default `.Machine$double.eps`) tolerance
#' @returns bool with `TRUE` if G p = eta is feasible
#'
#' @keywords internal
eta_isFeasible <- function(G, eta, tol = .Machine$double.eps){
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
  ## check dimensions
  if (NROW(G) != length(eta)){
    stop("`eta` must have one value per row of `G`")
  }
  eta_isFeasible_cpp(G, eta, tol)
}

#' Check which entries of two numeric vectors are (approximately) equal
#'
#' Checks which entries of two numeric vector are up to a difference of > `tol` equal.
#' Note that we follow the convention of Armadillo: if the difference of two numbers
#' x and y is smaller or equal to `tol` then the numbers are approximately equal.
#'
#' @param a a numeric vector
#' @param b a numeric vector of the same length as a
#' @param tol double value (default`.Machine$double.eps`) of the tolerance
#' @param reldiff logical (default `FALSE`) whether to compute the relative difference `TRUE`
#' or the absolute difference `FALSE`
#' @returns a logical vector with elements `TRUE` if the entries are approximately equal
#' or `FALSE` if not.
#' @keywords internal
which_approx_equal <- function(a, b, tol = .Machine$double.eps, reldiff = FALSE){
  a <- as.vector(a)
  if (!is.numeric(a) || !is.atomic(a)){
    stop("`a` is not a numeric vector")
  }
  b <- as.vector(b)
  if (!is.numeric(b) || !is.atomic(b)){
    stop("`b` is not a numeric vector")
  }
  if (length(a) != length(b)){
    stop("`b` has not the same number of entries as `a`")
  }
  which_approx_equal_cpp(a,b, tol, reldiff)
}
