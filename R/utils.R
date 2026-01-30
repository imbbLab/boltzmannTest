#' Check whether a matrix has full row rank
#'
#' Uses the QR decomposition to check whether a matrix `G` has full row rank
#'
#' @param G a numeric matrix with #rows <= #columns
#' @param tol double (default `.Machine$double.eps`) tolerance
#'
#' @returns bool with `TRUE` if G has full row rank and `FALSE` if not
#'
#' @export
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
#' @export
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
