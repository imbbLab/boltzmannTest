#' Compute exact multinomial p-values
#'
#' Computes the exact multinomial p-value using
#'
#' @param G matrix with `nrow(G)` constraints (does not require normalization)
#' for `ncol(G)` entities
#' @param eta vector with one value for each row of `G`
#' @param mu vector with the empirical values
#' @param v hypothesis distribution with one value for each column of `G`
#' @param N samples size
#'
#' @returns
#' the exact multinomial p-value
#'
#' @examples
#'
#' data <- data.frame(
#'   x = c(rep(-1, 4), rep(0, 4), rep(1, 2))
#' )
#' (bt <- boltzmann.test(~x, data = data))
#'
#' ## exact multinomial p-value
#' (exactPvalue <- multinomialPvalue(
#'   G = bt$coefficientMatrix[-1, , drop = FALSE],
#'   eta = bt$hypothesisExpectations[-1],
#'   mu = bt$alternativeExpectations[-1],
#'   v = bt$hypothesisDistribution,
#'   N = bt$sampleSize
#' ))
#'
#'
#' @export

multinomialPvalue <- function(G, eta, mu, v, N){

  if (!is.matrix(G) || !is.numeric(G)){
    stop("`G must be a numeric matrix")
  }


  if(!is.numeric(eta) || !is.vector(eta)){
    stop("`eta` must be a numeric vector")
  }

  if(!is.numeric(mu) || !is.vector(mu)){
    stop("`mu` must be a numeric vector")
  }

  if(!is.numeric(v) || !is.vector(v)){
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
