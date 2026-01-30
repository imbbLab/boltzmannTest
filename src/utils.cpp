#define RCPP_ARMADILLO_RETURN_COLVEC_AS_VECTOR
#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
bool matrix_hasFullRowRank_cpp(const arma::mat& G, double tol) {
  if (G.n_rows > G.n_cols) return false;
  arma::mat Q, R;
  arma::qr(Q, R, G);

  // Check if all diagonal elements of R are non-zero
  return arma::all(arma::abs(R.diag()) > tol);
}


// [[Rcpp::export]]
bool eta_isFeasible_cpp(const arma::mat& G, const arma::vec& eta, double tol){
  if (G.n_rows != eta.n_rows)
    return false;

  arma::vec row_min = arma::min(G, 1); // per-row minima
  arma::vec row_max = arma::max(G, 1); // per-row maxima

  return arma::all(eta >= row_min - tol) &&
    arma::all(eta <= row_max + tol);
}

// [[Rcpp::export]]
arma::imat startAndBars(arma::uword N, arma::uword k){
  arma::uword total = N + k - 1;
  arma::uword numberCombinations = R::choose(total, k - 1);
  arma::imat res(k, numberCombinations);
  arma::ivec bars(k - 1);
  bars = arma::regspace<arma::ivec>(0, k - 2);
  arma::uword dataSet = 0;
  int i;
  while(true){
    // convert bars -> counts
    res(0, dataSet) = bars(0);
    for (i = 1; i < k - 1; ++i)
      res(i, dataSet) = bars(i) - bars(i - 1) - 1;

    res(k - 1, dataSet) = total - bars(k - 2) - 1;

    dataSet++;

    // next combination of bars
    for (i = k - 2; i >= 0; --i) {
      if (bars(i) != i + total - (k - 1)) {
        bars(i)++;
        for (int j = i + 1; j < k - 1; ++j)
          bars(j) = bars(j - 1) + 1;
        break;
      }
    }
    if (i < 0) break;
  }

  return res;

}
