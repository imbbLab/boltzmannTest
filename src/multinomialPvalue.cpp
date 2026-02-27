#define RCPP_ARMADILLO_RETURN_COLVEC_AS_VECTOR
#include "RcppArmadillo.h"
#include <cmath>
#include <limits>
// [[Rcpp::depends(RcppArmadillo)]]


// exact multinomial p-value
double ldmultinom(const arma::vec& x, arma::uword N, const arma::vec& logV){
  double logCoeff = std::lgamma(N + 1);
  for(arma::uword i = 0; i < x.n_elem; ++i) {
    logCoeff -= std::lgamma(x(i) + 1);
  }
  // Log of probabilities: sum(xi * log(pi))
  double logProbs = arma::sum(x % logV);
  return logCoeff + logProbs;
}

// [[Rcpp::export]]
double multinomialPvalue_cpp(const arma::mat& G, const arma::vec& eta, const arma::vec& mu, const arma::vec& v, arma::uword N){
  // etamber of outcomes
  arma::uword k = G.n_cols;
  arma::uword total = N + k - 1;
  arma::vec bars = arma::regspace<arma::vec>(0, k - 2);
  arma::vec dataset(k);
  arma::vec mm;
  arma::vec logV = arma::log(v);
  arma::vec effectSize = arma::abs(mu - eta);
  double currentSamplingProbability;
  double pValue = 0;
  arma::sword i;

  while(true){
    // convert bars -> counts
    dataset(0) = bars(0);
    for (i = 1; i < k - 1; ++i)
      dataset(i) = bars(i) - bars(i - 1) - 1;

    dataset(k - 1) = total - bars(k - 2) - 1;
    // compute expectations
    mm = G * dataset / N;

    // compute multinomial probability and add

    if (arma::all(arma::abs(mm - eta) >= effectSize)){
      currentSamplingProbability = ldmultinom(dataset, N, logV);
      if (pValue > currentSamplingProbability){
        pValue = pValue + std::log1p(std::exp(currentSamplingProbability - pValue));
      } else {
        pValue = currentSamplingProbability + std::log1p(std::exp(pValue - currentSamplingProbability));
      }
    }


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
  return std::exp(pValue) - 1;
}

// exact multinomial p-value
double ldmultinom(const arma::uvec& x, arma::uword N, const arma::vec& logV){
  double logCoeff = std::lgamma(N + 1);
  for(arma::uword i = 0; i < x.n_elem; ++i) {
    logCoeff -= std::lgamma(x(i) + 1);
  }
  // Log of probabilities: sum(xi * log(pi))
  double logProbs = arma::sum(x % logV);
  return logCoeff + logProbs;
}



// [[Rcpp::export]]
double multinomialPvalue_incremental_cpp(const arma::mat& G,
                                         const arma::vec& eta,
                                         const arma::vec& mu,
                                         const arma::vec& v,
                                         arma::uword N)
{
  arma::uword k = G.n_cols;
  arma::uword total = N + k - 1;

  arma::uvec bars = arma::regspace<arma::uvec>(0, k - 2);
  arma::uvec counts(k);

  arma::vec effectSize = arma::abs(mu - eta);
  arma::vec logV = arma::log(v);

  arma::vec mm(G.n_rows, arma::fill::zeros);
  double pValueLog = -arma::datum::inf;

  // ---- initialize counts from bars ----
  counts(0) = bars(0);
  for (arma::uword i = 1; i < k - 1; ++i)
    counts(i) = bars(i) - bars(i - 1) - 1;
  counts(k - 1) = total - bars(k - 2) - 1;

  // ---- initialize mm ----
  for (arma::uword j = 0; j < k; ++j)
    mm += G.col(j) * counts(j);

  while (true) {

    arma::vec mm_norm = mm / N;

    if (arma::all(arma::abs(mm_norm - eta) >= effectSize)) {
      double currentLogSamplingProbability =
        ldmultinom(counts, N, logV);
      if (pValueLog > currentLogSamplingProbability){
        pValueLog = pValueLog + std::log1p(std::exp(currentLogSamplingProbability - pValueLog));
      } else {
        pValueLog = currentLogSamplingProbability + std::log1p(std::exp(pValueLog - currentLogSamplingProbability));
      }

    }

    // ---- next combination with incremental update ----
    arma::sword i;
    for (i = k - 2; i >= 0; --i) {

      if (bars(i) != i + total - (k - 1)) {

        // store old counts for bins that will change
        arma::uvec oldCounts = counts;

        bars(i)++;

        for (arma::uword j = i + 1; j < k - 1; ++j)
          bars(j) = bars(j - 1) + 1;

        // recompute only changed counts
        counts(0) = bars(0);
        for (arma::uword j = 1; j < k - 1; ++j)
          counts(j) = bars(j) - bars(j - 1) - 1;
        counts(k - 1) = total - bars(k - 2) - 1;

        // ---- incremental mm update ----
        for (arma::uword j = i; j < k; ++j) {
          arma::sword delta =
            static_cast<arma::sword>(counts(j)) -
            static_cast<arma::sword>(oldCounts(j));

          if (delta != 0)
            mm += delta * G.col(j);
        }

        break;
      }
    }

    if (i < 0) break;
  }

  return std::exp(pValueLog);
}


