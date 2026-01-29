#define RCPP_ARMADILLO_RETURN_COLVEC_AS_VECTOR
#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List iProjector_cpp(const arma::mat& G, const arma::vec& eta, const arma::vec& v, arma::uword maxit, double convTolerance){
   // initialize at the reference distribution
   arma::vec iProjection = v;
   arma::vec newProjection = v;
   // jacobian
   arma::mat jacobian;
   // inverse jacobian
   arma::mat inverse;
   // update
   arma::mat update;

   // convergence
   arma::uword converged = 0;
   // error message
   std::string error = "";
   // transposed G
   arma::mat transposedG = G.t();
   arma::uword iter;
   for (iter = 0; iter < maxit; iter++){
     // compute the jacobian
     jacobian = G * arma::diagmat(iProjection) * transposedG;
     // invert the jacobian
     if (arma::inv_sympd(inverse, jacobian, arma::inv_opts::allow_approx)){
       update = transposedG* inverse * (G * iProjection - eta);
       newProjection = iProjection % arma::exp(-update);
       // if this fails then abort
     } else{
       error = "Jacobian is poorly conditioned.";
       break;
     }
     // check for convergence
     // \sum_x newProjection * -update
     if (arma::approx_equal(newProjection, iProjection, "absdiff", convTolerance)){
       iProjection = newProjection;
       // check whether the constraints are met
       if (arma::approx_equal(G * iProjection, eta, "reldiff", convTolerance)){
         converged = 1;
       } else{
         error += " Probabilities not changing anymore, but constraints not met.";
         converged = 2;
       }
       break;
     }
     iProjection = newProjection;
   }
   return  Rcpp::List::create(
     Rcpp::Named("p", iProjection),
     Rcpp::Named("converged", converged),
     Rcpp::Named("iter", iter),
     Rcpp::Named("error", error)
   );

 }
