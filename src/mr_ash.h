
#ifndef MR_ASH_H
#define MR_ASH_H
#include <math.h>
#include <RcppArmadillo.h>

// This depends statement is needed to tell R where to find the
// additional header files.
//
// [[Rcpp::depends(RcppArmadillo)]]
//

// This is 
// [[Rcpp::export]]
arma::uvec random_order (int p, int numiter) {
  arma::uvec o(p * numiter);
  for (int i = 0 ; i < numiter; i++) {
    o.subvec(i * p, (i+1) * p - 1) = arma::randperm(p);
  }
  return o;
}

// This is 
arma::mat outerAddition   (const arma::vec& a, const arma::vec& b) {
  arma::mat A(a.n_elem, b.n_elem);
  A.fill(0);
  A.each_row()          += b.t();
  A.each_col()          += a;
  return A;
}

// This is
void updatebetaj       (const arma::vec& xj, double wj,
                        double& betaj, arma::vec& r,
                        arma::vec& piold, arma::vec& pi,
                        double sigma2, const arma::vec& sa2,
                        const arma::vec& s2inv,
                        double& a1, double& a2,
                        int j, int p,
                        double epstol, 
                        double arma::vec& xtomegaj) { //CHAGED!
  
  // calculate b (CHANGED!)
  double bj           = dot(r, xtomegaj) / wj + betaj ; 
  
  // update r first step
  r                    += xj * betaj ; 
  
  // calculate muj (CHANGED!)
  arma::vec muj         = bj * sa2 * s2inv;
  muj(0)                = 0;
  
  // calculate phij
  arma::vec phij        = log(piold + epstol) + log(s2inv)/2 - (bj * bj * s2inv / 2 );
  phij                  = exp(phij - max(phij));
  phij                  = phij / sum(phij);
  
  
  // pinew
  pi                   += phij / p;
  
  // update betaj
  betaj                 = dot(phij, muj);
  
  // update r second step
  r                    += -xj * betaj;
  
  // precalculate for M-step
  a1                   += dot(phij, muj * muj + sa2 * wj *s2inv); //CHANGED!!!
  a2                   += dot(phij, log(phij + epstol));
  phij(0)               = 0;
  a2                   += -dot(phij, log(s2inv)) / 2;
  
  return;
}

#endif
