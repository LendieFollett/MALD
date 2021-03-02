#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <RcppTN.h>

//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::depends(RcppTN)]]
using namespace Rcpp;

namespace lrf{
std::vector<int> csample_int( std::vector<int> x, 
                              int size,
                              bool replace, 
                              NumericVector prob) {
  std::vector<int> ret = RcppArmadillo::sample(x, size, replace, prob) ;
  return ret ;
}



std::vector<int> repc(std::vector<int> x, NumericVector y) {
  int n = y.size();
  std::vector<int> myvector(sum(y));
  int ind=0;
  for (int i=0; i < n; ++i) {
    std::fill(myvector.begin()+ind, myvector.begin()+ind+y[i], x[i]);
    ind += y[i];
  }
  return myvector;
}

const double log2pi = std::log(2.0 * M_PI);

arma::rowvec upperTriangular(arma::mat X){
  int n = X.n_rows;
  arma::rowvec x(n*(n-1)/2);
  int k = -1;
  for (int i = 0; i < n-1; i++){
    for (int j = i + 1; j < n; j++){
      k = k + 1;
      x[k] = X(i,j);
    }
  }
  return x;
}





arma::vec Mahalanobis(arma::mat x, arma::rowvec center, arma::mat cov) {
  int n = x.n_rows;
  arma::mat x_cen;
  x_cen.copy_size(x);
  for (int i=0; i < n; i++) {
    x_cen.row(i) = x.row(i) - center;
  }
  return sum((x_cen * cov.i()) % x_cen, 1);    
}


arma::vec dmvnorm_arma(arma::mat x, arma::rowvec mean, arma::mat sigma, bool log = false) { 
  arma::vec distval = Mahalanobis(x,  mean, sigma);
  double logdet = sum(arma::log(arma::eig_sym(sigma)));
  arma::vec logretval = -( (x.n_cols * log2pi + logdet + distval)/2  ) ;
  
  if (log) { 
    return(logretval);
  } else { 
    return(exp(logretval));
  }
}

arma::mat armgetOmega(arma::mat X, int t, int k1, int k2){
  arma::vec v = arma::trans(X.submat( t, k1, t, k2 -1));
  arma::mat xsub = diagmat(exp(v));
  return xsub;
}





}


const double log2pi = std::log(2.0 * M_PI);
arma::vec Mahalanobis(arma::mat x, arma::rowvec center, arma::mat cov) {
  int n = x.n_rows;
  arma::mat x_cen;
  x_cen.copy_size(x);
  for (int i=0; i < n; i++) {
    x_cen.row(i) = x.row(i) - center;
  }
  return sum((x_cen * cov.i()) % x_cen, 1);    
}

arma::vec dmvnorm_arma(arma::mat x, arma::rowvec mean, arma::mat sigma, bool log = false) { 
  arma::vec distval = Mahalanobis(x,  mean, sigma);
  double logdet = sum(arma::log(arma::eig_sym(sigma)));
  arma::vec logretval = -( (x.n_cols * log2pi + logdet + distval)/2  ) ;
  
  if (log) { 
    return(logretval);
  } else { 
    return(exp(logretval));
  }
}

// [[Rcpp::export]]
arma::vec dmvnrm_arma(arma::mat x,  
                      arma::rowvec mean,  
                      arma::mat sigma, 
                      bool logd = false) { 
  int n = x.n_rows;
  int xdim = x.n_cols;
  arma::vec out(n);
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
  
  for (int i=0; i < n; i++) {
    arma::vec z = rooti * arma::trans( x.row(i) - mean) ;    
    out(i)      = constants - 0.5 * arma::sum(z%z) + rootisum;     
  }  
  
  if (logd == false) {
    out = exp(out);
  }
  return(out);
}

arma::mat armgetSigma(arma::vec omegat, arma::vec sigma_v, arma::vec rho){
  arma::mat Sigma(4,4);
  Sigma(0,0) = omegat(0);
  Sigma(1,1) = omegat(1);
  Sigma(2,2) = sigma_v(0) * sigma_v(0) * omegat(0);
  Sigma(3,3) = sigma_v(1) * sigma_v(1) * omegat(1);
  Sigma(0,1) = rho(0) * sqrt(omegat(0) * omegat(1));
  Sigma(1,0) = rho(0) * sqrt(omegat(0) * omegat(1));
  Sigma(0,2) = rho(2) * sigma_v(0) * omegat(0);
  Sigma(2,0) = rho(2) * sigma_v(0) * omegat(0);
  Sigma(0,3) = 0;
  Sigma(3,0) = 0;
  Sigma(1,2) = 0;
  Sigma(2,1) = 0;
  Sigma(1,3) = rho(3) * sigma_v(1) * omegat(1);
  Sigma(3,1) = rho(3) * sigma_v(1) * omegat(1);
  Sigma(2,3) = rho(1) * sigma_v(0) * sigma_v(1) * sqrt(omegat(0) * omegat(1));
  Sigma(3,2) = rho(1) * sigma_v(0) * sigma_v(1) * sqrt(omegat(0) * omegat(1));
  return Sigma;
}

double corDet(arma::vec rho){
  arma::mat Sigma(4,4);
  Sigma(0,0) = 1;
  Sigma(1,1) = 1;
  Sigma(2,2) = 1;
  Sigma(3,3) = 1;
  Sigma(0,1) = rho(0);
  Sigma(1,0) = rho(0);
  Sigma(0,2) = rho(2);
  Sigma(2,0) = rho(2);
  Sigma(0,3) = 0;
  Sigma(3,0) = 0;
  Sigma(1,2) = 0;
  Sigma(2,1) = 0;
  Sigma(1,3) = rho(3);
  Sigma(3,1) = rho(3);
  Sigma(2,3) = rho(1);
  Sigma(3,2) = rho(1);
  return det(Sigma);
}


std::vector<int> csample_int( std::vector<int> x, 
                              int size,
                              bool replace, 
                              NumericVector prob) {
  std::vector<int> ret = RcppArmadillo::sample(x, size, replace, prob) ;
  return ret ;
}

// [[Rcpp::export]]
IntegerVector update_delta(arma::mat y, arma::mat x, arma::mat omega, arma::vec xiy1, arma::vec xiy2, arma::mat xic, arma::vec mu, arma::vec theta, arma::vec phi, arma::vec sigma_v ,arma::vec rho, arma::vec lambda){
  int T = y.n_rows;
  IntegerVector N(T);
  double part1_N3=0;  double part1_N2=0; double part1_N1=0; double part1_N0=0;
  
  arma::rowvec epsilon(4);
  arma::mat sigma(4,4);
  arma::vec xic1 = xic.col(0);
  arma::vec xic2 = xic.col(1);
  
  double log_P1; double log_P0; double log_P2; double log_P3;
  
  std::vector<int> N_possible(4) ; 
  std::iota (std::begin(N_possible), std::end(N_possible), 0); 
  //N_possible[2] = 3;
  NumericVector probs(4);
  
  for (int t = 0; t < T; t++){
    sigma = armgetSigma(arma::trans(omega.row(t)),  sigma_v,  rho);
    epsilon[2] = omega(t+1,0) - (theta(0) + phi(0) * (omega(t,0) - theta(0)));
    epsilon[3] = omega(t+1,1) - (theta(1) + phi(1) * (omega(t,1) - theta(1)));
    
    epsilon[0] = y(t,0) - (x(t,0) + mu(0) + xiy1(t));
    epsilon[1] = y(t,1) - (x(t,1) + mu(1));
    part1_N0 = dmvnorm_arma(epsilon,  arma::trans(arma::zeros(4)),  sigma, true)[0];//-.5*diff[0]*pow(proposal,-2);
  
    epsilon[0] = y(t,0) - (x(t,0) + mu(0));
    epsilon[1] = y(t,1) - (x(t,1) + mu(1) + xiy2(t));
    part1_N1 = dmvnorm_arma(epsilon,  arma::trans(arma::zeros(4)),  sigma, true)[0];//-.5*diff[0]*pow(gamma2[i],-1); 
  
    epsilon[0] = y(t,0) - (x(t,0) + mu(0) + xic(t,0));
    epsilon[1] = y(t,1) - (x(t,1) + mu(1) + xic(t,1));
    part1_N2 = dmvnorm_arma(epsilon,  arma::trans(arma::zeros(4)),  sigma, true)[0];//-.5*diff[0]*pow(gamma2[i],-1); 
    
    epsilon[0] = y(t,0) - (x(t,0) + mu(0));
    epsilon[1] = y(t,1) - (x(t,1) + mu(1));
    part1_N3 = dmvnorm_arma(epsilon,  arma::trans(arma::zeros(4)),  sigma, true)[0];//-.5*diff[0]*pow(gamma2[i],-1); 
    
    
    //Rcout << part1_N2 << std::endl;
    log_P0 = part1_N0 + log(lambda[0]);
    log_P1 = part1_N1 + log(lambda[1]);
    log_P2 = part1_N2 + log(lambda[2]);
    log_P3 = part1_N3 + log(lambda[3]);    
    
    // for numerical stability
    double pmx = max(NumericVector::create(log_P0, log_P1, log_P2, log_P3));
    log_P0 = log_P0 - pmx;
    log_P1 = log_P1 - pmx;
    log_P2 = log_P2 - pmx;
    log_P3 = log_P3 - pmx;
    
    //Rcout << log_P0 << " "<<log_P1 << " " << log_P2 <<std::endl;
    probs[0] = exp(log_P0)/(exp(log_P0) + exp(log_P1) + exp(log_P2) + exp(log_P3)); // y1
    probs[1] = exp(log_P1)/(exp(log_P0) + exp(log_P1) + exp(log_P2) + exp(log_P3)); // y2
    probs[2] = exp(log_P2)/(exp(log_P0) + exp(log_P1) + exp(log_P2) + exp(log_P3)); // c
    probs[3] = exp(log_P3)/(exp(log_P0) + exp(log_P1) + exp(log_P2) + exp(log_P3)); // none
    //Rcout << probs << std::endl;
    N[t] = csample_int(N_possible, 
                         1,
                         TRUE, 
                         probs)[0];
    
  }
  
  
  return N;
}

// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  //std::cout << sigma << std::endl;
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

// [[Rcpp::export]]
arma::rowvec pgas_s_cpp(arma::vec xi, double w, double eta, arma::vec sprim, int N){
// Initializations
  int T = xi.n_rows;
  arma::mat s(N,T);
  s.zeros();
  arma::mat sfinal(N,T);
  s.zeros();

  IntegerVector ancestors(N);

  NumericVector logWeights(N);
  double maxlogW;
  NumericVector weights(N);
  NumericVector ancestorWeights(N);
  std::vector<int> Indices(N) ; 
  std::iota (std::begin(Indices), std::end(Indices), 0);  
    
  for (int i=0; i < N-1; i++){
    s(i,0) = R::rexp(1);//#draw from transition
  }
    s(N-1,0) = sprim[0];
// Weights update
  for (int i=0; i < N; i++){
    logWeights[i] = R::dnorm(xi[0], w * s(i, 0), eta * sqrt(s(i, 0)), true);
  }
    maxlogW = max(logWeights);
    weights = exp(logWeights - maxlogW);
      
// t >= 2
for (int t=1; t < T; t++) {
      //for (t in c(2:T)){
    ancestors = lrf::csample_int(Indices, 
                            N,
                            TRUE, 
                            weights);    
        
        //s[1:N-1,t] = rexp(N-1, rate = 1)//#draw from transition
   //s.submat(0,t,N-2,t)=R::rexp(N-1);
   for (int i=0; i < N-1; i++){
     s(i,t) = R::rexp(1);//#draw from transition
   }
   s(N-1,t) = sprim[t];

// Ancestor sampling: no change from weights
    // ancestorWeights = weights;
    // ancestors[N-1] = csample_int(Indices, 
    //                     1,
    //                     TRUE, 
    //                     ancestorWeights)[0]; 
          
    // Weights update
    for (int i=0; i < N; i++){
      logWeights[i] = R::dnorm(xi[t], w * s(i, t), eta * sqrt(s(i, t)), true);
    }
    double maxlogW = max(logWeights);
    weights = exp(logWeights - maxlogW);
    for (int i=0;i < N; i++){
        //s[,1:(t-1)] = s[ancestors,1:(t-1)]
        sfinal.submat(i,0,i,t-1) = s.submat(ancestors[i], 0, ancestors[i], t-1);
    }
    sfinal.col(t) = s.col(t);
    s = sfinal;
  
    }
      //ind = sample(c(1:N),size = 1,prob =weights)
      double ind = csample_int(Indices, 
                        1,
                        TRUE, 
                        weights)[0]; 
      return s.row(ind);
      
}

// [[Rcpp::export]]
arma::rowvec pgas_sc_cpp(arma::mat xi, arma::rowvec w_c, arma::mat Sigma_c, arma::vec sprim, int N){
  // Initializations
  int T = xi.n_rows;
  arma::mat s(N,T);
  s.zeros();
  arma::mat sfinal(N,T);
  s.zeros();
  arma::mat sigma_temp(2,2);
  sigma_temp.zeros();
  
  IntegerVector ancestors(N);
  
  NumericVector logWeights(N);
  double maxlogW;
  NumericVector weights(N);
  NumericVector ancestorWeights(N);
  std::vector<int> Indices(N) ; 
  std::iota (std::begin(Indices), std::end(Indices), 0);  
  
  // Sample v using CSMC
  // t = 2 (first time point for the xi variables)
  for (int i=0; i < N-1; i++){
    s(i,0) = R::rexp(1);//#draw from transition 
  }
  //#draw from transition
  s(N-1,0) = sprim[0];
  // Weights update
  //logWeights = logG_s(s.col(0), xi[0], w, eta,N);  
  for (int i=0; i < N; i++){
    sigma_temp(1,1) = Sigma_c(1,1)*s(i,0);
    sigma_temp(0,0) = Sigma_c(0,0)*s(i,0);
    sigma_temp(1,0) = Sigma_c(1,0)*s(i,0);
    sigma_temp(0,1) = Sigma_c(0,1)*s(i,0);
    logWeights[i] = dmvnorm_arma(xi.row(0),  s(i,0)*w_c,  sigma_temp, true)[0];
  }
  maxlogW = max(logWeights);
  weights = exp(logWeights - maxlogW);
  
  // t >= 2
  for (int t=1; t < T; t++) {
    //for (t in c(2:T)){
    ancestors = lrf::csample_int(Indices, 
                                 N,
                                 TRUE, 
                                 weights);    
    
    //s[1:N-1,t] = rexp(N-1, rate = 1)//#draw from transition
    //s.submat(0,t,N-2,t)=R::rexp(N-1);
    for (int i=0; i < N-1; i++){
      s(i,t) = R::rexp(1);//#draw from transition 
    }
    s(N-1,t) = sprim[t];
    
    // Ancestor sampling: no change from weights
    // ancestorWeights = weights;
    // ancestors[N-1] = csample_int(Indices, 
    //                              1,
    //                              TRUE, 
    //                              ancestorWeights)[0]; 
    
    // Weights update
    for (int i=0; i < N; i++){
      sigma_temp(1,1) = Sigma_c(1,1)*s(i,t);
      sigma_temp(0,0) = Sigma_c(0,0)*s(i,t);
      sigma_temp(1,0) = Sigma_c(1,0)*s(i,t);
      sigma_temp(0,1) = Sigma_c(0,1)*s(i,t);
      logWeights[i] = dmvnorm_arma(xi.row(t),  s(i,t)*w_c,  sigma_temp, true)[0];
    }
    maxlogW = max(logWeights);
    weights = exp(logWeights - maxlogW);
    for (int i=0;i < N; i++){
      //s[,1:(t-1)] = s[ancestors,1:(t-1)]
      sfinal.submat(i,0,i,t-1) = s.submat(ancestors[i], 0, ancestors[i], t-1);
    }
    sfinal.col(t) = s.col(t);
    s = sfinal;
    
  }
  //ind = sample(c(1:N),size = 1,prob =weights)
  double ind = csample_int(Indices, 
                           1,
                           TRUE, 
                           weights)[0]; 
  return s.row(ind);
  
}

// [[Rcpp::export]]
arma::rowvec pgas_xiy1_cpp(arma::mat y, arma::mat x, arma::mat omega, arma::vec mu, arma::vec theta, arma::vec phi, arma::vec sigma_v, arma::vec rho, arma::vec xiy1prim, arma::vec xiy2, arma::mat xic, arma::vec Ny1, arma::vec Ny2, arma::vec Nc, double w, double eta, arma::vec sy1, int N){
  // Initializations
  int T = xiy1prim.n_rows;
  arma::mat xiy1(N,T);
  xiy1.zeros();
  arma::mat xiy1final(N,T);
  xiy1final.zeros();
  arma::mat sigma(4,4);
  arma::vec resid(4);

  double Z;
  arma::mat m_dat, s_dat;
  double m_pri, s_pri;
  double dat_ll, pri_ll;
  double p = 0.5;

  arma::mat sigma12;
  arma::mat sigma22;
  arma::mat sigma21;
  arma::uvec idx = {1,2,3};

  IntegerVector ancestors(N);

  NumericVector logWeights(N);
  double maxlogW;
  NumericVector weights(N);
  NumericVector ancestorWeights(N);
  std::vector<int> Indices(N) ;
  std::iota (std::begin(Indices), std::end(Indices), 0);

  sigma = armgetSigma(arma::trans(omega.row(0)),  sigma_v,  rho);
  sigma12 = sigma.cols(idx);
  sigma12 = sigma12.row(0);
  sigma22 = sigma.rows(idx);
  sigma22 = sigma22.cols(idx);
  sigma21 = sigma.rows(idx);
  sigma21 = sigma21.col(0);
  resid(0) = y(0,1) - (x(0,1) + mu(1));
  resid(1) = omega(1,0) - (theta(0) + phi(0) * (omega(0,0) - theta(0)));
  resid(2) = omega(1,1) - (theta(1) + phi(1) * (omega(0,1) - theta(1)));
  resid(3) = y(0,0) - (x(0,0) + mu(0));
  s_dat = sigma(0,0) - sigma12 * arma::inv(sigma22) * sigma21;
  m_dat = resid(3) - sigma12 * arma::inv(sigma22) * resid.head(3);

  m_pri = w * sy1[0];
  s_pri = eta * sqrt(sy1[0]);
  for (int i=0; i < N-1; i++){
    Z = R::runif(0,1);
    if (Z > p){
      xiy1(i,0) = R::rnorm(m_dat(0,0), sqrt(s_dat(0,0)));//#draw from transition
    } else {
      xiy1(i,0) = R::rnorm(m_pri, s_pri);//#draw from transition
    }
  }
  xiy1(N-1,0) = xiy1prim[0];
  // Weights update
  for (int i=0; i < N; i++){
    dat_ll = R::dnorm(xiy1(i,0), m_dat(0,0), sqrt(s_dat(0,0)), true);
    pri_ll = R::dnorm(xiy1(i,0), m_pri, s_pri, true);
    logWeights(i) = Ny1[0] * dat_ll + pri_ll - std::max(dat_ll,pri_ll) - log(p + (1 - p) * exp(-std::abs(pri_ll - dat_ll)));
  }
  maxlogW = max(logWeights);
  weights = exp(logWeights - maxlogW);

  // t >= 2
  for (int t=1; t < T; t++) {
    //for (t in c(2:T)){
    ancestors = lrf::csample_int(Indices,
                                 N,
                                 TRUE,
                                 weights);

    sigma = armgetSigma(arma::trans(omega.row(t)),  sigma_v,  rho);
    sigma12 = sigma.cols(idx);
    sigma12 = sigma12.row(0);
    sigma22 = sigma.rows(idx);
    sigma22 = sigma22.cols(idx);
    sigma21 = sigma.rows(idx);
    sigma21 = sigma21.col(0);
    resid(0) = y(t,1) - (x(t,1) + mu(1));
    resid(1) = omega(t+1,0) - (theta(0) + phi(0) * (omega(t,0) - theta(0)));
    resid(2) = omega(t+1,1) - (theta(1) + phi(1) * (omega(t,1) - theta(1)));
    resid(3) = y(t,0) - (x(t,0) + mu(0));
    s_dat = sigma(0,0) - sigma12 * arma::inv(sigma22) * sigma21;
    m_dat = resid(3) - sigma12 * arma::inv(sigma22) * resid.head(3);

    m_pri = w * sy1[t];
    s_pri = eta * sqrt(sy1[t]);

    for (int i=0; i < N-1; i++){
      Z = R::runif(0,1);
      if (Z > p){
        xiy1(i,t) = R::rnorm(m_dat(0,0), sqrt(s_dat(0,0)));//#draw from transition
      } else {
        xiy1(i,t) = R::rnorm(m_pri, s_pri);//#draw from transition
      }
    }
    xiy1(N-1,t) = xiy1prim[t];

    // Ancestor sampling: No change to Weights
    // ancestorWeights = weights;
    // ancestors[N-1] = csample_int(Indices,
    //                              1,
    //                              TRUE,
    //                              ancestorWeights)[0];

    // Weights update
    for (int i=0; i < N; i++){
      dat_ll = R::dnorm(xiy1(i,t), m_dat(0,0), sqrt(s_dat(0,0)), true);
      pri_ll = R::dnorm(xiy1(i,t), m_pri, s_pri, true);
      logWeights(i) = Ny1[t] * dat_ll + pri_ll - std::max(dat_ll,pri_ll) - log(p + (1 - p) * exp(-std::abs(pri_ll - dat_ll)));
    }
    maxlogW = max(logWeights);
    weights = exp(logWeights - maxlogW);
    for (int i=0;i < N; i++){
      //s[,1:(t-1)] = s[ancestors,1:(t-1)]
      xiy1final.submat(i,0,i,t-1) = xiy1.submat(ancestors[i], 0, ancestors[i], t-1);
    }
    xiy1final.col(t) = xiy1.col(t);
    xiy1 = xiy1final;
  }
  //ind = sample(c(1:N),size = 1,prob =weights)
  double ind = csample_int(Indices,
                           1,
                           TRUE,
                           weights)[0];
  return xiy1.row(ind);
}

// [[Rcpp::export]]
arma::rowvec pgas_xiy2_cpp(arma::mat y, arma::mat x, arma::mat omega, arma::vec mu, arma::vec theta, arma::vec phi, arma::vec sigma_v, arma::vec rho, arma::vec xiy1, arma::vec xiy2prim, arma::mat xic, arma::vec Ny1, arma::vec Ny2, arma::vec Nc, double w, double eta, arma::vec sy2, int N){
  // Initializations
  int T = xiy2prim.n_rows;
  arma::mat xiy2(N,T);
  xiy2.zeros();
  arma::mat xiy2final(N,T);
  xiy2final.zeros();
  arma::mat sigma(4,4);
  arma::vec resid(4);

  double Z;
  arma::mat m_dat, s_dat;
  double m_pri, s_pri;
  double dat_ll, pri_ll;
  double p = 0.5;

  arma::mat sigma12;
  arma::mat sigma22;
  arma::mat sigma21;
  arma::uvec idx = {0,2,3};

  IntegerVector ancestors(N);

  NumericVector logWeights(N);
  double maxlogW;
  NumericVector weights(N);
  NumericVector ancestorWeights(N);
  std::vector<int> Indices(N) ;
  std::iota (std::begin(Indices), std::end(Indices), 0);

  sigma = armgetSigma(arma::trans(omega.row(0)),  sigma_v,  rho);
  sigma12 = sigma.cols(idx);
  sigma12 = sigma12.row(1);
  sigma22 = sigma.rows(idx);
  sigma22 = sigma22.cols(idx);
  sigma21 = sigma.rows(idx);
  sigma21 = sigma21.col(1);
  resid(0) = y(0,0) - (x(0,0) + mu(0));
  resid(1) = omega(1,0) - (theta(0) + phi(0) * (omega(0,0) - theta(0)));
  resid(2) = omega(1,1) - (theta(1) + phi(1) * (omega(0,1) - theta(1)));
  resid(3) = y(0,1) - (x(0,1) + mu(1));
  s_dat = sigma(1,1) - sigma12 * arma::inv(sigma22) * sigma21;
  m_dat = resid(3) - sigma12 * arma::inv(sigma22) * resid.head(3);

  m_pri = w * sy2[0];
  s_pri = eta * sqrt(sy2[0]);
  for (int i=0; i < N-1; i++){
    Z = R::runif(0,1);
    if (Z > p){
      xiy2(i,0) = R::rnorm(m_dat(0,0), sqrt(s_dat(0,0)));//#draw from transition
    } else {
      xiy2(i,0) = R::rnorm(m_pri, s_pri);//#draw from transition
    }
  }
  xiy2(N-1,0) = xiy2prim[0];
  // Weights update
  for (int i=0; i < N; i++){
    dat_ll = R::dnorm(xiy2(i,0), m_dat(0,0), sqrt(s_dat(0,0)), true);
    pri_ll = R::dnorm(xiy2(i,0), m_pri, s_pri, true);
    logWeights(i) = Ny2[0] * dat_ll + pri_ll - std::max(dat_ll,pri_ll) - log(p + (1 - p) * exp(-std::abs(pri_ll - dat_ll)));
  }
  maxlogW = max(logWeights);
  weights = exp(logWeights - maxlogW);

  // t >= 2
  for (int t=1; t < T; t++) {
    //for (t in c(2:T)){
    ancestors = lrf::csample_int(Indices,
                                 N,
                                 TRUE,
                                 weights);

    sigma = armgetSigma(arma::trans(omega.row(t)),  sigma_v,  rho);
    sigma12 = sigma.cols(idx);
    sigma12 = sigma12.row(1);
    sigma22 = sigma.rows(idx);
    sigma22 = sigma22.cols(idx);
    sigma21 = sigma.rows(idx);
    sigma21 = sigma21.col(1);
    resid(0) = y(t,0) - (x(t,0) + mu(0));
    resid(1) = omega(t+1,0) - (theta(0) + phi(0) * (omega(t,0) - theta(0)));
    resid(2) = omega(t+1,1) - (theta(1) + phi(1) * (omega(t,1) - theta(1)));
    resid(3) = y(t,1) - (x(t,1) + mu(1));
    s_dat = sigma(1,1) - sigma12 * arma::inv(sigma22) * sigma21;
    m_dat = resid(3) - sigma12 * arma::inv(sigma22) * resid.head(3);

    m_pri = w * sy2[t];
    s_pri = eta * sqrt(sy2[t]);

    for (int i=0; i < N-1; i++){
      Z = R::runif(0,1);
      if (Z > p){
        xiy2(i,t) = R::rnorm(m_dat(0,0), sqrt(s_dat(0,0)));//#draw from transition
      } else {
        xiy2(i,t) = R::rnorm(m_pri, s_pri);//#draw from transition
      }
    }
    xiy2(N-1,t) = xiy2prim[t];

    // Ancestor sampling: No change to Weights
    // ancestorWeights = weights;
    // ancestors[N-1] = csample_int(Indices,
    //                              1,
    //                              TRUE,
    //                              ancestorWeights)[0];

    // Weights update
    for (int i=0; i < N; i++){
      dat_ll = R::dnorm(xiy2(i,t), m_dat(0,0), sqrt(s_dat(0,0)), true);
      pri_ll = R::dnorm(xiy2(i,t), m_pri, s_pri, true);
      logWeights(i) = Ny2[t] * dat_ll + pri_ll - std::max(dat_ll,pri_ll) - log(p + (1 - p) * exp(-std::abs(pri_ll - dat_ll)));
    }
    maxlogW = max(logWeights);
    weights = exp(logWeights - maxlogW);
    for (int i=0;i < N; i++){
      //s[,1:(t-1)] = s[ancestors,1:(t-1)]
      xiy2final.submat(i,0,i,t-1) = xiy2.submat(ancestors[i], 0, ancestors[i], t-1);
    }
    xiy2final.col(t) = xiy2.col(t);
    xiy2 = xiy2final;
  }
  //ind = sample(c(1:N),size = 1,prob =weights)
  double ind = csample_int(Indices,
                           1,
                           TRUE,
                           weights)[0];
  return xiy2.row(ind);
}

// [[Rcpp::export]]
arma::mat pgas_xic_cpp(arma::mat y, arma::mat x, arma::mat omega, arma::vec mu, arma::vec theta, arma::vec phi, arma::vec sigma_v, arma::vec rho, arma::vec xiy1, arma::vec xiy2, arma::mat xicprim, arma::vec Ny1, arma::vec Ny2, arma::vec Nc, arma::vec wc, arma::vec sigmac, double rhoc, arma::vec sc, int N){
  // Initializations
  int T = xicprim.n_rows;
  arma::cube xic(T,2,N);
  xic.zeros();
  arma::cube xicfinal(T,2,N);
  xicfinal.zeros();
  arma::mat xi(1,2);
  xi.zeros();
  arma::mat sigma(4,4);
  arma::vec resid(4);
  arma::rowvec epsilon(2);

  double Z;
  double dat_ll, pri_ll;
  double p = 0.5;
  arma::mat m_dat, s_dat;
  arma::vec m_pri(2);
  arma::mat s_pri(2,2);
  arma::uvec idx0 = {0,1};
  arma::uvec idx1 = {2,3};

  arma::mat sigma11;
  arma::mat sigma12;
  arma::mat sigma22;
  arma::mat sigma21;

  IntegerVector ancestors(N);

  NumericVector logWeights(N);
  double maxlogW;
  NumericVector weights(N);
  NumericVector ancestorWeights(N);
  std::vector<int> Indices(N) ;
  std::iota (std::begin(Indices), std::end(Indices), 0);

  m_pri(0) = wc(0) * sc(0);
  m_pri(1) = wc(1) * sc(0);
  s_pri(0, 0) = sigmac(0) * sigmac(0) * sc(0);
  s_pri(1, 1) = sigmac(1) * sigmac(1) * sc(0);
  s_pri(0, 1) = rhoc * sigmac(0) * sigmac(1) * sc(0);
  s_pri(1, 0) = rhoc * sigmac(0) * sigmac(1) * sc(0);

  sigma = armgetSigma(arma::trans(omega.row(0)), sigma_v, rho);
  sigma11 = sigma.cols(idx0);
  sigma11 = sigma11.rows(idx0);
  sigma12 = sigma.cols(idx0);
  sigma12 = sigma12.rows(idx1);
  sigma22 = sigma.rows(idx1);
  sigma22 = sigma22.cols(idx1);
  sigma21 = sigma.rows(idx0);
  sigma21 = sigma21.cols(idx1);
  s_dat = sigma11 - sigma12 * arma::inv(sigma22) * sigma21;

  resid(0) = y(0,0) - (x(0,0) + mu(0));
  resid(1) = y(0,1) - (x(0,1) + mu(1));
  resid(2) = omega(1,0) - (theta(0) + phi(0) * (omega(0,0) - theta(0)));
  resid(3) = omega(1,1) - (theta(1) + phi(1) * (omega(0,1) - theta(1)));

  m_dat = resid.head(2) - sigma12 * arma::inv(sigma22) * resid.tail(2);

  for (int i=0; i < N-1; i++){
    Z = R::runif(0,1);
    if (Z > p){
      xi = mvrnormArma(1, m_dat, s_dat);
    } else {
      xi = mvrnormArma(1, m_pri, s_pri);
    }
    xic(0,0,i) = xi(0,0);
    xic(0,1,i) = xi(0,1);
  }
  xic(0,0,N-1) = xicprim(0,0);
  xic(0,1,N-1) = xicprim(0,1);
  // Weights update
  for (int i=0; i < N; i++){
    epsilon(0) = xic(0,0,i);
    epsilon(1) = xic(0,1,i);
    dat_ll = dmvnorm_arma(epsilon, arma::trans(m_dat), s_dat, true)[0];
    //std::cout << "dat_ll = " << dat_ll << std::endl;
    pri_ll = dmvnorm_arma(epsilon, arma::trans(m_pri), s_pri, true)[0];
    //std::cout << "pri_ll = " << pri_ll << std::endl;
    logWeights(i) = Nc[0] * dat_ll + pri_ll - std::max(dat_ll,pri_ll) - log(p + (1 - p) * exp(-std::abs(pri_ll - dat_ll)));
  }
  maxlogW = max(logWeights);
  weights = exp(logWeights - maxlogW);

  // t >= 2
  for (int t=1; t < T; t++) {
    //for (t in c(2:T)){
    ancestors = lrf::csample_int(Indices,
                                 N,
                                 TRUE,
                                 weights);

    m_pri(0) = wc(0) * sc(t);
    m_pri(1) = wc(1) * sc(t);
    s_pri(0, 0) = sigmac(0) * sigmac(0) * sc(t);
    s_pri(1, 1) = sigmac(1) * sigmac(1) * sc(t);
    s_pri(0, 1) = rhoc * sigmac(0) * sigmac(1) * sc(t);
    s_pri(1, 0) = rhoc * sigmac(0) * sigmac(1) * sc(t);

    sigma = armgetSigma(arma::trans(omega.row(t)), sigma_v, rho);
    sigma11 = sigma.cols(idx0);
    sigma11 = sigma11.rows(idx0);
    sigma12 = sigma.cols(idx0);
    sigma12 = sigma12.rows(idx1);
    sigma22 = sigma.rows(idx1);
    sigma22 = sigma22.cols(idx1);
    sigma21 = sigma.rows(idx0);
    sigma21 = sigma21.cols(idx1);
    s_dat = sigma11 - sigma12 * arma::inv(sigma22) * sigma21;

    resid(0) = y(t,0) - (x(t,0) + mu(0));
    resid(1) = y(t,1) - (x(t,1) + mu(1));
    resid(2) = omega(t+1,0) - (theta(0) + phi(0) * (omega(t,0) - theta(0)));
    resid(3) = omega(t+1,1) - (theta(1) + phi(1) * (omega(t,1) - theta(1)));

    m_dat = resid.head(2) - sigma12 * arma::inv(sigma22) * resid.tail(2);

    for (int i=0; i < N-1; i++){
      Z = R::runif(0,1);
      if (Z > p){
        xi = mvrnormArma(1, m_dat, s_dat);
      } else {
        xi = mvrnormArma(1, m_pri, s_pri);
      }
      xic(t,0,i) = xi(0,0);
      xic(t,1,i) = xi(0,1);
    }
    xic(t,0,N-1) = xicprim(t,0);
    xic(t,1,N-1) = xicprim(t,1);

    // Ancestor sampling: No change for weights
    // ancestorWeights = weights;
    // ancestors[N-1] = csample_int(Indices,
    //                              1,
    //                              TRUE,
    //                              ancestorWeights)[0];

    // Weights update
    for (int i=0; i < N; i++){
      epsilon(0) = xic(t,0,i);
      epsilon(1) = xic(t,1,i);
      dat_ll = dmvnorm_arma(epsilon, arma::trans(m_dat), s_dat, true)[0];
      //std::cout << "dat_ll = " << dat_ll << std::endl;
      pri_ll = dmvnorm_arma(epsilon, arma::trans(m_pri), s_pri, true)[0];
      //std::cout << "pri_ll = " << pri_ll << std::endl;
      logWeights(i) = Nc[t] * dat_ll + pri_ll - std::max(dat_ll,pri_ll) - log(p + (1 - p) * exp(-std::abs(pri_ll - dat_ll)));
    }
    maxlogW = max(logWeights);
    weights = exp(logWeights - maxlogW);
    for (int i=0;i < N; i++){
      //s[,1:(t-1)] = s[ancestors,1:(t-1)]
      xicfinal.subcube(0,0,i,t-1,1,i) = xic.subcube(0,0,ancestors[i],t-1,1,ancestors[i]);
    }
    xicfinal.row(t) = xic.row(t);
    xic = xicfinal;
  }
  //ind = sample(c(1:N),size = 1,prob =weights)
  double ind = csample_int(Indices,
                           1,
                           TRUE,
                           weights)[0];
  return xic.slice(ind);
}

// // [[Rcpp::export]]
// arma::vec update_xiy1(arma::mat y, arma::mat x, arma::mat omega, arma::vec mu, arma::vec theta, arma::vec phi, arma::vec sigma_v, arma::vec rho, arma::vec Ny1, double w, double eta, arma::vec sy1){
//   // Initializations
//   int T = y.n_rows;
//   arma::vec xiy1(T);
//   xiy1.zeros();
// 
//   arma::mat m_dat, s_dat;
//   double m_pri, s_pri, xi, p;
//   
//   arma::vec resid(4);
//   arma::mat sigma(4,4);
//   arma::mat sigma12;
//   arma::mat sigma22;
//   arma::mat sigma21;
//   arma::uvec idx = {1,2,3};
//   
//   for (int t=0; t<T; t++){
//     m_pri = w * sy1[t];
//     s_pri = eta * sqrt(sy1[t]);
//     xi = R::rnorm(m_pri,s_pri);
//     if (Ny1(t) == 0){
//       xiy1(t) = xi;
//     } else {
//       sigma = armgetSigma(arma::trans(omega.row(t)),  sigma_v,  rho);
//       sigma12 = sigma.cols(idx);
//       sigma12 = sigma12.row(0);
//       sigma22 = sigma.rows(idx);
//       sigma22 = sigma22.cols(idx);
//       sigma21 = sigma.rows(idx);
//       sigma21 = sigma21.col(0);
//       resid(0) = y(t,1) - (x(t,1) + mu(1));
//       resid(1) = omega(t+1,0) - (theta(0) + phi(0) * (omega(t,0) - theta(0)));
//       resid(2) = omega(t+1,1) - (theta(1) + phi(1) * (omega(t,1) - theta(1)));
//       resid(3) = y(t,0) - (x(t,0) + mu(0));
//       s_dat = sigma(0,0) - sigma12 * arma::inv(sigma22) * sigma21;
//       m_dat = resid(3) - sigma12 * arma::inv(sigma22) * resid.head(3);
//       p = 1 / s_dat(0,0) / (1 / s_dat(0,0) + 1 / s_pri);
//       xiy1(t) = p * R::rnorm(m_dat(0,0),sqrt(s_dat(0,0))) + (1 - p) * xi;
//     }
//   }
//   
//   return xiy1;
// }
// 
// // [[Rcpp::export]]
// arma::vec update_xiy2(arma::mat y, arma::mat x, arma::mat omega, arma::vec mu, arma::vec theta, arma::vec phi, arma::vec sigma_v, arma::vec rho, arma::vec Ny2, double w, double eta, arma::vec sy2){
//   // Initializations
//   int T = y.n_rows;
//   arma::vec xiy2(T);
//   xiy2.zeros();
//   
//   arma::mat m_dat, s_dat;
//   double m_pri, s_pri, xi, p;
//   
//   arma::vec resid(4);
//   arma::mat sigma(4,4);
//   arma::mat sigma12;
//   arma::mat sigma22;
//   arma::mat sigma21;
//   arma::uvec idx = {0,2,3};
//   
//   for (int t=0; t<T; t++){
//     m_pri = w * sy2[t];
//     s_pri = eta * sqrt(sy2[t]);
//     xi = R::rnorm(m_pri,s_pri);
//     if (Ny2(t) == 0){
//       xiy2(t) = xi;
//     } else {
//       sigma = armgetSigma(arma::trans(omega.row(t)),  sigma_v,  rho);
//       sigma12 = sigma.cols(idx);
//       sigma12 = sigma12.row(1);
//       sigma22 = sigma.rows(idx);
//       sigma22 = sigma22.cols(idx);
//       sigma21 = sigma.rows(idx);
//       sigma21 = sigma21.col(1);
//       resid(0) = y(t,0) - (x(t,0) + mu(0));
//       resid(1) = omega(t+1,0) - (theta(0) + phi(0) * (omega(t,0) - theta(0)));
//       resid(2) = omega(t+1,1) - (theta(1) + phi(1) * (omega(t,1) - theta(1)));
//       resid(3) = y(t,1) - (x(t,1) + mu(1));
//       s_dat = sigma(1,1) - sigma12 * arma::inv(sigma22) * sigma21;
//       m_dat = resid(3) - sigma12 * arma::inv(sigma22) * resid.head(3);
//       p = 1 / s_dat(0,0) / (1 / s_dat(0,0) + 1 / s_pri);
//       xiy2(t) = p * R::rnorm(m_dat(0,0),sqrt(s_dat(0,0))) + (1 - p) * xi;
//     }
//   }
//   
//   return xiy2;
// }
// 
// // [[Rcpp::export]]
// arma::mat update_xic(arma::mat y, arma::mat x, arma::mat omega, arma::vec mu, arma::vec theta, arma::vec phi, arma::vec sigma_v, arma::vec rho, arma::vec Nc, arma::vec wc, arma::vec sigmac, double rhoc, arma::vec sc){
//   // Initializations
//   int T = y.n_rows;
//   arma::mat xic(T,2);
//   xic.zeros();
//   arma::mat xi(1,2);
//   xi.zeros();
//   arma::mat sigma(4,4);
//   arma::vec resid(4);
//   
//   arma::mat m_dat, s_dat;
//   arma::vec m_pri(2);
//   arma::mat s_pri(2,2);
//   arma::uvec idx0 = {0,1};
//   arma::uvec idx1 = {2,3};
//   
//   arma::mat p(2,2);
//   arma::mat sigma11;
//   arma::mat sigma12;
//   arma::mat sigma22;
//   arma::mat sigma21;
//   
//   for (int t=0;t<T;t++){
//     m_pri(0) = wc(0) * sc(t);
//     m_pri(1) = wc(1) * sc(t);
//     s_pri(0, 0) = sigmac(0) * sigmac(0) * sc(t);
//     s_pri(1, 1) = sigmac(1) * sigmac(1) * sc(t);
//     s_pri(0, 1) = rhoc * sigmac(0) * sigmac(1) * sc(t);
//     s_pri(1, 0) = rhoc * sigmac(0) * sigmac(1) * sc(t);
//     
//     xi = mvrnormArma(1, m_pri, s_pri);
//     if (Nc(t) == 0){
//       xic.row(t) = xi;
//     } else {
//       sigma = armgetSigma(arma::trans(omega.row(t)),  sigma_v,  rho);
//       sigma11 = sigma.cols(idx0);
//       sigma11 = sigma11.rows(idx0);
//       sigma12 = sigma.cols(idx0);
//       sigma12 = sigma12.rows(idx1);
//       sigma22 = sigma.rows(idx1);
//       sigma22 = sigma22.cols(idx1);
//       sigma21 = sigma.rows(idx0);
//       sigma21 = sigma21.cols(idx1);
//       resid(0) = y(t,0) - (x(t,0) + mu(0));
//       resid(1) = y(t,1) - (x(t,1) + mu(1));
//       resid(2) = omega(t+1,0) - (theta(0) + phi(0) * (omega(t,0) - theta(0)));
//       resid(3) = omega(t+1,1) - (theta(1) + phi(1) * (omega(t,1) - theta(1)));
//       s_dat = sigma11 - sigma12 * arma::inv(sigma22) * sigma21;
//       m_dat = resid.head(2) - sigma12 * arma::inv(sigma22) * resid.tail(2);
//       
//       p = arma::inv(s_dat) * arma::inv(arma::inv(s_dat) + arma::inv(s_pri));
//       xic.row(t) = arma::trans(p * arma::trans(xi) + (1 - p) * arma::trans(mvrnormArma(1, m_dat, s_dat)));
//     }
//   }
//   
//   return xic;
// }

// [[Rcpp::export]]
arma::mat pgas_v_cpp(arma::mat y, arma::mat x, arma::mat omega, arma::mat J, arma::vec mu, arma::vec theta, arma::vec phi, arma::vec sigma_v, arma::vec rho, int N){
  // Initializations
  int T = y.n_rows;
  arma::cube v(T+1,2,N);
  v.zeros();
  arma::cube vfinal(T+1,2,N);
  vfinal.zeros();
  arma::vec vv(2);
  arma::mat sigma(4,4);
  arma::mat sigma_cond(2,2);
  arma::vec resid(4);
  arma::rowvec epsilon(2);

  arma::mat m_dat, s_dat;
  arma::uvec idx0 = {0,1};
  arma::uvec idx1 = {2,3};

  arma::mat sigma11;
  arma::mat sigma12;
  arma::mat sigma22;
  arma::mat sigma21;

  IntegerVector ancestors(N);

  NumericVector logancestorWeights(N);
  NumericVector logWeights(N);
  double maxlogW;
  NumericVector weights(N);
  NumericVector ancestorWeights(N);
  std::vector<int> Indices(N) ;
  std::iota (std::begin(Indices), std::end(Indices), 0);

  for (int i=0; i < N-1; i++){
    v(0,0,i) = R::rnorm(theta(0),sigma_v(0)/sqrt(1-pow(phi(0),2)));
    v(0,1,i) = R::rnorm(theta(1),sigma_v(1)/sqrt(1-pow(phi(1),2)));
  }
  v(0,0,N-1) = omega(0,0);
  v(0,1,N-1) = omega(0,1);
  // Weights update
  for (int i=0; i < N; i++){
    if ((v(0,0,i) < 0) | (v(0,1,i) < 0)){
      logWeights(i) = R_NegInf;
    } else {
      vv(0) = v(0,0,i);
      vv(1) = v(0,1,i);
      sigma = armgetSigma(vv, sigma_v, rho);
      sigma11 = sigma.rows(idx0);
      sigma11 = sigma11.cols(idx0);
      epsilon(0) = y(0,0) - (x(0,0) + mu(0) + J(0,0));
      epsilon(1) = y(0,1) - (x(0,1) + mu(1) + J(0,1));
      logWeights(i) = dmvnorm_arma(epsilon, arma::trans(arma::zeros(2)), sigma11, true)[0];
      sigma_cond(0,0) = sigma_v(0) * sigma_v(0) / (1 - phi(0) * phi(0));
      sigma_cond(1,1) = sigma_v(0) * sigma_v(0) / (1 - phi(0) * phi(0));
      sigma_cond(0,1) = rho(1) * sigma_v(0) * sigma_v(1) / sqrt((1 - phi(0) * phi(0)) * (1 - phi(1) * phi(1)));
      sigma_cond(1,0) = rho(1) * sigma_v(0) * sigma_v(1) / sqrt((1 - phi(0) * phi(0)) * (1 - phi(1) * phi(1)));
      epsilon(0) = v(0,0,i) - theta(0);
      epsilon(1) = v(0,1,i) - theta(1);
      logWeights(i) += dmvnorm_arma(epsilon, arma::trans(arma::zeros(2)), sigma_cond, true)[0];
      logWeights(i) += -R::dnorm(v(0,0,i),theta(0),sigma_v(0)/sqrt(1-pow(phi(0),2)),true);
      logWeights(i) += -R::dnorm(v(0,1,i),theta(1),sigma_v(1)/sqrt(1-pow(phi(1),2)),true);
    }
  }
  maxlogW = max(logWeights);
  weights = exp(logWeights - maxlogW);

  // t >= 2
  for (int t=1; t <= T; t++) {
    //for (t in c(2:T)){
    ancestors = lrf::csample_int(Indices,
                                 N,
                                 TRUE,
                                 weights);

    for (int i=0; i < N-1; i++){
      v(t,0,i) = R::rnorm(theta(0) + phi(0) * (v(t-1,0,ancestors(i)) - theta(0)) + rho(2) * sigma_v(0) * (y(t-1,0) - (x(t-1,0) + mu(0) + J(t-1,0))),
                          sigma_v(0) * sqrt(v(t-1,0,ancestors(i)) * (1 - rho(2) * rho(2))));
      v(t,1,i) = R::rnorm(theta(1) + phi(1) * (v(t-1,1,ancestors(i)) - theta(1)) + rho(3) * sigma_v(1) * (y(t-1,1) - (x(t-1,1) + mu(1) + J(t-1,1))),
                          sigma_v(1) * sqrt(v(t-1,1,ancestors(i)) * (1 - rho(3) * rho(3))));
    }
    v(t,0,N-1) = omega(t,0);
    v(t,1,N-1) = omega(t,1);

    // Ancestor sampling
    for (int i=0; i < N; i++){
      if ((v(t-1,0,i) < 0) | (v(t-1,1,i) < 0)){
        logancestorWeights(i) = R_NegInf;
      } else {
        vv(0) = v(t-1,0,i);
        vv(1) = v(t-1,1,i);
        resid(0) = y(t-1,0) - (x(t-1,0) + mu(0) + J(t-1,0));
        resid(1) = y(t-1,1) - (x(t-1,1) + mu(1) + J(t-1,1));
        resid(2) = omega(t,0) - (theta(0) + phi(0) * (vv(0) - theta(0)));
        resid(3) = omega(t,1) - (theta(1) + phi(1) * (vv(1) - theta(1)));
        sigma = armgetSigma(vv, sigma_v, rho);
        sigma11 = sigma.cols(idx1);
        sigma11 = sigma11.rows(idx1);
        sigma12 = sigma.cols(idx0);
        sigma12 = sigma12.rows(idx1);
        sigma22 = sigma.rows(idx0);
        sigma22 = sigma22.cols(idx0);
        sigma21 = sigma.rows(idx0);
        sigma21 = sigma21.cols(idx1);
        s_dat = sigma11 - sigma12 * arma::inv(sigma22) * sigma21;
        m_dat = resid.tail(2) - sigma12 * arma::inv(sigma22) * resid.head(2);
        epsilon(0) = m_dat(0,0);
        epsilon(1) = m_dat(1,0);
        logancestorWeights(i) = logWeights(i) + dmvnorm_arma(epsilon, arma::trans(arma::zeros(2)), s_dat, true)[0];
      }
    }
    ancestorWeights = exp(logancestorWeights - max(logancestorWeights));
    ancestors[N-1] = csample_int(Indices,
                                 1,
                                 TRUE,
                                 ancestorWeights)[0];

    for (int i=0;i < N; i++){
      //s[,1:(t-1)] = s[ancestors,1:(t-1)]
      vfinal.subcube(0,0,i,t-1,1,i) = v.subcube(0,0,ancestors[i],t-1,1,ancestors[i]);
    }
    vfinal.row(t) = v.row(t);
    v = vfinal;
    
    // Weights update
    for (int i=0; i < N; i++){
      if ((v(t,0,i) < 0) | (v(t,1,i) < 0)){
        logWeights(i) = R_NegInf;
      } else if (t == T){
        logWeights(i) = 0.0;
      } else {
        vv(0) = v(t-1,0,i);
        vv(1) = v(t-1,1,i);
        resid(0) = y(t-1,0) - (x(t-1,0) + mu(0) + J(t-1,0));
        resid(1) = y(t-1,1) - (x(t-1,1) + mu(1) + J(t-1,1));
        resid(2) = v(t,0,i) - (theta(0) + phi(0) * (vv(0) - theta(0)));
        resid(3) = v(t,1,i) - (theta(1) + phi(1) * (vv(1) - theta(1)));
        sigma = armgetSigma(vv, sigma_v, rho);
        sigma11 = sigma.cols(idx1);
        sigma11 = sigma11.rows(idx1);
        sigma12 = sigma.cols(idx0);
        sigma12 = sigma12.rows(idx1);
        sigma22 = sigma.rows(idx0);
        sigma22 = sigma22.cols(idx0);
        sigma21 = sigma.rows(idx0);
        sigma21 = sigma21.cols(idx1);
        s_dat = sigma11 - sigma12 * arma::inv(sigma22) * sigma21;
        m_dat = resid.tail(2) - sigma12 * arma::inv(sigma22) * resid.head(2);
        epsilon(0) = m_dat(0,0);
        epsilon(1) = m_dat(1,0);
        logWeights(i) = dmvnorm_arma(epsilon, arma::trans(arma::zeros(2)), s_dat, true)[0];
        vv(0) = v(t,0,i);
        vv(1) = v(t,1,i);
        sigma = armgetSigma(vv, sigma_v, rho);
        sigma11 = sigma.rows(idx0);
        sigma11 = sigma11.cols(idx0);
        epsilon(0) = y(t,0) - (x(t,0) + mu(0) + J(t,0));
        epsilon(1) = y(t,1) - (x(t,1) + mu(1) + J(t,1));
        logWeights(i) += dmvnorm_arma(epsilon, arma::trans(arma::zeros(2)), sigma11, true)[0];
        logWeights(i) += -R::dnorm(v(t,0,i),
                                  theta(0) + phi(0) * (v(t-1,0,i) - theta(0)) + rho(2) * sigma_v(0) * (y(t-1,0) - (x(t-1,0) + mu(0) + J(t-1,0))),
                                  sigma_v(0) * sqrt(v(t-1,0,i) * (1 - rho(2) * rho(2))),true);
        logWeights(i) += -R::dnorm(v(t,1,i),
                                   theta(1) + phi(1) * (v(t-1,1,i) - theta(1)) + rho(3) * sigma_v(1) * (y(t-1,1) - (x(t-1,1) + mu(1) + J(t-1,1))),
                                   sigma_v(1) * sqrt(v(t-1,1,i) * (1 - rho(3) * rho(3))),true);
      }
    }
    maxlogW = max(logWeights);
    weights = exp(logWeights - maxlogW);
  }
  //ind = sample(c(1:N),size = 1,prob =weights)
  double ind = csample_int(Indices,
                           1,
                           TRUE,
                           weights)[0];
  return v.slice(ind);
}

// [[Rcpp::export]]
arma::mat update_mu(arma::mat y, arma::mat x, arma::mat omega, arma::mat J, arma::vec theta, arma::vec phi, arma::vec sigma_v, arma::vec rho) {
  int T = y.n_rows;
  //N(0,1) prior
  arma::vec S(2); // prior mean over prior variance
  S.zeros();
  arma::mat W(2,2); // prior precision (inverse of prior variance)
  W(0,0) = 1;
  W(1,1) = 1;
  W(0,1) = 0;
  W(1,0) = 0;
  arma::mat sigma(4,4);
  arma::mat sigma_cond(2,2);
  arma::vec y_resid(2);
  arma::vec v_resid(2);
  arma::mat mu(2,1);
  arma::vec J1 = J.col(0);
  arma::vec J2 = J.col(1);
  for (int t = 0; t < T; t++) {
    sigma = armgetSigma(arma::trans(omega.row(t)),sigma_v,rho);
    sigma_cond = sigma.submat(0,0,1,1) - sigma.submat(0,2,1,3) * arma::inv(sigma.submat(2,2,3,3)) * sigma.submat(2,0,3,1);
    W += arma::inv(sigma_cond);
    y_resid(0) = y(t,0) - (x(t,0) + J1(t));
    y_resid(1) = y(t,1) - (x(t,1) + J2(t));
    v_resid(0) = omega(t+1,0) - (theta(0) + phi(0) * (omega(t,0) - theta(0)));
    v_resid(1) = omega(t+1,1) - (theta(1) + phi(1) * (omega(t,1) - theta(1)));
    S += arma::inv(sigma_cond) * (y_resid - sigma.submat(0,2,1,3) * arma::inv(sigma.submat(2,2,3,3)) * v_resid);
  }
  arma::vec m = arma::inv(W) * S;
  mu = mvrnormArma(1, m, arma::inv(W));
  return mu;
}

// [[Rcpp::export]]
arma::vec update_theta(arma::mat y, arma::mat x, arma::mat omega, arma::mat J, arma::vec mu, arma::vec theta, arma::vec phi, arma::vec sigma_v, arma::vec rho) {
  int T = y.n_rows;
  //N(0,10) prior
  double S = 0; // prior mean over prior variance
  double W = pow(10,-1); // inverse of prior variance
  arma::mat sigma(4,4);
  arma::mat sigma12;
  arma::mat sigma22;
  arma::mat sigma21;
  arma::mat mu_cond(1,1);
  arma::mat sigma_cond(1,1);
  arma::vec resid(4);
  arma::vec J1 = J.col(0);
  arma::vec J2 = J.col(1);
  arma::uvec idx0 = {0,1,3};
  arma::uvec idx1 = {0,1,2};
  for (int t = 0; t < T; t++) {
    sigma = armgetSigma(arma::trans(omega.row(t)),sigma_v,rho);
    sigma12 = sigma.cols(idx0);
    sigma12 = sigma12.row(2);
    sigma22 = sigma.rows(idx0);
    sigma22 = sigma22.cols(idx0);
    sigma21 = sigma.rows(idx0);
    sigma21 = sigma21.col(2);
    sigma_cond = sigma(2,2) - sigma12 * arma::inv(sigma22) * sigma21;
    W += pow(1 - phi(0),2) / sigma_cond(0,0);
    resid(0) = (y(t,0) - (x(t,0) + mu(0) + J1(t)));
    resid(1) = (y(t,1) - (x(t,1) + mu(1) + J2(t)));
    resid(2) = (omega(t+1,1) - (theta(1) + phi(1) * (omega(t,1) - theta(1))));
    resid(3) = (omega(t+1,0) - phi(0) * omega(t,0));
    mu_cond = (resid(3) - sigma12 * arma::inv(sigma22) * resid.head(3));
    S += (1 - phi(0)) * mu_cond(0,0) / sigma_cond(0,0);
  }
  theta(0) = RcppTN::rtn1(S / W, sqrt(1 / W), 0.0, INFINITY);
  S = 0;
  W = pow(10,-1);
  for (int t = 0; t < T; t++) {
    sigma = armgetSigma(arma::trans(omega.row(t)),sigma_v,rho);
    sigma12 = sigma.cols(idx1);
    sigma12 = sigma12.row(3);
    sigma22 = sigma.rows(idx1);
    sigma22 = sigma22.cols(idx1);
    sigma21 = sigma.rows(idx1);
    sigma21 = sigma21.col(3);
    sigma_cond = sigma(3,3) - sigma12 * arma::inv(sigma22) * sigma21;
    W += pow(1 - phi(1),2) / sigma_cond(0,0);
    resid(0) = (y(t,0) - (x(t,0) + mu(0) + J1(t)));
    resid(1) = (y(t,1) - (x(t,1) + mu(1) + J2(t)));
    resid(2) = (omega(t+1,0) - (theta(0) + phi(0) * (omega(t,0) - theta(0))));
    resid(3) = (omega(t+1,1) - phi(1) * omega(t,1));
    mu_cond = (resid(3) - sigma12 * arma::inv(sigma22) * resid.head(3));
    S += (1 - phi(1)) * mu_cond(0,0) / sigma_cond(0,0);
  }
  theta(1) = RcppTN::rtn1(S / W, sqrt(1 / W), 0.0, INFINITY);
  return theta;
}

// [[Rcpp::export]]
arma::vec update_phi(arma::mat y, arma::mat x, arma::mat omega, arma::mat J, arma::vec mu, arma::vec theta, arma::vec phi, arma::vec sigma_v, arma::vec rho) {
  int T = y.n_rows;
  //N(1,.25) prior
  arma::mat S(1,1); // prior mean over prior variance
  arma::mat W(1,1);
  S(0,0) = 1 / .25;
  W(0,0) = 1 / .25;
  double m, s2;
  arma::mat sigma(4,4);
  arma::mat sigma12;
  arma::mat sigma22;
  arma::mat sigma21;
  arma::mat sigma_cond(1,1);
  arma::vec resid(3);
  arma::vec J1 = J.col(0);
  arma::vec J2 = J.col(1);
  arma::uvec idx0 = {0,1,3};
  arma::uvec idx1 = {0,1,2};
  for (int t = 0; t < T; t++) {
    sigma = armgetSigma(arma::trans(omega.row(t)),sigma_v,rho);
    sigma12 = sigma.cols(idx0);
    sigma12 = sigma12.row(2);
    sigma22 = sigma.rows(idx0);
    sigma22 = sigma22.cols(idx0);
    sigma21 = sigma.rows(idx0);
    sigma21 = sigma21.col(2);
    sigma_cond = sigma(2,2) - sigma12 * arma::inv(sigma22) * sigma21;
    sigma_cond *= pow(omega(t,0) - theta(0),-2);
    W += arma::inv(sigma_cond);
    resid(0) = (y(t,0) - (x(t,0) + mu(0) + J1(t))) / (omega(t,0) - theta(0));
    resid(1) = (y(t,1) - (x(t,1) + mu(1) + J2(t))) / (omega(t,0) - theta(0));
    resid(2) = (omega(t+1,1) - (theta(1) + phi(1) * (omega(t,1) - theta(1)))) / (omega(t,0) - theta(0));
    S += arma::inv(sigma_cond) * ((omega(t+1,0) - theta(0)) / (omega(t,0) - theta(0)) - sigma12 * arma::inv(sigma22) * resid);
  }
  m = S(0,0) / W(0,0);
  s2 = 1 / W(0,0);
  phi(0) = RcppTN::rtn1(m, sqrt(s2), 0.0, 1.0);
  S(0,0) = 1 / .25;
  W(0,0) = 1 / .25;
  for (int t = 0; t < T; t++) {
    sigma = armgetSigma(arma::trans(omega.row(t)),sigma_v,rho);
    sigma12 = sigma.cols(idx1);
    sigma12 = sigma12.row(3);
    sigma22 = sigma.rows(idx1);
    sigma22 = sigma22.cols(idx1);
    sigma21 = sigma.rows(idx1);
    sigma21 = sigma21.col(3);
    sigma_cond = sigma(3,3) - sigma12 * arma::inv(sigma22) * sigma21;
    sigma_cond *= pow(omega(t,1) - theta(1),-2);
    W += arma::inv(sigma_cond);
    resid(0) = (y(t,0) - (x(t,0) + mu(0) + J1(t))) / (omega(t,1) - theta(1));
    resid(1) = (y(t,1) - (x(t,1) + mu(1) + J2(t))) / (omega(t,1) - theta(1));
    resid(2) = (omega(t+1,0) - (theta(0) + phi(0) * (omega(t,0) - theta(0)))) / (omega(t,1) - theta(1));
    S += arma::inv(sigma_cond) * ((omega(t+1,1) - theta(1)) / (omega(t,1) - theta(1)) - sigma12 * arma::inv(sigma22) * resid);
  }
  m = S(0,0) / W(0,0);
  s2 = 1 / W(0,0);
  phi(1) = RcppTN::rtn1(m, sqrt(s2), 0.0, 1.0);
  return phi;
}

// [[Rcpp::export]]
arma::vec update_sigma_v(arma::mat y, arma::mat x, arma::mat omega, arma::mat J, arma::vec mu, arma::vec theta, arma::vec phi, arma::vec sigma_v, arma::vec rho, double tune_sigmasq) {
  int T = y.n_rows;
  double part1_try=0; double part1_old=0;
  double part2_try=0; double part2_old=0;
  arma::mat sigma_try(4,4);
  arma::mat sigma_old(4,4);
  arma::vec proposal(2);
  double a;
  arma::rowvec resid(4);
  arma::vec J1 = J.col(0);
  arma::vec J2 = J.col(1);
  for (int k = 0; k < 2; k++){
    part1_try = 0;
    part1_old = 0;
    proposal = sigma_v;
    proposal(k) = sigma_v(k) + R::rnorm(0,tune_sigmasq);
    if (proposal(k) < 0){
      a = R_NegInf;
    }else{
      for (int t = 0; t < T; t++){
        resid(0) = y(t,0) - (x(t,0) + mu(0) + J1(t));
        resid(1) = y(t,1) - (x(t,1) + mu(1) + J2(t));
        resid(2) = omega(t+1,0) - (theta(0) + phi(0) * (omega(t,0) - theta(0)));
        resid(3) = omega(t+1,1) - (theta(1) + phi(1) * (omega(t,1) - theta(1)));
        sigma_try = armgetSigma(arma::trans(omega.row(t)), proposal, rho);
        sigma_old = armgetSigma(arma::trans(omega.row(t)), sigma_v, rho);
        part1_old += dmvnorm_arma(resid, arma::trans(arma::zeros(4)), sigma_old, true)[0];
        part1_try += dmvnorm_arma(resid, arma::trans(arma::zeros(4)), sigma_try, true)[0];
      }
      //normal(0,.1) prior
      part2_try = R::dnorm(proposal(k), 0, .1,true);//pdf_invgamma(proposal, 3, 2);//-log(M_PI*(1+pow(proposal, 2)));  
      part2_old = R::dnorm(sigma_v(k), 0, .1,true);//pdf_invgamma(sigma_v*sigma_v, 3, 2);//-log(M_PI*(1+sigmasq));
      a = part1_try + part2_try - part1_old - part2_old;
      //Rcout << exp(a) << std::endl;
    }
    if (a > log(R::runif(0,1)))  {
      sigma_v = (proposal);
    }
  }
  //Rcout << proposal << std::endl;
  return sigma_v;
}

// [[Rcpp::export]]
arma::vec update_rho(arma::mat y, arma::mat x, arma::mat omega, arma::mat J, arma::vec mu, arma::vec theta, arma::vec phi, arma::vec sigma_v, arma::vec rho, double tune_sigmasq) {
  int T = y.n_rows;
  double part1_try=0; double part1_old=0;
  double part2_try=0; double part2_old=0;
  arma::mat sigma_try(4,4);
  arma::mat sigma_old(4,4);
  arma::vec proposal(4);
  double a;
  arma::rowvec resid(4);
  arma::vec J1 = J.col(0);
  arma::vec J2 = J.col(1);
  for (int k = 0; k < 4; k++){
    part1_try = 0;
    part1_old = 0;
    proposal = rho;
    proposal(k) = rho(k) + R::rnorm(0,tune_sigmasq);
    if ((proposal(k) < -1) | (proposal(k) > 1)){
      a = R_NegInf;
    }else if (corDet(proposal) < 0){
      a = R_NegInf;
    } else {
      for (int t = 0; t < T; t++){
        resid(0) = y(t,0) - (x(t,0) + mu(0) + J1(t));
        resid(1) = y(t,1) - (x(t,1) + mu(1) + J2(t));
        resid(2) = omega(t+1,0) - (theta(0) + phi(0) * (omega(t,0) - theta(0)));
        resid(3) = omega(t+1,1) - (theta(1) + phi(1) * (omega(t,1) - theta(1)));
        sigma_try = armgetSigma(arma::trans(omega.row(t)), sigma_v, proposal);
        sigma_old = armgetSigma(arma::trans(omega.row(t)), sigma_v, rho);
        part1_old += dmvnorm_arma(resid, arma::trans(arma::zeros(4)), sigma_old, true)[0];
        part1_try += dmvnorm_arma(resid, arma::trans(arma::zeros(4)), sigma_try, true)[0];
      }
      //normal(0,.1) prior
      part2_try = 0;// Uniform(-1,1) prior 
      part2_old = 0;
      a = part1_try + part2_try - part1_old - part2_old;
      //Rcout << exp(a) << std::endl;
    }
    if (a > log(R::runif(0,1)))  {
      rho = (proposal);
    }
  }
  //Rcout << proposal << std::endl;
  return rho;
}

// [[Rcpp::export]]
double update_w(arma::vec xi, double w, double eta, arma::vec s, double tune_wsd){
  int T = xi.n_rows;
  double prior_try = 0;
  double prior_old = 0;
  double part1_try = 0;
  double part1_old = 0;
  double proposal;
  double a;
  proposal = R::rnorm(w,tune_wsd);
  // Cauchy(0,5) prior
  prior_try = R::dcauchy(proposal, 0, 5, true);
  prior_old = R::dcauchy(w, 0, 5, true);
  for (int j = 0; j < T; j++){
    part1_try += R::dnorm(xi(j), proposal * s(j), eta * sqrt(s(j)),true);
    part1_old += R::dnorm(xi(j), w * s(j), eta * sqrt(s(j)),true);
  }
  a = part1_try + prior_try - part1_old - prior_old;
  if (a > log(R::runif(0, 1))) {
    w = proposal;
  }
  return w;
}

// [[Rcpp::export]]
double update_eta(arma::vec xi, double w, double eta, arma::vec s, double tune_wsd){
  int T = xi.n_rows;
  double prior_try = 0;
  double prior_old = 0;
  double part1_try = 0;
  double part1_old = 0;
  double proposal;
  double a;
  proposal = R::rnorm(eta, tune_wsd);
  if (proposal < 0){
    a = R_NegInf;
  } else {
  // Cauchy(0,2.5) prior
  prior_try = R::dcauchy(proposal, 0, 5, true);
  prior_old = R::dcauchy(eta, 0, 5, true);
  for (int j = 0; j < T; j++){
    part1_try += R::dnorm(xi(j), w * s(j), proposal * sqrt(s(j)),true);
    part1_old += R::dnorm(xi(j), w * s(j), eta * sqrt(s(j)),true);
  }
  a = part1_try + prior_try - part1_old - prior_old;
  }
  if (a > log(R::runif(0, 1))) {
    eta = proposal;
  }
  return eta;
}

// [[Rcpp::export]]
arma::vec update_w_c(arma::mat xi_c, arma::vec w_c, arma::vec sigma_c, double rho_c, arma::vec s_c, double tune_wsd) {
  int T = xi_c.n_rows;
  double prior_try = 0;
  double prior_old = 0;
  double part1_try = 0;
  double part1_old = 0;
  double proposal;
  double a;
  proposal = R::rnorm(w_c(0),tune_wsd);
  // Cauchy(0,5) prior
  prior_try = R::dcauchy(proposal, 0, 5, true);
  prior_old = R::dcauchy(w_c(0), 0, 5, true);
  for (int j = 0; j < T; j++){
    part1_try += R::dnorm(xi_c(j,0), proposal * s_c(j) + rho_c * sigma_c(0) / sigma_c(1) * (xi_c(j,1) - w_c(1) * s_c(j)), sigma_c(0) * sqrt(s_c(j) * (1 - pow(rho_c,2))),true);
    part1_try += R::dnorm(xi_c(j,1), w_c(1) * s_c(j), sigma_c(1) * sqrt(s_c(j)),true);
    part1_old += R::dnorm(xi_c(j,0), w_c(0) * s_c(j) + rho_c * sigma_c(0) / sigma_c(1) * (xi_c(j,1) - w_c(1) * s_c(j)), sigma_c(0) * sqrt(s_c(j) * (1 - pow(rho_c,2))),true);
    part1_old += R::dnorm(xi_c(j,1), w_c(1) * s_c(j), sigma_c(1) * sqrt(s_c(j)),true);
  }
  a = part1_try + prior_try - part1_old - prior_old;
  if (a > log(R::runif(0, 1))) {
    w_c(0) = proposal;
  }
  
  part1_try = 0;
  part1_old = 0;
  proposal = R::rnorm(w_c(1),tune_wsd);
  // Cauchy(0,5) prior
  prior_try = R::dcauchy(proposal, 0, 5, true);
  prior_old = R::dcauchy(w_c(1), 0, 5, true);
  for (int j = 0; j < T; j++){
    part1_try += R::dnorm(xi_c(j,1), proposal * s_c(j) + rho_c * sigma_c(1) / sigma_c(0) * (xi_c(j,0) - w_c(0) * s_c(j)), sigma_c(1) * sqrt(s_c(j) * (1 - pow(rho_c,2))),true);
    part1_try += R::dnorm(xi_c(j,0), w_c(0) * s_c(j), sigma_c(0) * sqrt(s_c(j)),true);
    part1_old += R::dnorm(xi_c(j,1), w_c(1) * s_c(j) + rho_c * sigma_c(1) / sigma_c(0) * (xi_c(j,0) - w_c(0) * s_c(j)), sigma_c(1) * sqrt(s_c(j) * (1 - pow(rho_c,2))),true);
    part1_old += R::dnorm(xi_c(j,0), w_c(0) * s_c(j), sigma_c(0) * sqrt(s_c(j)),true);
  }
  a = part1_try + prior_try - part1_old - prior_old;
  if (a > log(R::runif(0, 1))) {
    w_c(1) = proposal;
  }
  return w_c;
}

// [[Rcpp::export]]
arma::vec update_sigma_c(arma::mat xi_c, arma::vec w_c, arma::vec sigma_c, double rho_c, arma::vec s_c, double tune_wsd) {
  int T = xi_c.n_rows;
  double prior_try = 0;
  double prior_old = 0;
  double part1_try = 0;
  double part1_old = 0;
  double proposal;
  double a;
  proposal = R::rnorm(sigma_c(0),tune_wsd);
  if (proposal < 0){
    a = R_NegInf;
  }else{
  // Cauchy(0,5) prior
  prior_try = R::dcauchy(proposal, 0, 5, true);
  prior_old = R::dcauchy(sigma_c(0), 0, 5, true);
  for (int j = 0; j < T; j++){
    part1_try += R::dnorm(xi_c(j,0), w_c(0) * s_c(j) + rho_c * proposal / sigma_c(1) * (xi_c(j,1) - w_c(1) * s_c(j)), proposal * sqrt(s_c(j) * (1 - pow(rho_c,2))),true);
    part1_try += R::dnorm(xi_c(j,1), w_c(1) * s_c(j), sigma_c(1) * sqrt(s_c(j)),true);
    part1_old += R::dnorm(xi_c(j,0), w_c(0) * s_c(j) + rho_c * sigma_c(0) / sigma_c(1) * (xi_c(j,1) - w_c(1) * s_c(j)), sigma_c(0) * sqrt(s_c(j) * (1 - pow(rho_c,2))),true);
    part1_old += R::dnorm(xi_c(j,1), w_c(1) * s_c(j), sigma_c(1) * sqrt(s_c(j)),true);
  }
  a = part1_try + prior_try - part1_old - prior_old;
  }
  if (a > log(R::runif(0, 1))) {
    sigma_c(0) = proposal;
  }
  
  part1_try = 0;
  part1_old = 0;
  proposal = R::rnorm(sigma_c(1),tune_wsd);
  if (proposal < 0){
    a = R_NegInf;
  }else{
  // Cauchy(0,5) prior
  prior_try = R::dcauchy(proposal, 0, 5, true);
  prior_old = R::dcauchy(w_c(1), 0, 5, true);
  for (int j = 0; j < T; j++){
    part1_try += R::dnorm(xi_c(j,1), w_c(1) * s_c(j) + rho_c * proposal / sigma_c(0) * (xi_c(j,0) - w_c(0) * s_c(j)), proposal * sqrt(s_c(j) * (1 - pow(rho_c,2))),true);
    part1_try += R::dnorm(xi_c(j,0), w_c(0) * s_c(j), sigma_c(0) * sqrt(s_c(j)),true);
    part1_old += R::dnorm(xi_c(j,1), w_c(1) * s_c(j) + rho_c * sigma_c(1) / sigma_c(0) * (xi_c(j,0) - w_c(0) * s_c(j)), sigma_c(1) * sqrt(s_c(j) * (1 - pow(rho_c,2))),true);
    part1_old += R::dnorm(xi_c(j,0), w_c(0) * s_c(j), sigma_c(0) * sqrt(s_c(j)),true);
  }
  a = part1_try + prior_try - part1_old - prior_old;
  }
  if (a > log(R::runif(0, 1))) {
    sigma_c(1) = proposal;
  }
  return sigma_c;
}

// [[Rcpp::export]]
double update_rho_c(arma::mat xi_c, arma::vec w_c, arma::vec sigma_c, double rho_c, arma::vec s_c, double tune_wsd) {
  int T = xi_c.n_rows;
  double prior_try = 0;
  double prior_old = 0;
  double part1_try = 0;
  double part1_old = 0;
  double proposal;
  double a;
  proposal = R::rnorm(rho_c, tune_wsd);
  if ((proposal < -1) | (proposal > 1)){
    a = R_NegInf;
  }else{
    // Uniform(-1,1) prior
    prior_try = 0;
    prior_old = 0;
    for (int j = 0; j < T; j++){
      part1_try += R::dnorm(xi_c(j,0), w_c(0) * s_c(j) + proposal * sigma_c(0) / sigma_c(1) * (xi_c(j,1) - w_c(1) * s_c(j)), sigma_c(0) * sqrt(s_c(j) * (1 - pow(proposal,2))),true);
      part1_try += R::dnorm(xi_c(j,1), w_c(1) * s_c(j), sigma_c(1) * sqrt(s_c(j)),true);
      part1_old += R::dnorm(xi_c(j,0), w_c(0) * s_c(j) + rho_c * sigma_c(0) / sigma_c(1) * (xi_c(j,1) - w_c(1) * s_c(j)), sigma_c(0) * sqrt(s_c(j) * (1 - pow(rho_c,2))),true);
      part1_old += R::dnorm(xi_c(j,1), w_c(1) * s_c(j), sigma_c(1) * sqrt(s_c(j)),true);
    }
    a = part1_try + prior_try - part1_old - prior_old;
  }
  if (a > log(R::runif(0, 1))) {
    rho_c = proposal;
  }
  return rho_c;
}

// [[Rcpp::export]]
arma::vec update_lambda(arma::vec sumN, arma::vec prior_probs){
   double a0 = sumN(0) + prior_probs(0);
   double a1 = sumN(1) + prior_probs(1);
   double a2 = sumN(2) + prior_probs(2);
   double a3 = sumN(3) + prior_probs(3);
 
   arma::vec lambda(4);
   lambda(0) = R::rgamma(a0, 1);
   lambda(1) = R::rgamma(a1, 1);
   lambda(2) = R::rgamma(a2, 1);
   lambda(3) = R::rgamma(a3, 1);
 
   return lambda/sum(lambda);
}


