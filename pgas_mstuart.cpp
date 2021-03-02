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

arma::mat armgetSigma(double omegat, double sigma_v, double rho){
  arma::mat Sigma(2,2);
  Sigma(0,0) = omegat;
  Sigma(1,1) = sigma_v*sigma_v*omegat;
  Sigma(0,1) = rho*sigma_v*omegat;
  Sigma(1,0) = rho*sigma_v*omegat;
  return Sigma;
}

// [[Rcpp::export]]
double pdf_norm(double x, double mean, double sd) {

    return R::dnorm(x, 0, 2,false);

}



std::vector<int> csample_int( std::vector<int> x, 
                              int size,
                              bool replace, 
                              NumericVector prob) {
  std::vector<int> ret = RcppArmadillo::sample(x, size, replace, prob) ;
  return ret ;
}

// [[Rcpp::export]]
IntegerVector update_delta(arma::vec y, arma::vec x, arma::vec omega, arma::vec xiy, double mu, double theta, double phi, double sigma_v ,double rho, arma::vec lambda){
  int T = y.n_rows;
  IntegerVector N(T);
  double part1_N1=0; double part1_N0=0;
  
  arma::vec diffy_jump(T);
  arma::vec diffomega_jump(T);  
  
  arma::vec diffy_jumpy(T);
  arma::vec diffomega_jumpy(T);
  
  arma::rowvec epsilon(2);
  arma::mat sigma(2,2);
  
  double log_P1; double log_P0;
  
  std::vector<int> N_possible(2) ; 
  std::iota (std::begin(N_possible), std::end(N_possible), 0); 
  //N_possible[2] = 3;
  NumericVector probs(2);
  //N equal 0;  jumps in just y
  diffomega_jumpy = omega.subvec(1,T) - (theta+phi*(omega.subvec(0,T-1)-theta));
  diffy_jumpy = y - x - mu - xiy;
  
  //N equal 1; no jumps
  diffomega_jump =omega.subvec(1,T) - (theta+phi*(omega.subvec(0,T-1)-theta));
  diffy_jump = y -x-mu ;
  
  for (int j = 0; j < T; j++){
    sigma = armgetSigma(omega[j],  sigma_v,  rho);
    
    epsilon[0] = diffy_jumpy(j);
    epsilon[1] = diffomega_jumpy(j);
    part1_N0 = dmvnorm_arma(epsilon,  arma::trans(arma::zeros(2)),  sigma, true)[0];//-.5*diff[0]*pow(proposal,-2);
  
    epsilon[0] = diffy_jump(j);
    epsilon[1] = diffomega_jump(j);
    part1_N1 = dmvnorm_arma(epsilon,  arma::trans(arma::zeros(2)),  sigma, true)[0];//-.5*diff[0]*pow(gamma2[i],-1); 
    
    
    //Rcout << part1_N2 << std::endl;
    log_P0 = part1_N0 + log(lambda[0]);
    log_P1 = part1_N1 + log(lambda[1]);   
    
    // for numerical stability
    double pmx = max(NumericVector::create(log_P1, log_P0));
    log_P1 = log_P1 - pmx;
    log_P0 = log_P0 - pmx;

    //Rcout << log_P0 << " "<<log_P1 << " " << log_P2 <<std::endl;
    probs[0] = exp(log_P0)/(exp(log_P0) + exp(log_P1)); // y
    probs[1] = exp(log_P1)/(exp(log_P0) + exp(log_P1)); // none
    //Rcout << probs << std::endl;
    N[j] = csample_int(N_possible, 
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
    weights = weights/sum(weights);

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

// Ancestor sampling: Weights are the same
//    ancestorWeights = weights;
//    ancestors[N-1] = csample_int(Indices,
//                        1,
//                        TRUE,
//                        ancestorWeights)[0];

    // Weights update
    for (int i=0; i < N; i++){
      logWeights[i] = R::dnorm(xi[t], w * s(i, t), eta * sqrt(s(i, t)), true);
    }
    double maxlogW = max(logWeights);
    weights = exp(logWeights - maxlogW);
    weights = weights/sum(weights);
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

// // [[Rcpp::export]]
// arma::rowvec pgas_sc_cpp(arma::mat xi, arma::rowvec w_c, arma::mat Sigma_c, arma::vec sprim, int N){
//   // Initializations
//   int T = xi.n_rows;
//   arma::mat s(N,T);
//   s.zeros();
//   arma::mat sfinal(N,T);
//   s.zeros();
//   arma::mat sigma_temp(2,2);
//   sigma_temp.zeros();
// 
//   IntegerVector ancestors(N);
// 
//   NumericVector logWeights(N);
//   double maxlogW;
//   NumericVector weights(N);
//   NumericVector ancestorWeights(N);
//   std::vector<int> Indices(N) ;
//   std::iota (std::begin(Indices), std::end(Indices), 0);
// 
//   // Sample v using CSMC
//   // t = 2 (first time point for the xi variables)
//   for (int i=0; i < N-1; i++){
//     s(i,0) = R::rexp(1);//#draw from transition
//   }
//   //#draw from transition
//   s(N-1,0) = sprim[0];
//   // Weights update
//   //logWeights = logG_s(s.col(0), xi[0], w, eta,N);
//   for (int i=0; i < N; i++){
//     sigma_temp(1,1) = Sigma_c(1,1)*s(i,0);
//     sigma_temp(0,0) = Sigma_c(0,0)*s(i,0);
//     sigma_temp(1,0) = Sigma_c(1,0)*s(i,0);
//     sigma_temp(0,1) = Sigma_c(0,1)*s(i,0);
//     logWeights[i] = dmvnorm_arma(xi.row(0),  s(i,0)*w_c,  sigma_temp, true)[0];
//   }
//   maxlogW = max(logWeights);
//   weights = exp(logWeights - maxlogW);
//   weights = weights/sum(weights);
// 
//   // t >= 2
//   for (int t=1; t < T; t++) {
//     //for (t in c(2:T)){
//     ancestors = lrf::csample_int(Indices,
//                                  N,
//                                  TRUE,
//                                  weights);
// 
//     //s[1:N-1,t] = rexp(N-1, rate = 1)//#draw from transition
//     //s.submat(0,t,N-2,t)=R::rexp(N-1);
//     for (int i=0; i < N-1; i++){
//       s(i,t) = R::rexp(1);//#draw from transition
//     }
//     s(N-1,t) = sprim[t];
// 
//     // Ancestor sampling: Weights are the same
//     //    ancestorWeights = weights;
//     //    ancestors[N-1] = csample_int(Indices,
//     //                        1,
//     //                        TRUE,
//     //                        ancestorWeights)[0];
// 
//     // Weights update
//     for (int i=0; i < N; i++){
//       sigma_temp(1,1) = Sigma_c(1,1)*s(i,t);
//       sigma_temp(0,0) = Sigma_c(0,0)*s(i,t);
//       sigma_temp(1,0) = Sigma_c(1,0)*s(i,t);
//       sigma_temp(0,1) = Sigma_c(0,1)*s(i,t);
//       logWeights[i] = dmvnorm_arma(xi.row(t),  s(i,t)*w_c,  sigma_temp, true)[0];
//     }
//     maxlogW = max(logWeights);
//     weights = exp(logWeights - maxlogW);
//     weights = weights/sum(weights);
//     for (int i=0;i < N; i++){
//       //s[,1:(t-1)] = s[ancestors,1:(t-1)]
//       sfinal.submat(i,0,i,t-1) = s.submat(ancestors[i], 0, ancestors[i], t-1);
//     }
//     sfinal.col(t) = s.col(t);
//     s = sfinal;
// 
//   }
//   //ind = sample(c(1:N),size = 1,prob =weights)
//   double ind = csample_int(Indices,
//                            1,
//                            TRUE,
//                            weights)[0];
//   return s.row(ind);
// 
// }

// [[Rcpp::export]]
arma::rowvec pgas_xiy_cpp(arma::vec y, arma::vec x, arma::vec omega, double mu, double theta, double phi, double sigma_v, double rho, arma::vec xiyprim,  arma::vec Ny, double w, double eta, arma::vec sy, int N){
  // Initializations
  int T = xiyprim.n_rows;
  arma::mat xiy(N,T);
  xiy.zeros();
  arma::mat xiyfinal(N,T);
  xiyfinal.zeros();


  double m_dat, m_pri, s_dat, s_pri;
  double p = 0.5;
  double Z;
  double dat_ll;
  double pri_ll;

  IntegerVector ancestors(N);

  NumericVector logWeights(N);
  double maxlogW;
  NumericVector weights(N);
  NumericVector ancestorWeights(N);
  std::vector<int> Indices(N) ;
  std::iota (std::begin(Indices), std::end(Indices), 0);

  for (int i=0; i < N-1; i++){
    m_dat = y[0] - (x[0] + mu + rho/sigma_v*(omega[1] - (theta+phi*(omega[0]-theta))));
    m_pri = w * sy[0];
    s_dat = sqrt(omega[0] * (1-rho*rho));
    s_pri = eta * sqrt(sy[0]);
    Z = R::runif(0,1);
    if (Z > p){
      xiy(i,0) = R::rnorm(m_dat, s_dat);//#draw from transition
    } else {
      xiy(i,0) = R::rnorm(m_pri, s_pri);//#draw from transition
    }

  }
  xiy(N-1,0) = xiyprim[0];
  // Weights update
  for (int i=0; i < N; i++){
    dat_ll = R::dnorm(xiy(i,0), m_dat, s_dat, true);
    pri_ll = R::dnorm(xiy(i,0), m_pri, s_pri, true);
    logWeights(i) = Ny[0] * dat_ll + pri_ll - std::max(dat_ll,pri_ll) - log(p + (1 - p) * exp(-std::abs(pri_ll - dat_ll)));
  }
  maxlogW = max(logWeights);
  weights = exp(logWeights - maxlogW);
  weights = weights/sum(weights);

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
      m_dat = y[t] - (x[t] + mu + rho/sigma_v*(omega[t+1] - (theta+phi*(omega[t]-theta))));
      m_pri = w * sy[t];
      s_dat = sqrt(omega[t] * (1 - rho * rho));
      s_pri = eta * sqrt(sy[t]);
      Z = R::runif(0,1);
      if (Z > p){
        xiy(i,t) = R::rnorm(m_dat, s_dat);//#draw from transition
      } else {
        xiy(i,t) = R::rnorm(m_pri, s_pri);//#draw from transition
      }
    }
    xiy(N-1,t) = xiyprim[t];

    // Ancestor sampling: Weights are the same
    //    ancestorWeights = weights;
    //    ancestors[N-1] = csample_int(Indices,
    //                        1,
    //                        TRUE,
    //                        ancestorWeights)[0];

    // Weights update
    for (int i=0; i < N; i++){
      dat_ll = R::dnorm(xiy(i,t), m_dat, s_dat, true);
      pri_ll = R::dnorm(xiy(i,t), m_pri, s_pri, true);
      logWeights(i) = Ny[t] * dat_ll + pri_ll - std::max(dat_ll,pri_ll) - log(p + (1 - p) * exp(-std::abs(pri_ll - dat_ll)));
    }
    maxlogW = max(logWeights);
    weights = exp(logWeights - maxlogW);
    weights = weights/sum(weights);
    for (int i=0;i < N; i++){
      //s[,1:(t-1)] = s[ancestors,1:(t-1)]
      xiyfinal.submat(i,0,i,t-1) = xiy.submat(ancestors[i], 0, ancestors[i], t-1);
    }
    xiyfinal.col(t) = xiy.col(t);
    xiy = xiyfinal;

  }
  //ind = sample(c(1:N),size = 1,prob =weights)
  double ind = csample_int(Indices,
                           1,
                           TRUE,
                           weights)[0];
  return xiy.row(ind);

}

// // [[Rcpp::export]]
// arma::rowvec pgas_xiv_cpp(arma::vec y, arma::vec x, arma::vec omega, double mu, double theta, double phi, double sigma_v, double rho, arma::vec xiy, arma::vec xivprim, arma::mat xic, arma::vec Ny, arma::vec Nv, arma::vec Nc, double w, double eta, arma::vec sv, int N){
//   // Initializations
//   int T = xivprim.n_rows;
//   arma::mat xiv(N,T);
//   xiv.zeros();
//   arma::mat xivfinal(N,T);
//   xivfinal.zeros();
// 
//   double m_dat, m_pri, s_dat, s_pri;
//   double p = 0.5;
//   double Z;
//   double dat_ll;
//   double pri_ll;
// 
//   IntegerVector ancestors(N);
// 
//   NumericVector logWeights(N);
//   double maxlogW;
//   NumericVector weights(N);
//   NumericVector ancestorWeights(N);
//   std::vector<int> Indices(N) ;
//   std::iota (std::begin(Indices), std::end(Indices), 0);
// 
//   for (int i=0; i < N-1; i++){
//     m_dat = omega[1] - (theta + phi * (omega[0] - theta) + Nc[0] * xic(0,1) + rho * sigma_v * (y[0] - (x[0] + mu + Ny[0] * xiy[0] + Nc[0] * xic(0,1))));
//     m_pri = w * sv[0];
//     s_dat = sigma_v * sqrt(omega[0] * (1 - rho * rho));
//     s_pri = eta * sqrt(sv[0]);
//     Z = R::runif(0,1);
//     if (Z > p){
//       xiv(i,0) = R::rnorm(m_dat, s_dat);//#draw from transition
//     } else {
//       xiv(i,0) = R::rnorm(m_pri, s_pri);//#draw from transition
//     }
//   }
//   xiv(N-1,0) = xivprim[0];
//   // Weights update
//   for (int i=0; i < N; i++){
//     dat_ll = R::dnorm(xiv(i,0), m_dat, s_dat, true);
//     pri_ll = R::dnorm(xiv(i,0), m_pri, s_pri, true);
//     logWeights(i) = Nv[0] * dat_ll + pri_ll - std::max(dat_ll,pri_ll) - log(p + (1 - p) * exp(-std::abs(pri_ll - dat_ll)));
//   }
//   maxlogW = max(logWeights);
//   weights = exp(logWeights - maxlogW);
//   weights = weights/sum(weights);
// 
//   // t >= 2
//   for (int t=1; t < T; t++) {
//     //for (t in c(2:T)){
//     ancestors = lrf::csample_int(Indices,
//                                  N,
//                                  TRUE,
//                                  weights);
// 
//     //s[1:N-1,t] = rexp(N-1, rate = 1)//#draw from transition
//     //s.submat(0,t,N-2,t)=R::rexp(N-1);
//     for (int i=0; i < N-1; i++){
//       m_dat = omega[t+1] - (theta + phi * (omega[t] - theta) + Nc[t] * xic(t,1) + rho * sigma_v * (y[t] - (x[t] + mu + Ny[t] * xiy[t] + Nc[t] * xic(t,1))));
//       m_pri = w * sv[t];
//       s_dat = sigma_v * sqrt(omega[t] * (1 - rho * rho));
//       s_pri = eta * sqrt(sv[t]);
//       Z = R::runif(0,1);
//       if (Z > p){
//         xiv(i,t) = R::rnorm(m_dat, s_dat);//#draw from transition
//       } else {
//         xiv(i,t) = R::rnorm(m_pri, s_pri);//#draw from transition
//       }
//     }
//     xiv(N-1,t) = xivprim[t];
// 
//     // Ancestor sampling
//     ancestorWeights = weights;
//     ancestors[N-1] = csample_int(Indices,
//                                  1,
//                                  TRUE,
//                                  ancestorWeights)[0];
// 
//     // Weights update
//     for (int i=0; i < N; i++){
//       dat_ll = R::dnorm(xiv(i,t), m_dat, s_dat, true);
//       pri_ll = R::dnorm(xiv(i,t), m_pri, s_pri, true);
//       logWeights(i) = Nv[t] * dat_ll + pri_ll - std::max(dat_ll,pri_ll) - log(p + (1 - p) * exp(-std::abs(pri_ll - dat_ll)));
//     }
//     maxlogW = max(logWeights);
//     weights = exp(logWeights - maxlogW);
//     weights = weights/sum(weights);
//     for (int i=0;i < N; i++){
//       //s[,1:(t-1)] = s[ancestors,1:(t-1)]
//       xivfinal.submat(i,0,i,t-1) = xiv.submat(ancestors[i], 0, ancestors[i], t-1);
//     }
//     xivfinal.col(t) = xiv.col(t);
//     xiv = xivfinal;
// 
//   }
//   //ind = sample(c(1:N),size = 1,prob =weights)
//   double ind = csample_int(Indices,
//                            1,
//                            TRUE,
//                            weights)[0];
//   return xiv.row(ind);
// }
// 
// // [[Rcpp::export]]
// arma::mat pgas_xic_cpp(arma::vec y, arma::vec x, arma::vec omega, double mu, double theta, double phi, double sigma_v, double rho, arma::vec xiy, arma::vec xiv, arma::mat xicprim, arma::vec Ny, arma::vec Nv, arma::vec Nc,  arma::vec wc, arma::vec sigmac, double rhoc, arma::vec sc, int N){
//   // Initializations
//   int T = xicprim.n_rows;
//   arma::cube xic(T,2,N);
//   xic.zeros();
//   arma::cube xicfinal(T,2,N);
//   xicfinal.zeros();
//   arma::mat xi(1,2);
//   xi.zeros();
// 
//   double Z;
//   double p = 0.5;
//   double dat_ll;
//   double pri_ll;
//   arma::vec Mu_dat(2);
//   Mu_dat.zeros();
//   arma::mat Sigma_dat(2,2);
//   Sigma_dat.zeros();
//   arma::vec Mu_pri(2);
//   Mu_pri.zeros();
//   arma::mat Sigma_pri(2,2);
//   Sigma_pri.zeros();
// 
//   arma::rowvec epsilon(2);
// 
//   IntegerVector ancestors(N);
// 
//   NumericVector logWeights(N);
//   double maxlogW;
//   NumericVector weights(N);
//   NumericVector ancestorWeights(N);
//   std::vector<int> Indices(N) ;
//   std::iota (std::begin(Indices), std::end(Indices), 0);
// 
//   Mu_dat(0) = y[0] - (x[0] + mu + Ny[0] + xiy[0]);
//   Mu_dat(1) = omega[1] - (theta + phi * (omega[0] - theta) + Nv[0] * xiv[0]);
//   //std::cout << "Mu_dat = " << Mu_dat << std::endl;
// 
//   Mu_pri(0) = wc(0) * sc(0);
//   Mu_pri(1) = wc(1) * sc(0);
//   //std::cout << "Mu_pri = " << Mu_pri << std::endl;
// 
//   Sigma_dat = armgetSigma(omega[0], sigma_v, rho);
//   //std::cout << "Sigma_dat = " << Sigma_dat << std::endl;
// 
//   Sigma_pri(0, 0) = sigmac(0) * sigmac(0) * sc(0);
//   Sigma_pri(1, 1) = sigmac(1) * sigmac(1) * sc(0);
//   Sigma_pri(0, 1) = rhoc * sigmac(0) * sigmac(1) * sc(0);
//   Sigma_pri(1, 0) = rhoc * sigmac(0) * sigmac(1) * sc(0);
//   //std::cout << "Sigma_pri = " << Sigma_pri << std::endl;
//   for (int i=0; i < N-1; i++){
//     Z = R::runif(0,1);
//     if (Z > p){
//       xi = mvrnormArma(1, Mu_pri, Sigma_pri);
//     } else {
//       xi = mvrnormArma(1, Mu_dat, Sigma_dat);
//     }
//     xic(0,0,i) = xi(0,0);
//     xic(0,1,i) = xi(0,1);
//   }
//   xic(0,0,N-1) = xicprim(0,0);
//   xic(0,1,N-1) = xicprim(0,1);
//   // Weights update
//   for (int i=0; i < N; i++){
//     epsilon(0) = xic(0,0,i);
//     epsilon(1) = xic(0,1,i);
//     dat_ll = dmvnorm_arma(epsilon, arma::trans(Mu_dat), Sigma_dat, true)[0];
//     //std::cout << "dat_ll = " << dat_ll << std::endl;
//     pri_ll = dmvnorm_arma(epsilon, arma::trans(Mu_pri), Sigma_pri, true)[0];
//     //std::cout << "pri_ll = " << pri_ll << std::endl;
//     logWeights(i) = Nc[0] * dat_ll + pri_ll - std::max(dat_ll,pri_ll) - log(p + (1 - p) * exp(-std::abs(pri_ll - dat_ll)));
//   }
//   maxlogW = max(logWeights);
//   weights = exp(logWeights - maxlogW);
//   weights = weights/sum(weights);
// 
//   // t >= 2
//   for (int t=1; t < T; t++) {
//     //for (t in c(2:T)){
//     ancestors = lrf::csample_int(Indices,
//                                  N,
//                                  TRUE,
//                                  weights);
// 
//     Mu_dat(0) = y[t] - (x[t] + mu + Ny[t] + xiy[t]);
//     Mu_dat(1) = omega[t+1] - (theta + phi * (omega[t] - theta) + Nv[t] * xiv[t]);
// 
//     Mu_pri(0) = wc(0) * sc(t);
//     Mu_pri(1) = wc(1) * sc(t);
// 
//     Sigma_dat = armgetSigma(omega[t], sigma_v, rho);
// 
//     Sigma_pri(0, 0) = sigmac(0) * sigmac(0) * sc(t);
//     Sigma_pri(1, 1) = sigmac(1) * sigmac(1) * sc(t);
//     Sigma_pri(0, 1) = rhoc * sigmac(0) * sigmac(1) * sc(t);
//     Sigma_pri(1, 0) = rhoc * sigmac(0) * sigmac(1) * sc(t);
//     for (int i=0; i < N-1; i++){
//       Z = R::runif(0,1);
//       if (Z > p){
//         xi = mvrnormArma(1, Mu_pri, Sigma_pri);
//       } else {
//         xi = mvrnormArma(1, Mu_dat, Sigma_dat);
//       }
//       xic(t,0,i) = xi(0,0);
//       xic(t,1,i) = xi(0,1);
//     }
//     xic(t,0,N-1) = xicprim(t,0);
//     xic(t,1,N-1) = xicprim(t,1);
// 
//     // Ancestor sampling
//     ancestorWeights = weights;
//     ancestors[N-1] = csample_int(Indices,
//                                  1,
//                                  TRUE,
//                                  ancestorWeights)[0];
// 
//     // Weights update
//     for (int i=0; i < N; i++){
//       epsilon(0) = xic(t,0,i);
//       epsilon(1) = xic(t,1,i);
//       dat_ll = dmvnorm_arma(epsilon, arma::trans(Mu_dat), Sigma_dat, true)[0];
//       //std::cout << "dat_ll = " << dat_ll << std::endl;
//       pri_ll = dmvnorm_arma(epsilon, arma::trans(Mu_pri), Sigma_pri, true)[0];
//       //std::cout << "pri_ll = " << pri_ll << std::endl;
//       logWeights(i) = Nc[t] * dat_ll + pri_ll - std::max(dat_ll,pri_ll) - log(p + (1 - p) * exp(-std::abs(pri_ll - dat_ll)));
//     }
//     maxlogW = max(logWeights);
//     weights = exp(logWeights - maxlogW);
//     weights = weights/sum(weights);
//     for (int i=0;i < N; i++){
//       //s[,1:(t-1)] = s[ancestors,1:(t-1)]
//       xicfinal.subcube(0,0,i,t-1,1,i) = xic.subcube(0,0,ancestors[i],t-1,1,ancestors[i]);
//       xicfinal.subcube(t,0,i,t,1,i) = xic.subcube(t,0,i,t,1,i);
//     }
//     xic = xicfinal;
//   }
//   //ind = sample(c(1:N),size = 1,prob =weights)
//   double ind = csample_int(Indices,
//                            1,
//                            TRUE,
//                            weights)[0];
//   return xic.slice(ind);
// }

// [[Rcpp::export]]
arma::rowvec pgas_v_cpp(arma::vec y, arma::vec x, arma::vec omega, arma::vec J, double mu, double theta, double phi, double sigma_v, double rho, int N){
  int T = y.n_rows;
  arma::mat v(N,T+1);
  v.zeros();
  arma::mat vfinal(N,T+1);
  vfinal.zeros();
  
  arma::rowvec epsilon(2);
  IntegerVector ancestors(N);
  
  NumericVector logWeights(N);
  double maxlogW;
  //double tmp;
  NumericVector weights(N);
  NumericVector logancestorWeights(N);
  NumericVector ancestorWeights(N);
  std::vector<int> Indices(N) ;
  std::iota (std::begin(Indices), std::end(Indices), 0);
  
  for (int i = 0; i < N - 1; i++){
    v(i,0) = RcppTN::rtn1(theta, sigma_v/sqrt(1-phi*phi),0,1000000);
  }
  v(N-1,0) = omega(0);
  for (int i = 0; i < N; i++){
    if (v(i,0) < 0){
      logWeights(i) = R_NegInf;
    } else {
      logWeights(i) = R::dnorm(y(0), x(0) + mu + J(0), sqrt(v(i,0)), true);
    }
  }
  //std::cout << v.col(0) << std::endl;
  //std::cout << logWeights << std::endl;
  maxlogW = max(logWeights);
  weights = exp(logWeights - maxlogW);
  weights = weights/sum(weights);
  
  for (int t = 1; t <= T; t++){
    // Simulate Ancestors
    ancestors = lrf::csample_int(Indices,
                                 N,
                                 TRUE,
                                 weights);
    // Simulate New Particles
    for (int i = 0; i < N-1; i++){
      v(i,t) = RcppTN::rtn1(theta + phi*(v(ancestors(i),t-1) - theta) + rho * sigma_v * (y(t-1) - (x(t-1) + mu + J(t-1))),
                        sigma_v * sqrt(v(ancestors(i),t-1) * (1 - rho * rho)),0,1000000);
    }
    v(N-1,t) = omega(t);
    //std::cout << v.col(t) << std::endl;
    // Simulate Final Ancestor
    for (int i = 0; i < N; i++){
      if (v(i,t-1) < 0){
        logancestorWeights(i) = R_NegInf;
      } else {
        logancestorWeights(i) = logWeights(i) + R::dnorm(omega(t),
                           theta + phi*(v(i,t-1) - theta) + rho * sigma_v * (y(t-1) - (x(t-1) + mu + J(t-1))),
                           sigma_v * sqrt(v(i,t-1) * (1 - rho * rho)),true);
      }
    }
    ancestorWeights = exp(logancestorWeights - max(logancestorWeights));
    ancestors[N-1] = csample_int(Indices,
                                 1,
                                 TRUE,
                                 ancestorWeights)[0];
    // Calculate New Weights
    for (int i = 0; i < N; i++){
      if (v(i,t) < 0){
        logWeights(i) = R_NegInf;
      } else if (t==T){
        logWeights(i) = 0;
      } else {
        logWeights(i) = R::dnorm(y(t), x(t) + mu + J(t), sqrt(v(i,t)), true);
      }
    }
    //std::cout << logWeights << std::endl;
    maxlogW = max(logWeights);
    weights = exp(logWeights - maxlogW);
    weights = weights/sum(weights);
    for (int i=0;i < N; i++){
      //s[,1:(t-1)] = s[ancestors,1:(t-1)]
      vfinal.submat(i,0,i,t-1) = v.submat(ancestors[i], 0, ancestors[i], t-1);
    }
    vfinal.col(t) = v.col(t);
    v = vfinal;
  }
  //ind = sample(c(1:N),size = 1,prob =weights)
  double ind = csample_int(Indices,
                           1,
                           TRUE,
                           weights)[0];
  return v.row(ind);
}

// [[Rcpp::export]]
double update_mu(arma::vec y, arma::vec x, arma::vec omega, arma::vec J, double theta, double phi, double sigma_v, double rho) {
  int T = y.n_rows;
  //N(0,1) prior
  double S = 0.0; // prior mean over prior variance
  double W = pow(1.0, -2); // prior precision (inverse of prior variance)
  double y_resid, v_resid;
  double mu;
  for (int j = 0; j < T; j++) {
    W += 1.0 / omega(j) / (1 - pow(rho, 2.0));
    y_resid = y(j) - (x(j) + J(j));
    v_resid = omega(j + 1) - (theta + phi * (omega(j) - theta));
    S += 1.0 / omega(j) / (1 - pow(rho, 2.0)) * (y_resid - rho / sigma_v * v_resid);
  }
  mu = R::rnorm(S / W, sqrt(1 / W));
  return mu;
}

// [[Rcpp::export]]
double update_theta(arma::vec y, arma::vec x, arma::vec omega, arma::vec J, double mu, double phi, double sigma_v, double rho) {
  int T = y.n_rows;
  //truncated (positive) N(0,10) prior
  double S = 0.0; // prior mean over prior variance
  double W = pow(10.0, -1.0); // prior precision (inverse of prior variance)
  double y_resid, v_resid;
  double theta;
  for (int j = 0; j < T; j++) {
    W += pow(1.0 - phi, 2.0) /(omega(j) * pow(sigma_v, 2.0) * (1 - pow(rho, 2.0)));
    y_resid = y(j) - (x(j) + mu + J(j));
    v_resid = omega(j + 1) - (phi * omega(j));
    S += (1.0 - phi)*(v_resid - rho * sigma_v * y_resid) / (omega(j) *pow(sigma_v, 2.0) * (1 - pow(rho, 2.0))) ;
  }
  //Rcout << S/W <<std::endl;
  //Rcout << 1/W <<std::endl;
  theta = RcppTN::rtn1(S / W, sqrt(1 / W), 0.0, INFINITY);
  return theta;
}

// [[Rcpp::export]]
double update_phi(arma::vec y, arma::vec x, arma::vec omega, arma::vec J, double mu, double theta, double sigma_v, double rho){
  int T = y.n_rows;
  //truncated (positive) N(0,1) prior
  double S = 1.0/pow(.25, 2); // prior mean over prior variance 
  double W = pow(.25, -2.0); // prior precision (inve rse of prior variance)
  double y_resid, v_resid;
  double phi;
  for (int j = 0; j < T; j++){
    W += pow(pow(sigma_v, 2.0)*(1 - pow(rho, 2.0)),-1)*pow(theta - omega(j), 2.0) / omega(j)  ;
    y_resid = y(j) - (x(j) + mu + J(j));
    v_resid =omega(j + 1) - theta;//omega(j + 1) - (omega(j) + J2(j));
    S += (omega(j)-theta)* (v_resid - rho * sigma_v * y_resid) / (omega(j)* (pow(sigma_v, 2.0) * (1 - pow(rho, 2.0)))) ; // theat - omega(j)
  }
  //Rcout << S/W <<std::endl;
  //Rcout << 1/W <<std::endl;
  phi = RcppTN::rtn1(S / W, sqrt(1 / W), 0.0, 1.0);
  return phi;
}

// [[Rcpp::export]]
double update_sigma_v(arma::vec y, arma::vec x, arma::vec omega, arma::vec J, double mu, double theta, double phi, double sigma_v, double rho, double tune_sigmasq){
  int T = y.n_rows;
  double part1_try=0; double part1_old=0;
  double part2_try=0; double part2_old=0;
  double proposal;
  double a;
  arma::rowvec resid(2);
  double y_resid;
  double v_resid;
  proposal = sigma_v + R::rnorm(0,tune_sigmasq);
  if (proposal < 0){
    a = R_NegInf;
  }else{
  for (int t = 0; t < T; t++){
      y_resid = y(t) - (x(t) + mu +J(t)) ;//diffy(j-1);
      v_resid = omega(t+1) - (theta + phi * (omega(t) - theta));
      //Rcout << v_resid << std::endl;
      //Rcout << epsilon <<std::endl;
      part1_old += R::dnorm(v_resid,  rho*sigma_v*y_resid,  sigma_v*sqrt(omega(t)*(1-rho*rho)), true);//-.5*diff[0]*pow(proposal,-2);
      part1_try += R::dnorm(v_resid,  rho*proposal*y_resid,  proposal*sqrt(omega(t)*(1-rho*rho)), true);//-.5*diff[0]*pow(gamma2[i],-1);
    }
    //normal(0,.1) prior
    part2_try = R::dnorm(proposal, 0, .1,true);//pdf_invgamma(proposal, 3, 2);//-log(M_PI*(1+pow(proposal, 2)));
    part2_old = R::dnorm(sigma_v, 0, .1,true);//pdf_invgamma(sigma_v*sigma_v, 3, 2);//-log(M_PI*(1+sigmasq));
    a = part1_try + part2_try - part1_old - part2_old;
    //Rcout << exp(a) << std::endl;
  }
  if (a > log(R::runif(0,1)))  {
    sigma_v = (proposal);
  }
  //Rcout << proposal << std::endl;
  return sigma_v;
}

// [[Rcpp::export]]
double update_rho(arma::vec y, arma::vec x, arma::vec omega, arma::vec J, double mu, double theta, double phi, double sigma_v, double rho, double tune_rhosd) {
  int T = y.n_rows;
  double part1_try = 0;
  double part1_old = 0;
  double part2_try = 0;
  double part2_old = 0;
  double prop_try = 0;
  double prop_old = 0;

  double proposal;
  double a;

  proposal = R::rnorm(rho,tune_rhosd);
  if ((proposal < -1) | (proposal > 1)){
    a = R_NegInf;
  }else{
    for (int t = 0; t < T; t++) {
      part1_try += R::dnorm(omega(t+1),theta + phi * (omega(t) - theta),sigma_v * sqrt(omega(t)),true);
      part1_try += R::dnorm(y(t),x(t) + mu + J(t) + proposal / sigma_v * (omega(t+1) - (theta + phi * (omega(t) - theta))),sqrt(omega(t) * (1 - proposal * proposal)),true);
      part1_old += R::dnorm(omega(t+1),theta + phi * (omega(t) - theta),sigma_v * sqrt(omega(t)),true);
      part1_old += R::dnorm(y(t),x(t) + mu + J(t) + rho / sigma_v * (omega(t+1) - (theta + phi * (omega(t) - theta))),sqrt(omega(t) * (1 - rho * rho)),true);
    }

    part2_try = 0; // log(1) (unif(0,1) density)
    part2_old = 0;
    a = part1_try + part2_try - prop_try - part1_old - part2_old + prop_old;
  }
  if (a > log(R::runif(0, 1))) {
    rho = proposal;
  }

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

// // [[Rcpp::export]]
// arma::vec update_w_c(arma::mat xi_c, arma::vec w_c, arma::vec sigma_c, double rho_c, arma::vec s_c, double tune_wsd) {
//   int T = xi_c.n_rows;
//   double prior_try = 0;
//   double prior_old = 0;
//   double part1_try = 0;
//   double part1_old = 0;
//   double proposal;
//   double a;
//   proposal = R::rnorm(w_c(0),tune_wsd);
//   // Cauchy(0,5) prior
//   prior_try = R::dcauchy(proposal, 0, 5, true);
//   prior_old = R::dcauchy(w_c(0), 0, 5, true);
//   for (int j = 0; j < T; j++){
//     part1_try += R::dnorm(xi_c(j,0), proposal * s_c(j) + rho_c * sigma_c(0) / sigma_c(1) * (xi_c(j,1) - w_c(1) * s_c(j)), sigma_c(0) * sqrt(s_c(j) * (1 - pow(rho_c,2))),true);
//     part1_try += R::dnorm(xi_c(j,1), w_c(1) * s_c(j), sigma_c(1) * sqrt(s_c(j)),true);
//     part1_old += R::dnorm(xi_c(j,0), w_c(0) * s_c(j) + rho_c * sigma_c(0) / sigma_c(1) * (xi_c(j,1) - w_c(1) * s_c(j)), sigma_c(0) * sqrt(s_c(j) * (1 - pow(rho_c,2))),true);
//     part1_old += R::dnorm(xi_c(j,1), w_c(1) * s_c(j), sigma_c(1) * sqrt(s_c(j)),true);
//   }
//   a = part1_try + prior_try - part1_old - prior_old;
//   if (a > log(R::runif(0, 1))) {
//     w_c(0) = proposal;
//   }
// 
//   part1_try = 0;
//   part1_old = 0;
//   proposal = R::rnorm(w_c(1),tune_wsd);
//   // Cauchy(0,5) prior
//   prior_try = R::dcauchy(proposal, 0, 5, true);
//   prior_old = R::dcauchy(w_c(1), 0, 5, true);
//   for (int j = 0; j < T; j++){
//     part1_try += R::dnorm(xi_c(j,1), proposal * s_c(j) + rho_c * sigma_c(1) / sigma_c(0) * (xi_c(j,0) - w_c(0) * s_c(j)), sigma_c(1) * sqrt(s_c(j) * (1 - pow(rho_c,2))),true);
//     part1_try += R::dnorm(xi_c(j,0), w_c(0) * s_c(j), sigma_c(0) * sqrt(s_c(j)),true);
//     part1_old += R::dnorm(xi_c(j,1), w_c(1) * s_c(j) + rho_c * sigma_c(1) / sigma_c(0) * (xi_c(j,0) - w_c(0) * s_c(j)), sigma_c(1) * sqrt(s_c(j) * (1 - pow(rho_c,2))),true);
//     part1_old += R::dnorm(xi_c(j,0), w_c(0) * s_c(j), sigma_c(0) * sqrt(s_c(j)),true);
//   }
//   a = part1_try + prior_try - part1_old - prior_old;
//   if (a > log(R::runif(0, 1))) {
//     w_c(1) = proposal;
//   }
//   return w_c;
// }
// 
// // [[Rcpp::export]]
// arma::vec update_sigma_c(arma::mat xi_c, arma::vec w_c, arma::vec sigma_c, double rho_c, arma::vec s_c, double tune_wsd) {
//   int T = xi_c.n_rows;
//   double prior_try = 0;
//   double prior_old = 0;
//   double part1_try = 0;
//   double part1_old = 0;
//   double proposal;
//   double a;
//   proposal = R::rnorm(sigma_c(0),tune_wsd);
//   if (proposal < 0){
//     a = R_NegInf;
//   }else{
//   // Cauchy(0,5) prior
//   prior_try = R::dcauchy(proposal, 0, 5, true);
//   prior_old = R::dcauchy(sigma_c(0), 0, 5, true);
//   for (int j = 0; j < T; j++){
//     part1_try += R::dnorm(xi_c(j,0), w_c(0) * s_c(j) + rho_c * proposal / sigma_c(1) * (xi_c(j,1) - w_c(1) * s_c(j)), proposal * sqrt(s_c(j) * (1 - pow(rho_c,2))),true);
//     part1_try += R::dnorm(xi_c(j,1), w_c(1) * s_c(j), sigma_c(1) * sqrt(s_c(j)),true);
//     part1_old += R::dnorm(xi_c(j,0), w_c(0) * s_c(j) + rho_c * sigma_c(0) / sigma_c(1) * (xi_c(j,1) - w_c(1) * s_c(j)), sigma_c(0) * sqrt(s_c(j) * (1 - pow(rho_c,2))),true);
//     part1_old += R::dnorm(xi_c(j,1), w_c(1) * s_c(j), sigma_c(1) * sqrt(s_c(j)),true);
//   }
//   a = part1_try + prior_try - part1_old - prior_old;
//   }
//   if (a > log(R::runif(0, 1))) {
//     sigma_c(0) = proposal;
//   }
// 
//   part1_try = 0;
//   part1_old = 0;
//   proposal = R::rnorm(sigma_c(1),tune_wsd);
//   if (proposal < 0){
//     a = R_NegInf;
//   }else{
//   // Cauchy(0,5) prior
//   prior_try = R::dcauchy(proposal, 0, 5, true);
//   prior_old = R::dcauchy(w_c(1), 0, 5, true);
//   for (int j = 0; j < T; j++){
//     part1_try += R::dnorm(xi_c(j,1), w_c(1) * s_c(j) + rho_c * proposal / sigma_c(0) * (xi_c(j,0) - w_c(0) * s_c(j)), proposal * sqrt(s_c(j) * (1 - pow(rho_c,2))),true);
//     part1_try += R::dnorm(xi_c(j,0), w_c(0) * s_c(j), sigma_c(0) * sqrt(s_c(j)),true);
//     part1_old += R::dnorm(xi_c(j,1), w_c(1) * s_c(j) + rho_c * sigma_c(1) / sigma_c(0) * (xi_c(j,0) - w_c(0) * s_c(j)), sigma_c(1) * sqrt(s_c(j) * (1 - pow(rho_c,2))),true);
//     part1_old += R::dnorm(xi_c(j,0), w_c(0) * s_c(j), sigma_c(0) * sqrt(s_c(j)),true);
//   }
//   a = part1_try + prior_try - part1_old - prior_old;
//   }
//   if (a > log(R::runif(0, 1))) {
//     sigma_c(1) = proposal;
//   }
//   return sigma_c;
// }
// 
// // [[Rcpp::export]]
// double update_rho_c(arma::mat xi_c, arma::vec w_c, arma::vec sigma_c, double rho_c, arma::vec s_c, double tune_wsd) {
//   int T = xi_c.n_rows;
//   double prior_try = 0;
//   double prior_old = 0;
//   double part1_try = 0;
//   double part1_old = 0;
//   double proposal;
//   double a;
//   proposal = R::rnorm(rho_c, tune_wsd);
//   if ((proposal < -1) | (proposal > 1)){
//     a = R_NegInf;
//   }else{
//     // Uniform(-1,1) prior
//     prior_try = 0;
//     prior_old = 0;
//     for (int j = 0; j < T; j++){
//       part1_try += R::dnorm(xi_c(j,0), w_c(0) * s_c(j) + proposal * sigma_c(0) / sigma_c(1) * (xi_c(j,1) - w_c(1) * s_c(j)), sigma_c(0) * sqrt(s_c(j) * (1 - pow(proposal,2))),true);
//       part1_try += R::dnorm(xi_c(j,1), w_c(1) * s_c(j), sigma_c(1) * sqrt(s_c(j)),true);
//       part1_old += R::dnorm(xi_c(j,0), w_c(0) * s_c(j) + rho_c * sigma_c(0) / sigma_c(1) * (xi_c(j,1) - w_c(1) * s_c(j)), sigma_c(0) * sqrt(s_c(j) * (1 - pow(rho_c,2))),true);
//       part1_old += R::dnorm(xi_c(j,1), w_c(1) * s_c(j), sigma_c(1) * sqrt(s_c(j)),true);
//     }
//     a = part1_try + prior_try - part1_old - prior_old;
//   }
//   if (a > log(R::runif(0, 1))) {
//     rho_c = proposal;
//   }
//   return rho_c;
// }

// [[Rcpp::export]]
arma::vec update_lambda(arma::vec sumN, arma::vec prior_probs){
   double a0 = sumN(0) + prior_probs(0);
   double a1 = sumN(1) + prior_probs(1);

   arma::vec lambda(2);
   lambda(0) = R::rgamma(a0, 1);
   lambda(1) = R::rgamma(a1, 1);
   
   return lambda/sum(lambda);
}

// [[Rcpp::export]]
double pg_lnp_y_mid_theta(arma::vec y, arma::vec x, double mu, double theta, double phi, double sigma_v, double rho, double xi_yw, double xi_yeta, double lambda, int N){
  int T = y.n_rows;
  double lnp_y_mid_theta = 0;
  std::vector<int> N_possible(2) ; 
  std::iota (std::begin(N_possible), std::end(N_possible), 0); 
  arma::mat v(N,T);
  arma::mat xi_y(N,T);
  arma::mat delta(N,T);
  arma::mat xi_ys(N,T);
  
  arma::mat v_final(N,T);
  arma::mat xi_y_final(N,T);
  arma::mat delta_final(N,T);
  arma::mat xi_ys_final(N,T);
  
  NumericVector probs(2);
  probs(0) = 1 - lambda;
  probs(1) = lambda;
  
  IntegerVector ancestors(N);
  
  //double maxlogW;
  //double sumW;
  NumericVector logWeights(N);
  NumericVector weights(N);
  std::vector<int> Indices(N);
  std::iota (std::begin(Indices), std::end(Indices), 0);
  
  for (int i=0; i < N;i++){
    v(i,0) = RcppTN::rtn1(theta,sigma_v/sqrt(1-phi*phi),0,1000000);
    xi_ys(i,0) = R::rexp(1);
    xi_y(i,0) = R::rnorm(xi_yw*xi_ys(i,0),xi_yeta*sqrt(xi_ys(i,0)));
    delta(i,0) = csample_int(N_possible, 
          1,
          TRUE, 
          probs)[0];
    if (v(i,0) < 0){
      logWeights(i) = R_NegInf;
    } else {
      logWeights(i) = R::dnorm(y(0), x(0) + mu + delta(i,0)*xi_y(i,0), sqrt(v(i,0)), true);
    }
  }
  weights = exp(logWeights);
  lnp_y_mid_theta += log(mean(weights));
  std::cout << lnp_y_mid_theta << std::endl;
  weights = weights/sum(weights);
  
  for (int t=1; t < T;t++){
    std::cout << t << std::endl;
    ancestors = lrf::csample_int(Indices,
                                 N,
                                 TRUE,
                                 weights);
    for (int i=0; i < N; i++){
      v(i,t) = RcppTN::rtn1(theta + phi*(v(ancestors(i),t-1) - theta) + rho * sigma_v * (y(t-1) - (x(t-1) + mu + delta(ancestors(i),t-1)*xi_y(ancestors(i),t-1))),
                        sigma_v * sqrt(v(ancestors(i),t-1) * (1 - rho * rho)),0,1000000);
      xi_ys(i,t) = R::rexp(1);
      xi_y(i,t) = R::rnorm(xi_yw*xi_ys(i,t),xi_yeta*sqrt(xi_ys(i,t)));
      delta(i,t) = csample_int(N_possible, 
            1,
            TRUE, 
            probs)[0];
      if (v(i,t) < 0){
        logWeights(i) = R_NegInf;
      } else {
        logWeights(i) = R::dnorm(y(t), x(t) + mu + delta(i,t)*xi_y(i,t), sqrt(v(i,t)), true);
      }
      v_final.submat(i,0,i,t-1) = v.submat(ancestors(i), 0, ancestors(i), t-1);
      xi_y_final.submat(i,0,i,t-1) = xi_y.submat(ancestors(i), 0, ancestors(i), t-1);
      delta_final.submat(i,0,i,t-1) = delta.submat(ancestors(i), 0, ancestors(i), t-1);
    }
    v_final.col(t) = v.col(t);
    xi_y_final.col(t) = xi_y.col(t);
    delta_final.col(t) = delta.col(t);
    v = v_final;
    xi_y = xi_y_final;
    delta = delta_final;
    
    weights = exp(logWeights);
    lnp_y_mid_theta += log(mean(weights));
    std::cout << lnp_y_mid_theta << std::endl;
    weights = weights/sum(weights);
  }
  weights = exp(logWeights);
  return lnp_y_mid_theta;
}


