
#information criterion calculations FOR 1 D ONLY
##---> RDIC and DIC
#where RDIC = D(theta-hat) + 2tr(I(theta-hat)*VCOV(theta-hat))
#and   DIC7 = -4E(ln(p(y|z,theta))) + 2ln(p(y|z-hat, theta-hat))
#and   DIC1 = -4E(ln(p(y|theta))) + 2ln(p(y| theta-hat))
library(ggplot2)
library(grid)
library(gridExtra)
library(MASS)
library(tmvtnorm)
library(ald)
library(ggplot2)
library(dplyr)
library(truncnorm)
library(mvtnorm)
library(tmvtnorm)
library(Rcpp)
library(MCMCpack)
library(quantmod)
library(RcppTN)

# DATA AND MCMC SAMPLES ------------
#download full MCMC results (obtained from run_mcmc.R with fix = FALSE)
keepsBTC <- readRDS("keepsBTC.rds")
keepsSP <- readRDS("keepsSP.rds")
#need to re-run MCMC with parameters fixed at MAP estimates (find code in empirical_data_study.R)
#this is needed for D(theta-hat), Information matrix

#download full MCMC results (obtained from run_mcmc.R with fix = TRUE)
keepsBTC_MAP <- readRDS("keepsBTC_MAP.rds")
keepsSP_MAP  <- readRDS("keepsSP_MAP.rds")

#download data 
# getSymbols("BTC-USD",from = "2014-09-17",to = "2020-09-30")
# BTC <- as.data.frame(`BTC-USD`)
# BTC$Date <- rownames(BTC)
# getSymbols("^GSPC",from = "2014-09-17",to = "2020-09-30")
# SP500 <- as.data.frame(`GSPC`)
# SP500$Date <- rownames(SP500)
# S <- merge(BTC,SP500)
# T <- nrow(S) - 1
# y <- data.frame(
# BTC = 100*log(S$`GSPC.Close`),
# SP = 100*log(S$`GSPC.Close`)
# )

#-----collect MAP estimates---------
#(used for additional MCMC run to get p(Y|theta-hat) & 1st, 2nd derivatives)

#1d model
THETA_map1d <- list()

THETA_map1d$lambda<-  cbind(keepsBTC$lambda,keepsSP$lambda) %>%  apply(2, mean) 
THETA_map1d$sigma_v <- cbind(keepsBTC$sigma_v,keepsSP$sigma_v) %>% apply(2, mean)
THETA_map1d$phi <-cbind(keepsBTC$phi,keepsSP$phi) %>%apply(2, mean)
THETA_map1d$theta<- cbind(keepsBTC$theta,keepsSP$theta) %>%apply(2, mean)
THETA_map1d$mu <- cbind(keepsBTC$mu,keepsSP$mu) %>%apply(2, mean)
THETA_map1d$rho<- cbind(keepsBTC$rho,keepsSP$rho) %>%apply(2, mean)
THETA_map1d$xi_yeta <- cbind(keepsBTC$xi_yeta,keepsSP$xi_yeta) %>%apply(2, mean)
THETA_map1d$xi_yw<- cbind(keepsBTC$xi_yw,keepsSP$xi_yw) %>%apply(2, mean)
THETA_map1d$v <- cbind(apply(keepsBTC$v,2,mean),apply(keepsSP$v,2,mean))
THETA_map1d$J <- cbind(apply(keepsBTC$J,2,mean),apply(keepsSP$J,2,mean))


#RDIC CALCULATIONS  (1)------- D(theta-hat) = 2ln(p(y|theta-hat)) = 2*(Q-H)----

#ln(p(y|theta-hat)) approximated using weights created in PGAS (volatility)

Elnpy_mid_theta <- keepsBTC$lnp_y_mid_theta + keepsSP$lnp_y_mid_theta #(not sure what they are named)

#ln(p(y|theta-hat)) approximated using weights created in PGAS (volatility)

Elnpy_mid_ztheta <- keepsBTC_MAP$lnp_y_mid_ztheta + keepsSP_MAP$lnp_y_mid_ztheta
  
#RDIC CALCULATIONS  (2)------- 2tr(I(theta-hat)*V(theta-hat))----

#-----get VCOV matrices

#THIS NEEDS TO BE IN THE SAME ORDER AS THE DERIVATIVES MATRIX!
SPsamps <- cbind(as.data.frame(keepsSP[c("lambda", "sigma_v", "phi", "theta", "mu", "rho", "xi_yw", "xi_yeta")]) %>%mutate(kappa = 1-phi) )
colnames(SPsamps) <- paste0(colnames(SPsamps), "SP")

my_order <- c("rho","rhoSP", "mu","muSP", "theta","thetaSP", "kappa","kappaSP", "sigma_v","sigma_vSP",
              "xi_yw","xi_ywSP", "xi_yeta", "xi_yetaSP", "lambda", "lambdaSP")


VCOV_1d <- keepsBTC[c("lambda", "sigma_v", "phi", "theta", "mu", "rho", "xi_yw", "xi_yeta")] %>%
  as.data.frame() %>%
  mutate(kappa = 1-phi) %>%
  cbind(SPsamps) %>%
  as.data.frame() %>%
  dplyr::select(my_order)%>% 
  cov()


#-----get Information matrices
#the Information matrix components are calculated WRT p(z|theta-hat, Y)
#i.e., don't use posterior means of z - take posterior means of deritavies
#using the MCMC samps where theta is fixed at theta-hat
#these use 

#---FIRST DERIVATIVE FUNCTIONS for I(theta-hat)-----
#MCMC = list of MCMC (particle) draws where theta HAS BEEN HELD FIXED AT MAPS
#MAPs = list of MAP values at which theta was held (e.g., THETA_map1d, THETA_map2d)
#functions work on ONE ITERATION r AT A TIME

#many of the functions need y, v residuals
#parameters (e.g., mu, theta, etc...) = MAP estimate
#latent variables (e.g., J, v, etc...) = r^th iteration's conditional posterior sample
#y_idx = which y are we looking at (SP = 2 or BTC = 1)

get_y_resid1d <- function(y_idx, MAPs, MCMCs, r){
 y[,y_idx] - ( MAPs$mu[y_idx] + MCMCs$J[r,] ) #output Tx1 (vector)
}

get_y_resid2d <- function(MAPs,MCMCs, r){
  y - (MAPs$mu + MCMCs$J[r,]) #output Tx2 matrix
}
 
get_v_resid1d <- function(y_idx, MAPs,MCMCs, r){
  MCMCs$v[r,-1] - (MAPs$theta[y_idx] + MAPs$phi[y_idx]*(MCMCs$v[r, -T] - MAPs$theta[y_idx]))#output Tx1 (vector)
}

get_v_resid2d <- function(MAPs, MCMCs,r){
  MCMCs$v[r,-1,] - (MAPs$theta + MAPs$phi*(MCMCs$v[r,-T,] - MAPs$theta))#output Tx2 matrix
}



#RHO 1
dL_drho_1d <- function(y_idx, MAPs,MCMCs,r){
  rho <- MAPs$rho[y_idx]
  sigma_v <- MAPs$sigma_v[y_idx]
  epsilon_y <- get_y_resid1d(y_idx, MAPs, MCMCs,r)
  epsilon_v <- get_v_resid1d( MAPs, MCMCs,r)
  v <- MCMCs$v[r,]
  #derivative return
  T*rho/((1-rho^2)) - 
    sum(rho*epsilon_y^2/(v[-T]*(1-rho^2)) - 
          (1+ rho^2)*epsilon_y*epsilon_v/(sigma_v*v[-T]*(1-rho^2)^2)  +
          rho*epsilon_v^2/(sigma_v^2 * v[-T]*(1-rho^2)^2) )
}

#MU 1
dL_dmu_1d <- function(y_idx, MAPs,MCMCs,r){
  rho <- MAPs$rho[y_idx]
  sigma_v <- MAPs$sigma_v[y_idx]
  epsilon_y <- get_y_resid1d(y_idx, MAPs, MCMCs,r)
  epsilon_v <- get_v_resid1d(MAPs, MCMCs,r)
  v <- MCMCs$v[r,]
  #derivative return
    -sum(-epsilon_y/(v[-T]*(1-rho^2)) + 
           rho*epsilon_v/(v[-T]*(1-rho^2)*sigma_v) )
}
            
#THETA 1
dL_dtheta_1d <- function(y_idx, MAPs, MCMCs, r){
  rho <- MAPs$rho[y_idx]
  sigma_v <- MAPs$sigma_v[y_idx]
  kappa <- 1-MAPs$phi[y_idx]
  epsilon_y <- get_y_resid1d(y_idx, MAPs, MCMCs,r)
  epsilon_v <- get_v_resid1d(y_idx,MAPs, MCMCs,r)
  #derivative return  
  -sum(-kappa*epsilon_v/(v[-T]*(1-rho^2)*sigma_v^2) + 
         rho*kappa*epsilon_y/(v[-T]*(1-rho^2)*sigma_v))

  }

#KAPPA 1
dL_dkappa_1d <- function(y_idx, MAPs, MCMCs, r){
  rho <- MAPs$rho[y_idx]
  sigma_v <- MAPs$sigma_v[y_idx]
  kappa <- 1-MAPs$phi[y_idx]
  theta <- MAPs$theta[y_idx]
  epsilon_y <- get_y_resid1d(y_idx, MAPs, MCMCs,r)
  epsilon_v <- get_v_resid1d(y_idx,MAPs, MCMCs,r)
  #derivative return  
  -sum(-(theta - v[-T])*epsilon_v/(sigma_v^2*v[-T]*(1-rho^2)) +
         rho*(theta - v[-T])*epsilon_y/(sigma_v*v[-T]*(1-rho^2)))
  
}

#SIGMA_V 1
dL_dsigma_v_1d <- function(y_idx, MAPs, MCMCs, r){
  rho <- MAPs$rho[y_idx]
  sigma_v <- MAPs$sigma_v[y_idx]
  kappa <- 1-MAPs$phi[y_idx]
  theta <- MAPs$theta[y_idx]
  epsilon_y <- get_y_resid1d(y_idx, MAPs, MCMCs,r)
  epsilon_v <- get_v_resid1d(y_idx,MAPs, MCMCs,r)
  #derivative return  
  -T/sigma_v - 
    sum(rho*epsilon_y*epsilon_v/(sigma_v*sigma_v*v[-T]*(1-rho*rho)) - 
            epsilon_v*epsilon_v/(sigma_v*sigma_v*sigma_v*v[-T]*(1-rho*rho)))
  
}

#MU_Y (xi_yw)
dL_dmuy_1d <- function(y_idx, MAPs, MCMCs, r){
  B <- MCMCs$xi_ys[r,] #NEED TO SAVE XI_YS IN run_mcmc.R
  mu_y <- MAPs$xi_yw[y_idx]
  v_y <- MAPs$xi_yeta[y_idx]
  xi_y <- MCMCs$xi_y #NEED TO SAVE XI_Y IN run_mcmc.R
  
  sum((xi_y - mu_y*B)/(v_y*v_y))
  
}

#MU_Y (xi_yw)
dL_dvy_1d <- function(y_idx, MAPs, MCMCs, r){
  B <- MCMCs$xi_ys[r,] #NEED TO SAVE XI_YS IN run_mcmc.R
  mu_y <- MAPs$xi_yw[y_idx]
  v_y <- MAPs$xi_yeta[y_idx]
  xi_y <- MCMCs$xi_y #NEED TO SAVE XI_Y IN run_mcmc.R
  
  -T/v_y +sum((xi_y - mu_y*B)^2/(v_y*v_y*v_y*B))
  
}

#LAMBDA
dL_dlambda_1d <- function(y_idx, MAPs, MCMCs, r){
  lambda <- MAPs$lambda[y_idx]
  N <- MCMCs$delta[r,]
  sum(N)/lambda - sum(1-N)/(1-lambda)
}



#--SECOND DERIVATIVE FUNCTIONS for I(theta-hat)-------

#RHO 2
d2L_drho2_1d <- function(y_idx, MAPs,MCMCs,r){
  rho <- MAPs$rho[y_idx]
  sigma_v <- MAPs$sigma_v[y_idx]
  epsilon_y <- get_y_resid1d(y_idx, MAPs, MCMCs,r)
  epsilon_v <- get_v_resid1d(y_idx, MAPs, MCMCs,r)
  v <- MCMCs$v[r,]
  #derivative return
  T*(1+rho^2)/((1-rho^2)^2) - 
    sum((1 + 3*rho^2)*epsilon_y^2/(v[-T]*(1-rho^2)^3) - 
          2*rho*(3 + rho^2)*epsilon_y*epsilon_v/(sigma_v*v[-T]*(1-rho^2)^3)  +
          (1 + 3*rho^2)*epsilon_v^2/(sigma_v^2 * v[-T]*(1-rho^2)^3) )
}

#MU 2
d2L_dmu2_1d <- function(y_idx, MAPs,MCMCs,r){
  rho <- MAPs$rho[y_idx]
  sigma_v <- MAPs$sigma_v[y_idx]
  epsilon_y <- get_y_resid1d(y_idx, MAPs, MCMCs,r)
  epsilon_v <- get_v_resid1d(y_idx,MAPs, MCMCs,r)
  v <- MCMCs$v[r,]
  #derivative return
  -sum(1/(v[-T]*(1-rho^2)))
}

#THETA 2
d2L_dtheta2_1d <- function(y_idx, MAPs, MCMCs, r){
  rho <- MAPs$rho[y_idx]
  sigma_v <- MAPs$sigma_v[y_idx]
  kappa <- 1-MAPs$phi[y_idx]
  epsilon_y <- get_y_resid1d(y_idx, MAPs, MCMCs,r)
  epsilon_v <- get_v_resid1d(y_idx, MAPs, MCMCs,r)
  #derivative return  
  -sum(kappa*kappa/(v[-T]*(1-rho^2)*sigma_v^2) )
  
}

#KAPPA 2
d2L_dkappa2_1d <- function(y_idx, MAPs, MCMCs, r){
  rho <- MAPs$rho[y_idx]
  sigma_v <- MAPs$sigma_v[y_idx]
  kappa <- 1-MAPs$phi[y_idx]
  theta <- MAPs$theta[y_idx]
  epsilon_y <- get_y_resid1d(y_idx, MAPs, MCMCs,r)
  epsilon_v <- get_v_resid1d(y_idx,MAPs, MCMCs,r)
  #derivative return  
  -sum((theta - v[-T])*(theta - v[-T])/(sigma_v^2*v[-T]*(1-rho^2)))
  
}

#SIGMA_V 2
d2L_dsigma_v2_1d <- function(y_idx, MAPs, MCMCs, r){
  rho <- MAPs$rho[y_idx]
  sigma_v <- MAPs$sigma_v[y_idx]
  kappa <- 1-MAPs$phi[y_idx]
  theta <- MAPs$theta[y_idx]
  epsilon_y <- get_y_resid1d(y_idx, MAPs, MCMCs,r)
  epsilon_v <- get_v_resid1d(y_idx,MAPs, MCMCs,r)
  #derivative return  
  T/sigma_v^2 + 
    sum(2*rho*epsilon_y*epsilon_v/((sigma_v^3)*v[-T]*(1-rho*rho)) - 
          6*epsilon_v*epsilon_v/((sigma_v^4)*v[-T]*(1-rho*rho)))
  
}

#MU_Y (xi_yw) 2
d2L_dmuy2_1d <- function(y_idx, MAPs, MCMCs, r){
  B <- MCMCs$xi_ys[r,] #NEED TO SAVE XI_YS IN run_mcmc.R
  mu_y <- MAPs$xi_yw[y_idx]
  v_y <- MAPs$xi_yeta[y_idx]
  xi_y <- MCMCs$xi_y #NEED TO SAVE XI_Y IN run_mcmc.R
  
  sum(-(B)/(v_y*v_y))
  
}

#MU_Y (xi_yw) 2
d2L_dvy2_1d <- function(y_idx, MAPs, MCMCs, r){
  B <- MCMCs$xi_ys[r,] #NEED TO SAVE XI_YS IN run_mcmc.R
  mu_y <- MAPs$xi_yw[y_idx]
  v_y <- MAPs$xi_yeta[y_idx]
  xi_y <- MCMCs$xi_y #NEED TO SAVE XI_Y IN run_mcmc.R
  
  T/(v_y^2) -sum(3*(xi_y - mu_y*B)^2/((v_y^4)*B))
  
}

#LAMBDA
d2L_dlambda2_1d <- function(y_idx, MAPs, MCMCs, r){
  lambda <- MAPs$lambda[y_idx]
  N <- MCMCs$delta[r,]
  -sum(N)/(lambda^2) - sum(1-N)/(1-lambda)^2
}

#--d mu d "other"----------#
#D MU D THETA
d2L_dmudtheta_1d <- function(y_idx, MAPs,MCMCs,r){
  rho <- MAPs$rho[y_idx]
  kappa <- 1-MAPs$phi[y_idx]
  sigma_v <- MAPs$sigma_v[y_idx]
  epsilon_y <- get_y_resid1d(y_idx, MAPs, MCMCs,r)
  epsilon_v <- get_v_resid1d(y_idx,MAPs, MCMCs,r)
  v <- MCMCs$v[r,]
  #derivative return
  sum(rho*kappa/(v[-T]*(1-rho^2)*sigma_v) )
}

#D MU D KAPPA
d2L_dmudkappa_1d <- function(y_idx, MAPs,MCMCs,r){
  rho <- MAPs$rho[y_idx]
  theta <- MAPs$theta[y_idx]
  sigma_v <- MAPs$sigma_v[y_idx]
  epsilon_y <- get_y_resid1d(y_idx, MAPs, MCMCs,r)
  epsilon_v <- get_v_resid1d(y_idx,MAPs, MCMCs,r)
  v <- MCMCs$v[r,]
  #derivative return
  -sum((theta - v[-T])*rho/(v[-T]*(1-rho^2)*sigma_v) )
}

#D MU D SIGMA_V
d2L_dmudsigma_v_1d <- function(y_idx, MAPs,MCMCs,r){
  rho <- MAPs$rho[y_idx]
  sigma_v <- MAPs$sigma_v[y_idx]
  epsilon_y <- get_y_resid1d(y_idx, MAPs, MCMCs,r)
  epsilon_v <- get_v_resid1d(y_idx,MAPs, MCMCs,r)
  v <- MCMCs$v[r,]
  #derivative return
  sum(rho*epsilon_v/(v[-T]*(1-rho^2)*sigma_v*sigma_v) )
}

#D MU D RHO
d2L_dmudrho_1d <- function(y_idx, MAPs,MCMCs,r){
  rho <- MAPs$rho[y_idx]
  sigma_v <- MAPs$sigma_v[y_idx]
  epsilon_y <- get_y_resid1d(y_idx, MAPs, MCMCs,r)
  epsilon_v <- get_v_resid1d(y_idx,MAPs, MCMCs,r)
  v <- MCMCs$v[r,]
  #derivative return
  sum(2*rho*epsilon_y/(v[-T]*(1-rho^2)^2) -
         (1+rho)*epsilon_v/(v[-T]*((1-rho^2)^2)*sigma_v) )
}


#--d theta d "other"----------#

#D THETA D KAPPA
d2L_dthetadkappa_1d <- function(y_idx, MAPs,MCMCs,r){
rho <- MAPs$rho[y_idx]
sigma_v <- MAPs$sigma_v[y_idx]
kappa <- 1-MAPs$phi[y_idx]
epsilon_y <- get_y_resid1d(y_idx, MAPs, MCMCs,r)
epsilon_v <- get_v_resid1d(y_idx,MAPs, MCMCs,r)
#derivative return  
sum(epsilon_v/(v[-T]*(1-rho^2)*sigma_v^2) - 
       rho*epsilon_y/(v[-T]*(1-rho^2)*sigma_v))
}

#D THETA D SIGMA_V
d2L_dthetadsigma_v_1d <- function(y_idx, MAPs,MCMCs,r){
  rho <- MAPs$rho[y_idx]
  sigma_v <- MAPs$sigma_v[y_idx]
  kappa <- 1-MAPs$phi[y_idx]
  epsilon_y <- get_y_resid1d(y_idx, MAPs, MCMCs,r)
  epsilon_v <- get_v_resid1d(y_idx,MAPs, MCMCs,r)
  #derivative return  
  sum(-2*kappa*epsilon_v/(v[-T]*(1-rho^2)*sigma_v^3) + 
         rho*kappa*epsilon_y/(v[-T]*(1-rho^2)*sigma_v^2))
}

#D THETA D RHO
d2L_dthetadrho_1d <- function(y_idx, MAPs,MCMCs,r){
  rho <- MAPs$rho[y_idx]
  sigma_v <- MAPs$sigma_v[y_idx]
  kappa <- 1-MAPs$phi[y_idx]
  epsilon_y <- get_y_resid1d(y_idx, MAPs, MCMCs,r)
  epsilon_v <- get_v_resid1d(y_idx,MAPs, MCMCs,r)
  #derivative return  
  sum(2*kappa*epsilon_v*rho/(v[-T]*((1-rho^2)^2)*sigma_v^2) - 
         (1+rho)*kappa*epsilon_y/(v[-T]*((1-rho^2)^2)*sigma_v))
}


#--d kappa d "other"----------#

d2L_dkappadsigma_v_1d <- function(y_idx, MAPs, MCMCs, r){
  rho <- MAPs$rho[y_idx]
  sigma_v <- MAPs$sigma_v[y_idx]
  kappa <- 1-MAPs$phi[y_idx]
  theta <- MAPs$theta[y_idx]
  epsilon_y <- get_y_resid1d(y_idx, MAPs, MCMCs,r)
  epsilon_v <- get_v_resid1d(y_idx,MAPs, MCMCs,r)
  #derivative return  
  sum(-2*(theta - v[-T])*epsilon_v/((sigma_v^3)*v[-T]*(1-rho^2)) +
         rho*(theta - v[-T])*epsilon_y/((sigma_v^2)*v[-T]*(1-rho^2)))
  
}

d2L_dkappadrho_1d <- function(y_idx, MAPs, MCMCs, r){
  rho <- MAPs$rho[y_idx]
  sigma_v <- MAPs$sigma_v[y_idx]
  kappa <- 1-MAPs$phi[y_idx]
  theta <- MAPs$theta[y_idx]
  epsilon_y <- get_y_resid1d(y_idx, MAPs, MCMCs,r)
  epsilon_v <- get_v_resid1d(y_idx,MAPs, MCMCs,r)
  #derivative return  
  sum(2*rho*(theta - v[-T])*epsilon_v/(sigma_v^2*v[-T]*(1-rho^2)^2) -
         (1+rho)*(theta - v[-T]*epsilon_y)/sigma_v*v[-T]*(1-rho^2)^2)
  
}


#--d SIGMA_V d "other"----------#
#D SIGMA_V D RHO
d2L_dsigma_vdrho_1d <- function(y_idx, MAPs, MCMCs, r){
  rho <- MAPs$rho[y_idx]
  sigma_v <- MAPs$sigma_v[y_idx]
  kappa <- 1-MAPs$phi[y_idx]
  theta <- MAPs$theta[y_idx]
  epsilon_y <- get_y_resid1d(y_idx, MAPs, MCMCs,r)
  epsilon_v <- get_v_resid1d(y_idx, MAPs, MCMCs,r)
  #derivative return  
    sum(-(1 +rho^2)*epsilon_y*epsilon_v/(sigma_v*sigma_v*v[-T]*(1-rho*rho)^2) + 
          2*rho*epsilon_v*epsilon_v/(sigma_v*sigma_v*sigma_v*v[-T]*(1-rho*rho)^2))
  
}

#--d mu_y d v_y"----------#
d2L_dmuydvy_1d <- function(y_idx, MAPs, MCMCs, r){
  B <- MCMCs$xi_ys[r,] #NEED TO SAVE XI_YS IN run_mcmc.R
  mu_y <- MAPs$xi_yw[y_idx]
  v_y <- MAPs$xi_yeta[y_idx]
  xi_y <- MCMCs$xi_y #NEED TO SAVE XI_Y IN run_mcmc.R
  
  sum(-2*(xi_y - mu_y*B)/(v_y^3))
  
}



# Need E(1st ds | y, theta-hat), E(1st ds SQUARED | y, theta-hat), E(2nd ds | y, theta-hat)-----

get_S <- function(r){
  x <- c(
    dL_drho_1d(y_idx = 1, MAPs = THETA_map1d, MCMCs=keepsBTC_MAP, r=r),
    dL_drho_1d(y_idx = 2, MAPs = THETA_map1d, MCMCs=keepsSP_MAP , r=r),
    
    dL_dmu_1d(y_idx = 1, MAPs = THETA_map1d, MCMCs=keepsBTC_MAP, r=r),
    dL_dmu_1d(y_idx = 2, MAPs = THETA_map1d, MCMCs=keepsSP_MAP , r=r),
    
    dL_dtheta_1d(y_idx = 1, MAPs = THETA_map1d, MCMCs=keepsBTC_MAP, r=r),
    dL_dtheta_1d(y_idx = 2, MAPs = THETA_map1d, MCMCs=keepsSP_MAP , r=r),
    
    dL_dkappa_1d(y_idx = 1, MAPs = THETA_map1d, MCMCs=keepsBTC_MAP, r=r),
    dL_dkappa_1d(y_idx = 2, MAPs = THETA_map1d, MCMCs=keepsSP_MAP , r=r),
    
    dL_dsigma_v_1d(y_idx = 1, MAPs = THETA_map1d, MCMCs=keepsBTC_MAP, r=r),
    dL_dsigma_v_1d(y_idx = 2, MAPs = THETA_map1d, MCMCs=keepsSP_MAP , r=r),
    
    dL_dmuy_1d(y_idx = 1, MAPs = THETA_map1d, MCMCs=keepsBTC_MAP, r=r),
    dL_dmuy_1d(y_idx = 2, MAPs = THETA_map1d, MCMCs=keepsSP_MAP , r=r),
    
    dL_dvy_1d(y_idx = 1, MAPs = THETA_map1d, MCMCs=keepsBTC_MAP, r=r),
    dL_dvy_1d(y_idx = 2, MAPs = THETA_map1d, MCMCs=keepsSP_MAP , r=r),
    
    dL_dlambda_1d(y_idx = 1, MAPs = THETA_map1d, MCMCs=keepsBTC_MAP, r=r),
    dL_dlambda_1d(y_idx = 2, MAPs = THETA_map1d, MCMCs=keepsSP_MAP , r=r)
  ) %>% as.vector()
  return(x)
}

get_L2 <- function(r){
  drho1d_ <- c(
    d2L_drho2_1d(y_idx = 1, MAPs = THETA_map1d, MCMCs=keepsBTC_MAP, r=r),
    0,
    d2L_dmudrho_1d(y_idx = 1, MAPs = THETA_map1d, MCMCs=keepsBTC_MAP, r=r),
    0,
    d2L_dthetadrho_1d(y_idx = 1, MAPs = THETA_map1d, MCMCs=keepsBTC_MAP, r=r),
    0,
    d2L_dkappadrho_1d(y_idx = 1, MAPs = THETA_map1d, MCMCs=keepsBTC_MAP, r=r),
    0,
    d2L_dsigma_vdrho_1d(y_idx = 1, MAPs = THETA_map1d, MCMCs=keepsBTC_MAP, r=r),
    0,
    0,#d2L_drhodmuy_1d
    0,
    0,#d2L_drhodvy_1d
    0,
    0,#d2L_dlambda_1d,
    0
  ) %>% as.vector()
  
  drho2d_ <- c(
    0,
    d2L_drho2_1d(y_idx = 2, MAPs = THETA_map1d, MCMCs=keepsSP_MAP, r=r),
    0,
    d2L_dmudrho_1d(y_idx = 2, MAPs = THETA_map1d, MCMCs=keepsSP_MAP, r=r),
    0,
    d2L_dthetadrho_1d(y_idx = 2, MAPs = THETA_map1d, MCMCs=keepsSP_MAP, r=r),
    0,
    d2L_dkappadrho_1d(y_idx = 2, MAPs = THETA_map1d, MCMCs=keepsSP_MAP, r=r),
    0,
    d2L_dsigma_vdrho_1d(y_idx = 2, MAPs = THETA_map1d, MCMCs=keepsSP_MAP, r=r),
    0,#d2L_drhodmuy_1d
    0,
    0,#d2L_drhodvy_1d
    0,
    0,#d2L_dlambda_1d,
    0
  ) %>% as.vector()
  
  dmu1d_ <-c(
    d2L_dmudrho_1d(y_idx = 1, MAPs = THETA_map1d, MCMCs=keepsBTC_MAP, r=r),
    0,
    d2L_dmu2_1d(y_idx = 1, MAPs = THETA_map1d, MCMCs=keepsBTC_MAP, r=r),
    0,
    d2L_dmudtheta_1d(y_idx = 1, MAPs = THETA_map1d, MCMCs=keepsBTC_MAP, r=r),
    0,
    d2L_dmudkappa_1d(y_idx = 1, MAPs = THETA_map1d, MCMCs=keepsBTC_MAP, r=r),
    0,
    d2L_dmudsigma_v_1d(y_idx = 1, MAPs = THETA_map1d, MCMCs=keepsBTC_MAP, r=r),
    0,
    0,#d2L_dmudmuy_1d
    0,
    0,#d2L_dmudvy_1d
    0,
    0,#d2L_dmudlambda_1d,
    0
  ) %>% as.vector()
  
  dmu2d_ <-c(
    0,
    d2L_dmudrho_1d(y_idx = 2, MAPs = THETA_map1d, MCMCs=keepsSP_MAP, r=r),
    0,
    d2L_dmu2_1d(y_idx = 2, MAPs = THETA_map1d, MCMCs=keepsSP_MAP, r=r),
    0,
    d2L_dmudtheta_1d(y_idx = 2, MAPs = THETA_map1d, MCMCs=keepsSP_MAP, r=r),
    0,
    d2L_dmudkappa_1d(y_idx = 2, MAPs = THETA_map1d, MCMCs=keepsSP_MAP, r=r),
    0,
    d2L_dmudsigma_v_1d(y_idx = 2, MAPs = THETA_map1d, MCMCs=keepsSP_MAP, r=r),
    0,#d2L_dmudmuy_1d
    0,
    0,#d2L_dmudvy_1d
    0,
    0,#d2L_dmudlambda_1d,
    0
  ) %>% as.vector()
  
  dtheta1d_ <-c(
    d2L_dthetadrho_1d(y_idx = 1, MAPs = THETA_map1d, MCMCs=keepsBTC_MAP, r=r),
    0,
    d2L_dmudtheta_1d(y_idx = 1, MAPs = THETA_map1d, MCMCs=keepsBTC_MAP, r=r),
    0,
    d2L_dtheta2_1d(y_idx = 1, MAPs = THETA_map1d, MCMCs=keepsBTC_MAP, r=r),
    0,
    d2L_dthetadkappa_1d(y_idx = 1, MAPs = THETA_map1d, MCMCs=keepsBTC_MAP, r=r),
    0,
    d2L_dthetadsigma_v_1d(y_idx = 1, MAPs = THETA_map1d, MCMCs=keepsBTC_MAP, r=r),
    0,
    0,#d2L_dmudmuy_1d
    0,
    0,#d2L_dmudvy_1d
    0,
    0,#d2L_dmudlambda_1d,
    0
  ) %>% as.vector()
  
  
  
  dtheta2d_ <-c(0,
                d2L_dthetadrho_1d(y_idx = 2, MAPs = THETA_map1d, MCMCs=keepsSP_MAP, r=r),
                0,
                d2L_dmudtheta_1d(y_idx = 2, MAPs = THETA_map1d, MCMCs=keepsSP_MAP, r=r),
                0,
                d2L_dtheta2_1d(y_idx = 2, MAPs = THETA_map1d, MCMCs=keepsSP_MAP, r=r),
                0,
                d2L_dthetadkappa_1d(y_idx = 2, MAPs = THETA_map1d, MCMCs=keepsSP_MAP, r=r),
                0,
                d2L_dthetadsigma_v_1d(y_idx = 2, MAPs = THETA_map1d, MCMCs=keepsSP_MAP, r=r),
                0,#d2L_dmudmuy_1d
                0,
                0,#d2L_dmudvy_1d
                0,
                0,#d2L_dmudlambda_1d,
                0
  ) %>% as.vector()
  
  
  dkappa1d_ <-c(
    d2L_dkappadrho_1d(y_idx = 1, MAPs = THETA_map1d, MCMCs=keepsBTC_MAP, r=r),
    0,
    d2L_dmudkappa_1d(y_idx = 1, MAPs = THETA_map1d, MCMCs=keepsBTC_MAP, r=r),
    0,
    d2L_dthetadkappa_1d(y_idx = 1, MAPs = THETA_map1d, MCMCs=keepsBTC_MAP, r=r),
    0,
    d2L_dkappa2_1d(y_idx = 1, MAPs = THETA_map1d, MCMCs=keepsBTC_MAP, r=r),
    0,
    d2L_dkappadsigma_v_1d(y_idx = 1, MAPs = THETA_map1d, MCMCs=keepsBTC_MAP, r=r),
    0,
    0,#d2L_dmudmuy_1d
    0,
    0,#d2L_dmudvy_1d
    0,
    0,#d2L_dmudlambda_1d,
    0
  ) %>% as.vector()
  
  
  dkappa2d_ <-c(0,
                d2L_dkappadrho_1d(y_idx = 2, MAPs = THETA_map1d, MCMCs=keepsSP_MAP, r=r),
                0,
                d2L_dmudkappa_1d(y_idx = 2, MAPs = THETA_map1d, MCMCs=keepsSP_MAP, r=r),
                0,
                d2L_dthetadkappa_1d(y_idx = 2, MAPs = THETA_map1d, MCMCs=keepsSP_MAP, r=r),
                0,
                d2L_dkappa2_1d(y_idx = 2, MAPs = THETA_map1d, MCMCs=keepsSP_MAP, r=r),
                0,
                d2L_dkappadsigma_v_1d(y_idx = 2, MAPs = THETA_map1d, MCMCs=keepsSP_MAP, r=r),
                0,#d2L_dmudmuy_1d
                0,
                0,#d2L_dmudvy_1d
                0,
                0,#d2L_dmudlambda_1d,
                0
  ) %>% as.vector()
  
  
  dsigma_v1d_ <-c(
    d2L_dsigma_vdrho_1d(y_idx = 1, MAPs = THETA_map1d, MCMCs=keepsBTC_MAP, r=r),
    0,
    d2L_dmudsigma_v_1d(y_idx = 1, MAPs = THETA_map1d, MCMCs=keepsBTC_MAP, r=r),
    0,
    d2L_dthetadsigma_v_1d(y_idx = 1, MAPs = THETA_map1d, MCMCs=keepsBTC_MAP, r=r),
    0,
    d2L_dkappadsigma_v_1d(y_idx = 1, MAPs = THETA_map1d, MCMCs=keepsBTC_MAP, r=r),
    0,
    d2L_dsigma_v2_1d(y_idx = 1, MAPs = THETA_map1d, MCMCs=keepsBTC_MAP, r=r),
    0,
    0,#d2L_dmudmuy_1d
    0,
    0,#d2L_dmudvy_1d
    0,
    0,#d2L_dmudlambda_1d,
    0
  ) %>% as.vector()
  
  
  dsigma_v2d_ <-c(0,
                  d2L_dsigma_vdrho_1d(y_idx = 2, MAPs = THETA_map1d, MCMCs=keepsSP_MAP, r=r),
                  0,
                  d2L_dmudsigma_v_1d(y_idx = 2, MAPs = THETA_map1d, MCMCs=keepsSP_MAP, r=r),
                  0,
                  d2L_dthetadsigma_v_1d(y_idx = 2, MAPs = THETA_map1d, MCMCs=keepsSP_MAP, r=r),
                  0,
                  d2L_dkappadsigma_v_1d(y_idx = 2, MAPs = THETA_map1d, MCMCs=keepsSP_MAP, r=r),
                  0,
                  d2L_dsigma_v2_1d(y_idx = 2, MAPs = THETA_map1d, MCMCs=keepsSP_MAP, r=r),
                  0,#d2L_dmudmuy_1d
                  0,
                  0,#d2L_dmudvy_1d
                  0,
                  0,#d2L_dmudlambda_1d,
                  0
  ) %>% as.vector()
  
  dmuy1d_ <-c(
    rep(0,10),
    d2L_dmuy2_1d(y_idx = 1, MAPs = THETA_map1d, MCMCs=keepsBTC_MAP, r=r),
    0,
    d2L_dmuydvy_1d(y_idx = 1, MAPs = THETA_map1d, MCMCs=keepsBTC_MAP, r=r),#d2L_dmuydvy_1d
    0,
    0,#d2L_dmuydlambda_1d,
    0
  ) %>% as.vector()
  
  dmuy2d_ <-c(
    rep(0,10),
    0,
    d2L_dmuy2_1d(y_idx = 2, MAPs = THETA_map1d, MCMCs=keepsSP_MAP, r=r),
    0,
    d2L_dmuydvy_1d(y_idx = 2, MAPs = THETA_map1d, MCMCs=keepsSP_MAP, r=r),#d2L_dmuydvy_1d
    0,#d2L_dmuydlambda_1d,
    0
  ) %>% as.vector()
  
  
  dvy1d_ <-c(
    rep(0,10),
    d2L_dmuydvy_1d(y_idx = 1, MAPs = THETA_map1d, MCMCs=keepsBTC_MAP, r=r),
    0,
    d2L_dvy2_1d(y_idx = 1, MAPs = THETA_map1d, MCMCs=keepsBTC_MAP, r=r),#d2L_dmuydvy_1d
    0,
    0,#d2L_dmuydlambda_1d,
    0
  ) %>% as.vector()
  
  dvy2d_ <-c(
    rep(0,10),
    0,
    d2L_dmuydvy_1d(y_idx = 2, MAPs = THETA_map1d, MCMCs=keepsSP_MAP, r=r),
    0,
    d2L_dvy2_1d(y_idx = 2, MAPs = THETA_map1d, MCMCs=keepsSP_MAP, r=r),#d2L_dmuydvy_1d
    0,#d2L_dmuydlambda_1d,
    0
  ) %>% as.vector()
  
  dlambda1d_ <-c(
    rep(0,14),
    d2L_dlambda2_1d(y_idx = 1, MAPs = THETA_map1d, MCMCs=keepsBTC_MAP, r=r),#d2L_dmuydlambda_1d,
    0
  ) %>% as.vector()
  dlambda2d_ <-c(
    rep(0,14),
    0,
    d2L_dlambda2_1d(y_idx = 2, MAPs = THETA_map1d, MCMCs=keepsSP_MAP, r=r)#d2L_dmuydlambda_1d,
    
  ) %>% as.vector()
  
  L2 <- rbind(drho1d_,drho2d_, 
              dmu1d_,dmu2d_,
              dtheta1d_,dtheta2d_,
              dkappa1d_,dkappa2d_,
              dsigma_v1d_,dsigma_v2d_,
              dmuy1d_,dmuy2d_,
              dvy1d_,dvy2d_,
              dlambda1d_,dlambda2d_)
  return(L2)
}


#----E(S), E(S^2), E(2nd derivative matrix)------

ES <- 0 #expected value of the S matrix (will square AFTER)
ES2 <- 0 #expected value of S^2
EL2 <- 0 #expected value of the matrix of second derivatives
for (r in 1:R){
  print(r)
  x <- getS(r)
  ES <- ES + x/R #page 16 RDIC paper
  ES2 <- ES2 + x%*%t(x)/R #page 16 RDIC paper
  EL2 <- L2 + get_L2(r)/R
}
ES %*% t(ES) #E(S|y)^2
ES2 #E(S^2)

#------FINALLY DO IT (calculate RDIC)-------
I <- -EL2 -ES2 + ES %*% t(ES)

#----DIC 7 calculations-----
#calculating DIC 7 from page 8 of Li, Zeng, Yu
lnpy_mid_zhatthhetahat <- sapply(1:2,function(i){
  dnorm(y[],
        x+mu+J+rho/sigma_v*c(v[-1]-theta-phi*(v[-(T+1)]-theta),0),
        sqrt(v*c(rep(1-rho^2,T-1),1)),log=TRUE) %>% sum
}) %>% sum
lnpy_mid_thetahat <- sapply(1:2,function(i){
  pg_lnp_y_mid_theta(y,x,mu,theta,phi,sigma_v,rho,xi_yw,xi_yeta,lambda,N=10000)
}) %>% sum

DIC7_1d = -4*Elnpy_mid_ztheta + 2*lnpy_mid_zhatthetahat
DIC7_1d
#16466.6

#----DIC 1 calculations-----
#calculating DIC 1 from page 8 of Li, Zeng, Yu
DIC1_1d <- -4*Elnpy_mid_theta + 2*lnpy_mid_thetahat #also anticlimactic...

# deviance evaluated at theta hat + 2* trace(I*V)
RDIC_1d = -2*lnpy_mid_thetahat + 2*((I%*%VCOV_1d) %>% diag() %>% sum())# how anticlimactic 


results_1d <- data.frame(estimate = c(DIC1_1d, DIC7_1d, RDIC_1d),
                         IC = c("DIC1", "DIC7", "RDIC"))
saveRDS(results1d)


