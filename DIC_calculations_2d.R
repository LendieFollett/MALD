
rm(list = ls())
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

#DIC7 = -4E(ln(p(y|z,theta))) + 2ln(p(y|z-hat, theta-hat))
filepath <- "/Users/000766412/OneDrive - Drake University/Documents/Research/Asymmetric Laplace Jumps/keeps/"
# DATA AND MCMC SAMPLES ------------
#download full MCMC results (obtained from run_mcmc.R with fix = FALSE)
keepsIND <- readRDS(paste0(filepath,"keepsBTCSP_IND.rds")) #independence
keepsBTCSP <- readRDS(paste0(filepath,"keepsBTCSP.rds")) #MALD jump;s
keepsBTCSP_MVN <- readRDS(paste0(filepath,"keepsBTCSP_MVN.rds")) #multivariate normal jumps
keepsBTCSP_LD <- readRDS(paste0(filepath,"keepsBTCSP_LD.rds")) #laplacian jumps
getSymbols("BTC-USD",from = "2014-09-15",to = "2020-09-30")
BTC <- as.data.frame(`BTC-USD`)
BTC$Date <- seq(as.Date("2014-09-17"),as.Date("2020-09-30"),by="days")
BTC$`BTC-USD.Close`[BTC$Date=="2020-04-17"] <- 7096.18
getSymbols("^GSPC",from = "2014-09-15",to = "2020-09-30")
SP500 <- as.data.frame(`GSPC`)
SP500$Date <- as.Date(rownames(SP500))

T <- nrow(S) - 1
S <- merge(BTC,SP500)
y <- as.matrix(100*(log(S[-1,c("BTC-USD.Close","GSPC.Close")]) - log(S[-nrow(S),c("BTC-USD.Close","GSPC.Close")])))

#----2 D FUNCTIONS-----

#DIC7 = -4E(ln(p(y|z,theta))) + 2ln(p(y|z-hat, theta-hat))
partial_likelihood1 <-  function(k, y){ #k for keeps, y for correct vector
  R <- dim(k$v[complete.cases(k$v[,1,1]),,])[1]
  total = 0
  for (r in 1:R){
    print(r)
    total_sub = 0
    for (t in 1:nrow(y)){
      Sigma11 <- matrix(c(k$v[r,t,1],
                          k$rho[r,1]*sqrt(prod(k$v[r,t,])),
                          k$rho[r,1]*sqrt(prod(k$v[r,t,])),
                          k$v[r,t,2]),nrow=2)
      Sigma22 <- matrix(c(k$sigma_v[r,1]^2*k$v[r,t,1],
                          k$rho[r,2]*prod(k$sigma_v[r,])*sqrt(prod(k$v[r,t,])),
                          k$rho[r,2]*prod(k$sigma_v[r,])*sqrt(prod(k$v[r,t,])),
                          k$sigma_v[r,2]^2*k$v[r,t,2]),nrow=2)
      Sigma12 <- diag(c(k$rho[r,3:4]*k$sigma_v[r,]*k$v[r,t,]))
      eps <- k$v[r,t+1,] - k$theta[r,] - k$phi[r,] * (k$v[r,t,] - k$theta[r,])
    total_sub = total_sub + (dmvnorm(y[t,], 
                                     k$mu[r,] + k$J[r,t,] + Sigma12 %*% solve(Sigma22) %*% eps,
                                     Sigma11 - Sigma12 %*% solve(Sigma22) %*% Sigma12, 
                                     log = TRUE)) 
    }
    total = total + total_sub/R #average over R draws
    
  }
  return(total)
}

#use plug-in posterior means of z, theta
partial_likelihood2 <-  function(k, y){ #k for keeps, y for correct vector

  mu_mean =  apply(k$mu, 2, mean, na.rm=TRUE)
  theta_mean = apply(k$theta, 2, mean, na.rm=TRUE)
  phi_mean = apply(k$phi, 2, mean, na.rm=TRUE)
  sigma_v_mean = apply(k$sigma_v, 2, mean, na.rm=TRUE)
  rho_mean = apply(k$rho, 2, mean, na.rm=TRUE)
  J_mean = apply(k$J,2:3, mean,na.rm=TRUE)
  v_mean = apply(k$v,2:3, mean,na.rm=TRUE)
  total_sub = 0
  for (t in 1:nrow(y)){
    print(t)
    Sigma11 <- matrix(c(v_mean[t,1],
                        rho_mean[1]*sqrt(prod(v_mean[t,])),
                        rho_mean[1]*sqrt(prod(v_mean[t,])),
                        v_mean[t,2]),nrow=2)
    Sigma22 <- matrix(c(sigma_v_mean[1]^2*v_mean[t,1],
                        rho_mean[2]*prod(sigma_v_mean)*sqrt(prod(v_mean[t,])),
                        rho_mean[2]*prod(sigma_v_mean)*sqrt(prod(v_mean[t,])),
                        sigma_v_mean[2]^2*v_mean[t,2]),nrow=2)
    Sigma12 <- diag(c(rho_mean[3:4]*sigma_v_mean*v_mean[t,]))
    eps <- v_mean[t+1,] - theta_mean - phi_mean * (v_mean[t,] - theta_mean)
    total_sub = total_sub + (dmvnorm(y[t,], 
                                     mu_mean + J_mean[t,] + Sigma12 %*% solve(Sigma22) %*% eps,
                                     Sigma11 - Sigma12 %*% solve(Sigma22) %*% Sigma12, 
                                     log = TRUE)) 
  }
  return(total_sub)
}



#----1 D FUNCTIONS----
partial_likelihood1_1d <-  function(k, y){ #k for keeps, y for correct vector
  R <- dim(k$v[complete.cases(k$v),])[1]
  total = 0
  for (r in 1:R){
    total = total + (dnorm(y, 
                           k$mu[r] + k$J[r,]+
                             (k$rho[r]/k$sigma_v[r])*(k$v[r,-1]-k$theta[r]-k$phi[r]*(k$v[r,-(T+1)]-k$theta[r])),
                           sqrt(k$v[r,-(T+1)]*(1-k$rho[r]^2)), log = TRUE) %>%sum)/R #average over R draws
  }
  return(total)
}

#use plug-in posterior means of z, theta
partial_likelihood2_1d  <-  function(k, y){ #k for keeps, y for correct vector
  v_mean = apply(k$v, 2, mean)
  theta_mean <- mean(k$theta)
  phi_mean <- mean(k$phi)
  T <- length(y)
  dnorm(y, 
        mean(k$mu) + apply(k$J,2, mean) + 
          (mean(k$rho)/mean(k$sigma_v))*(v_mean[-1]-theta_mean-phi_mean*(v_mean[-(T+1)]-theta_mean)),
        sqrt(v_mean[-(T+1)]*(1-mean(k$rho)^2)), log = TRUE) %>%sum
  
}



#-----CALCULATIONS------------
#----SVALD
#E(ln(p(y|z,theta)))
Elnpy_mid_ztheta <- partial_likelihood1(keepsIND, y)
#[1]
lnpy_mid_zhatthetahat <- partial_likelihood2(keepsIND, y) 
#[1] 

DIC7_1d = -4*Elnpy_mid_ztheta + 2*lnpy_mid_zhatthetahat
DIC7_1d
#[1] 

#----SVMALD
#E(ln(p(y|z,theta)))
Elnpy_mid_ztheta <- partial_likelihood1(keepsBTCSP, y) 
#[1]
lnpy_mid_zhatthetahat <- partial_likelihood2(keepsBTCSP, y) 
#[1] 
DIC7_2d = -4*Elnpy_mid_ztheta + 2*lnpy_mid_zhatthetahat
DIC7_2d
#[1] 
#compare to Xxxxx from the 1d model... prefer 1d

#----SVMVN

Elnpy_mid_ztheta <- partial_likelihood1(keepsBTCSP_MVN, y) 
#[1] 
lnpy_mid_zhatthetahat <- partial_likelihood2(keepsBTCSP_MVN, y) 
#[1] 
DIC7_MVN = -4*Elnpy_mid_ztheta + 2*lnpy_mid_zhatthetahat
DIC7_MVN
#[1] 

#----SVLD

Elnpy_mid_ztheta <- partial_likelihood1(keepsBTCSP_LD, y) 
#[1] 
lnpy_mid_zhatthetahat <- partial_likelihood2(keepsBTCSP_LD, y) 
#[1] 
DIC7_LD= -4*Elnpy_mid_ztheta + 2*lnpy_mid_zhatthetahat
DIC7_LD
#[1]  
