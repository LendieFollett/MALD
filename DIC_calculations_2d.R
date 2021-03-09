
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

# DATA AND MCMC SAMPLES ------------

getSymbols("BTC-USD",from = "2014-09-15",to = "2020-09-30")
BTC <- as.data.frame(`BTC-USD`)
BTC$Date <- seq(as.Date("2014-09-17"),as.Date("2020-09-30"),by="days")
BTC$`BTC-USD.Close`[BTC$Date=="2020-04-17"] <- 7096.18
getSymbols("^GSPC",from = "2014-09-15",to = "2020-09-30")
SP500 <- as.data.frame(`GSPC`)
SP500$Date <- as.Date(rownames(SP500))
Rsequence <- seq(1,nrow(keepsBTC$v), by = 100)
S <- merge(BTC,SP500)
T <- nrow(S) - 1
S <- merge(BTC,SP500)
y <- as.matrix(100*(log(S[-1,c("BTC-USD.Close","GSPC.Close")]) - log(S[-nrow(S),c("BTC-USD.Close","GSPC.Close")])))
#download full MCMC results (obtained from run_mcmc.R with fix = FALSE)
keepsBTCSP <- readRDS("/Users/000766412/Box Sync/ALD_Codes/keepsBTCSP.rds")

#DIC7 = -4E(ln(p(y|z,theta))) + 2ln(p(y|z-hat, theta-hat))
partial_likelihood1 <-  function(k, y){ #k for keeps, y for correct vector
  R <- dim(k$v[complete.cases(k$v[,1,1]),,])[1]
  total = 0
  for (r in 1:R){
    print(r)
    total_sub = 0
    for (t in 1:nrow(y)){
    total_sub = total_sub + (dmvnorm(y[t,], 
                           k$mu[r,] + k$J[r,t,],
                           diag(sqrt(k$v[r,t+1,]))%*%(k$rho[r]*(1-diag(2)) + diag(2))%*%diag(sqrt(k$v[r,t+1,])), 
                           log = TRUE)) 
    }
    total = total + total_sub/R #average over R draws
    
  }
  return(total)
}

#use plug-in posterior means of z, theta
partial_likelihood2 <-  function(k, y){ #k for keeps, y for correct vector

  mu_mean =  apply(k$mu, 2, mean, na.rm=TRUE)
  J_mean = apply(k$J,2:3, mean,na.rm=TRUE)
  v_mean = apply(k$v,2:3, mean,na.rm=TRUE)
  rho_mean = mean(k$rho, na.rm=TRUE)
  total_sub = 0
  for (t in 1:nrow(y)){
    total_sub = total_sub + (dmvnorm(y[t,], 
                                     mu_mean+ J_mean[t,],
                                     diag(sqrt(v_mean[r,t+1,]))%*%(rho_mean[2]*(1-diag(2)) + diag(2))%*%diag(sqrt(v_mean[r,t+1,])), 
                                     log = TRUE)) 
  }
  return(total_sub)
}


#E(ln(p(y|z,theta)))
Elnpy_mid_ztheta <- partial_likelihood1(keepsBTCSP, y) 

lnpy_mid_zhatthetahat <- partial_likelihood2(keepsBTCSP, y) 


DIC7_1d = -4*Elnpy_mid_ztheta + 2*lnpy_mid_zhatthetahat
DIC7_1d



