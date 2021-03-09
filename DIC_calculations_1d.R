
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
keepsBTC <- readRDS("/Users/000766412/Box Sync/ALD_Codes/keepsBTC.rds")
keepsSP <- readRDS("/Users/000766412/Box Sync/ALD_Codes/keepsSP.rds")


#DIC7 = -4E(ln(p(y|z,theta))) + 2ln(p(y|z-hat, theta-hat))
partial_likelihood1 <-  function(k, y){ #k for keeps, y for correct vector
  R <- dim(k$v[complete.cases(k$v),])[1]
  total = 0
  for (r in 1:R){
  total = total + (dnorm(y, 
        k$mu[r] + k$J[r,]+
          (k$rho[r]/k$sigma_v[r])*(k$v[r,-1]-k$theta[r]-k$phi[r]*(k$v[r,-(T+1)]-k$theta[r])),
        sqrt(k$v[r,-1]*(1-k$rho[r]^2)), log = TRUE) %>%sum)/R #average over R draws
  }
  return(total)
}

#use plug-in posterior means of z, theta
partial_likelihood2 <-  function(k, y){ #k for keeps, y for correct vector
  v_mean = apply(k$v, 2, mean)
  theta_mean <- mean(k$theta)
  phi_mean <- mean(k$phi)
  T <- nrow(y)
 dnorm(y, 
       mean(k$mu) + apply(k$J,2, mean) + 
         (mean(k$rho)/mean(k$sigma_v))*(v_mean[-1]-theta_mean-phi_mean*(v_mean[-(T+1)]-theta_mean)),
       sqrt(v_mean[-1]*(1-mean(k$rho)^2)), log = TRUE) %>%sum

}


#E(ln(p(y|z,theta)))
Elnpy_mid_ztheta <- partial_likelihood1(keepsBTC, y[,1]) + partial_likelihood1(keepsSP, y[,2])

lnpy_mid_zhatthetahat <- partial_likelihood2(keepsBTC, y[,1]) + partial_likelihood2(keepsSP, y[,2])



DIC7_1d = -4*Elnpy_mid_ztheta + 2*lnpy_mid_zhatthetahat
DIC7_1d
#17878.61


