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
partial_likelihood <-  function(k, y){ #k for keeps, y for correct vector
  R <- dim(k$v[complete.cases(k$v),])[1]
  total = 0
  for (r in 1:R){
  total = total + (dnorm(y, 
        k$mu[r] + k$J[r,],
        sqrt(k$v[r,-1]), log = TRUE) %>%sum)/R #average over R draws
  }
  return(total)
}

pl <- partial_likelihood(keepsBTC, y[,1]) + partial_likelihood(keepsSP, y[,2])


lnpy_mid_zhatthetahat <- sapply(1:2,function(i){
  dnorm(y[],
        x+mu+J+rho/sigma_v*c(v[-1]-theta-phi*(v[-(T+1)]-theta),0),
        sqrt(v*c(rep(1-rho^2,T-1),1)),log=TRUE) %>% sum
}) %>% sum



DIC7_1d = -4*Elnpy_mid_ztheta + 2*lnpy_mid_zhatthetahat
DIC7_1d



