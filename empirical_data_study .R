rm(list = ls())  ## DON'T FORGET TO SET WD!!
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


thin <- 5 #thinning param
B <- 10000 #how many burn in draws to throw away
R <- 100000 #how many draws to keep after burn in
n_chns <- 1 #how many chains to run
#load data
getSymbols("BTC-USD",from = "2014-09-15",to = "2020-09-30")
BTC <- as.data.frame(`BTC-USD`)
BTC$Date <- seq(as.Date("2014-09-17"),as.Date("2020-09-30"),by="days")
BTC$`BTC-USD.Close`[BTC$Date=="2020-04-17"] <- 7096.18
getSymbols("^GSPC",from = "2014-09-15",to = "2020-09-30")
SP500 <- as.data.frame(`GSPC`)
SP500$Date <- as.Date(rownames(SP500))
S <- merge(BTC,SP500)
T <- nrow(S) - 1


#################################################### 
# SVMALD MODEL ---------- LRF RUNS
#################################################### 
use_starting_values <- TRUE
sourceCpp("pgas_2d.cpp") #C++ updates
# #2-D MODEL MCMCb        cfv09
y <- as.matrix(100*(log(S[-1,c("BTC-USD.Close","GSPC.Close")]) - log(S[-nrow(S),c("BTC-USD.Close","GSPC.Close")])))
yprim <- array(0,dim=dim(y))
#source("starting_values_2d.R") #initialize values (performed within run_mcmc_2d.R)
exp_jumps <- norm_jumps <- ind <- FALSE
source("run_mcmc_2d.R") #R+B iterations of pgas.R and pgas.cpp updates
saveRDS(keeps,paste0("keeps_long/keepsBTCSP.rds"))

#################################################### 
# SVALD INDEPENDENCE (1d) MODEL ---------- LRF RUNS
#################################################### 
use_starting_values <- TRUE
sourceCpp("pgas_2d.cpp") #C++ updates
#initialize values, create space to save draws
# #2-D MODEL MCMC
y <- as.matrix(100*(log(S[-1,c("BTC-USD.Close","GSPC.Close")]) - log(S[-nrow(S),c("BTC-USD.Close","GSPC.Close")])))
yprim <- array(0,dim=dim(y))
#source("starting_values_2d.R") #initialize values
exp_jumps  <- FALSE 
norm_jumps <- FALSE 
ind <- TRUE #Bitcoin and S&P 500 have no relationship in return, volatility or jumps
source("run_mcmc_2d.R") #R+B iterations of pgas.R and pgas.cpp updates
saveRDS(keeps,paste0("keeps_long/keepsBTCSP_IND.rds"))

#################################################### 
# SVMVN MODEL ---------- MS RUNS
#################################################### 
sourceCpp("pgas_2d.cpp") #C++ updates
#initialize values, create space to save draws
# #2-D MODEL MCMC
y <- as.matrix(100*(log(S[-1,c("BTC-USD.Close","GSPC.Close")]) - log(S[-nrow(S),c("BTC-USD.Close","GSPC.Close")])))
yprim <- array(0,dim=dim(y))
source("starting_values_2d.R") #initialize values
exp_jumps <- FALSE
norm_jumps <- TRUE #NORMAL JUMPS SET TO TRUE SO B/s are set to 1
ind <- FALSE
source("run_mcmc_2d.R") #R+B iterations of pgas.R and pgas.cpp updates
saveRDS(keeps,paste0("keeps_long/keepsBTCSP_MVN.rds"))

#################################################### 
# SVLD MODEL ----------MS RUNS
#################################################### 
sourceCpp("pgas_2d.cpp") #C++ updates
#initialize values, create space to save draws
# #2-D MODEL MCMC
y <- as.matrix(100*(log(S[-1,c("BTC-USD.Close","GSPC.Close")]) - log(S[-nrow(S),c("BTC-USD.Close","GSPC.Close")])))
yprim <- array(0,dim=dim(y))
source("starting_values_2d.R") #initialize values
exp_jumps  <- TRUE #ASYMMETRY PARAMETERS W SET TO 0 (exponential, laplace distributed jumps)
norm_jumps <- FALSE 
ind <- FALSE
source("run_mcmc_2d.R") #R+B iterations of pgas.R and pgas.cpp updates
saveRDS(keeps,paste0("keeps_long/keepsBTCSP_LD.rds"))


#################################################### 
# CONVERGENCE CHECKS ----------
#################################################### 
library(LaplacesDemon)

total <- 20000 #number of mcmc iterations saved after burn-in, thinning
doESS <- function(x, total){
  R <- total
  if(!is.null(dim(x))){ #if it's a data frame
    return(apply(x[1:R,], 2, ESS))
  }else{
    return(ESS(x[1:R]))
  }
}


#SVMALD
lapply(keeps[c(4,6:17)], doESS, total = 20000) %>% str()
lapply(keepsBTCSP[c(4,6:17)], domean, total = 20000) %>% str()

plot(keeps$sigma_c[,1], type = "l")
plot(keeps$sigma_c[,2], type = "l")
plot(keeps$rhoc, type = "l")
plot(keeps$xi_cw[,1], type = "l")
plot(keeps$xi_cw[,2], type = "l")
plot(keeps$xi_y1eta, type = "l")
plot(keeps$xi_y2eta, type = "l")
plot(keeps$xi_y1w, type = "l")
plot(keeps$xi_y2w, type = "l")
#SVALD
lapply(keepsIND[c(4,6:17)], doESS, total = total) %>% str()
#SVMVN
lapply(keepsBTCSP_MVN[c(4,6:17)], doESS, total = total) %>% str()
#SVLD
lapply(keepsBTCSP_LD[c(4,6:17)], doESS, total = total) %>% str()



plot(keepsBTCSP$sigma_c[1:total, 1]);length(unique(keepsBTCSP$sigma_c[1:total, 1]))/total
plot(keepsBTCSP$sigma_c[1:total, 1])
plot(keepsBTCSP$rhoc[1:total])
plot(keepsBTCSP$xi_y2eta[1:total])

#use for starting values?

starting_values <- lapply(keeps, domean, total = 20000)
starting_values %>% str
starting_values$delta <- apply(starting_values$delta , 2, round)
starting_values$xi_y1 <- starting_values$xi_y1c <- starting_values$J[,1] + rnorm(length(starting_values$xi_y1), 0, .001)
starting_values$xi_y2 <- starting_values$xi_y2c <- starting_values$J[,2]+ rnorm(length(starting_values$xi_y1), 0, .001)
saveRDS(starting_values,"starting_values_MALD.rds")


