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

k <- 5 #thinning param
B <- 50000 #how many burn in draws to throw away
R <- 50000 #how many draws to keep after burn in
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
# SVMALD MODEL ----------
#################################################### 

sourceCpp("pgas_2d.cpp") #C++ updates
# #2-D MODEL MCMC
y <- as.matrix(100*(log(S[-1,c("BTC-USD.Close","GSPC.Close")]) - log(S[-nrow(S),c("BTC-USD.Close","GSPC.Close")])))
x <- array(0,dim=dim(y))
source("starting_values_2d.R") #initialize values
exp_jumps <- norm_jumps <- FALSE
source("run_mcmc_2d.R") #R+B iterations of pgas.R and pgas.cpp updates
saveRDS(keeps,paste0("/Users/000766412/Box Sync/ALD_Codes/keepsBTCSP.rds"))

#################################################### 
# SVALD INDEPENDENCE (1d) MODEL ----------
#################################################### 
sourceCpp("pgas_mstuart.cpp") #C++ updates
#for BTC
y <- 100*(log(S[-1,c("BTC-USD.Close")]) - log(S[-nrow(S),c("BTC-USD.Close")]))
x <- rep(0,length(y))
source("starting_values.R") #initialize values
source("run_mcmc.R") 
saveRDS(keeps,paste0("../keeps/keepsBTC.rds"))
#for SP
y <- 100*(log(S[-1,c("GSPC.Close")]) - log(S[-nrow(S),c("GSPC.Close")]))
x <- rep(0,length(y))
source("starting_values.R") #initialize values
source("run_mcmc.R") 
saveRDS(keeps,paste0("../keeps/keepsSP.rds"))


#################################################### 
# SVMVN MODEL ----------
#################################################### 
sourceCpp("pgas_2d.cpp") #C++ updates
#initialize values, create space to save draws
# #2-D MODEL MCMC
y <- as.matrix(100*(log(S[-1,c("BTC-USD.Close","GSPC.Close")]) - log(S[-nrow(S),c("BTC-USD.Close","GSPC.Close")])))
x <- array(0,dim=dim(y))
source("starting_values_2d.R") #initialize values
exp_jumps  <- FALSE
norm_jumps <- TRUE #NORMAL JUMPS SET TO TRUE SO B/s are set to 1
source("run_mcmc_2d.R") #R+B iterations of pgas.R and pgas.cpp updates
saveRDS(keeps,paste0("/Users/000766412/Box Sync/ALD_Codes/keepsBTCSP_MVN.rds"))

#################################################### 
# SVLD MODEL ----------
#################################################### 
sourceCpp("pgas_2d.cpp") #C++ updates
#initialize values, create space to save draws
# #2-D MODEL MCMC
y <- as.matrix(100*(log(S[-1,c("BTC-USD.Close","GSPC.Close")]) - log(S[-nrow(S),c("BTC-USD.Close","GSPC.Close")])))
x <- array(0,dim=dim(y))
source("starting_values_2d.R") #initialize values
exp_jumps  <- TRUE #ASYMMETRY PARAMETERS W SET TO 0 (exponential, laplace distributed jumps)
norm_jumps <- FALSE 
source("run_mcmc_2d.R") #R+B iterations of pgas.R and pgas.cpp updates
saveRDS(keeps,paste0("/Users/000766412/Box Sync/ALD_Codes/keepsBTCSP_LD.rds"))


#################################################### 
# CONVERGENCE CHECKS ----------
#################################################### 
library(LaplacesDemon)

total <- 10000 #number of mcmc iterations saved after burn-in, thinning
doESS <- function(x, total){
  R <- total
  if(!is.null(dim(x))){ #if it's a data frame
    return(apply(x[1:R,], 2, ESS))
  }else{
    return(ESS(x[1:R]))
  }
}

#SVMALD
lapply(keepsBTCSP[c(2:3,5:16)], doESS, total = total) %>% str()
#SVMVN
lapply(keepsBTCSP_MVN[c(2:3,5:16)], doESS, total = total) %>% str()
#SVLD
lapply(keepsBTCSP_LD[c(2:3,5:16)], doESS, total = total) %>% str()



plot(keepsBTCSP$sigma_c[1:total, 1]);length(unique(keepsBTCSP$sigma_c[1:total, 1]))/total
plot(keepsBTCSP$sigma_c[1:total, 1])
plot(keepsBTCSP$rhoc[1:total])
plot(keepsBTCSP$xi_y2eta[1:total])
