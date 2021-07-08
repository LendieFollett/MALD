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

###Data frame of model parameters
models <- data.frame(exp_jumps =  c(FALSE,   TRUE,  FALSE,     FALSE),
                     norm_jumps = c(FALSE,   FALSE, TRUE,      FALSE),
                     ind =        c(FALSE,   FALSE, FALSE,     TRUE),
                     model =      c("SVMALD", "SVLD", "SVMVN", "SVIND"))

##### LONG TIME PERIOD #####
fix <- FALSE
thin <- 5 #thinning param
B <- 10000 #how many burn in draws to throw away
R <- 100000 #how many draws to keep after burn in
n_chns <- 1 #how many chains to run
getSymbols("NVAX",from = "2014-09-15",to = "2020-09-30")
NVAX <- as.data.frame(NVAX)
NVAX$Date <- as.Date(rownames(NVAX))
getSymbols("^GSPC",from = "2014-09-15",to = "2020-09-30")
SP500 <- as.data.frame(`GSPC`)
SP500$Date <- as.Date(rownames(SP500))
S <- merge(NVAX,SP500)
T <- nrow(S) - 1
for (k in 1:nrow(models)){
  print(paste0("----- > Starting ", models$model[k], " model < -------"))
  use_starting_values <- FALSE
  sourceCpp("pgas_2d.cpp") #C++ updates
  # #2-D MODEL MCMCb        cfv09
  y <- as.matrix(100*(log(S[-1,c("NVAX.Close","GSPC.Close")]) - log(S[-nrow(S),c("NVAX.Close","GSPC.Close")])))
  yprim <- array(0,dim=dim(y))
  exp_jumps <- models$exp_jumps[k]
  norm_jumps <- models$norm_jumps[k]
  ind <- models$ind[k]
  source("run_mcmc_2d.R") #R+B iterations of pgas.R and pgas.cpp updates
  saveRDS(keeps,paste0("keeps_",models$model[k] ,"_NVAX_long.rds"))
}

##### SHORT TIME PERIOD #####
fix <- FALSE
thin <- 5 #thinning param
B <- 10000 #how many burn in draws to throw away
R <- 100000 #how many draws to keep after burn in
n_chns <- 1 #how many chains to run
getSymbols("NVAX",from = "2020-10-01",to = "2021-06-30")
NVAX <- as.data.frame(NVAX)
NVAX$Date <- as.Date(rownames(NVAX))
getSymbols("^GSPC",from = "2020-10-01",to = "2021-06-30")
SP500 <- as.data.frame(`GSPC`)
SP500$Date <- as.Date(rownames(SP500))
S <- merge(NVAX,SP500)
T <- nrow(S) - 1
for (k in 1:nrow(models)){
  print(paste0("----- > Starting ", models$model[k], " model < -------"))
  use_starting_values <- FALSE
  sourceCpp("pgas_2d.cpp") #C++ updates
  # #2-D MODEL MCMCb        cfv09
  y <- as.matrix(100*(log(S[-1,c("NVAX.Close","GSPC.Close")]) - log(S[-nrow(S),c("NVAX.Close","GSPC.Close")])))
  yprim <- array(0,dim=dim(y))
  exp_jumps <- models$exp_jumps[k]
  norm_jumps <- models$norm_jumps[k]
  ind <- models$ind[k]
  source("run_mcmc_2d.R") #R+B iterations of pgas.R and pgas.cpp updates
  saveRDS(keeps,paste0("keeps_",models$model[k] ,"_NVAX_short.rds"))
}