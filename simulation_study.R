rm(list = ls())  ## DON'T FORGET TO SET WD!!
library(ald)
library(ggplot2)
library(dplyr)
library(truncnorm)
library(tmvtnorm)
library(Rcpp)
library(RcppTN)
library(MCMCpack)

sourceCpp("pgas_mstuart.cpp") #C++ updates
#source("pgas_mstuart.R") #R updates
#initialize values, create space to save draws
n_chns <- 1
fix <- FALSE
B <- 5000 #how many burn in draws to throw away
R <- 5000 #how many draws to keep after burn in

for (s in 1:100){
  source("sim_data_alt_param.R") #simulate data based on simulation id s
  Y = y
  source("starting_values.R") #initialize values
  source("run_mcmc.R") #R+B iterations of pgas.R and pgas.cpp updates
  #save results! not sure how to do this... 
  #maybe save each keeps object as separate RDS file?
  saveRDS(keeps,paste0("keeps/keeps",s,".rds"))
}


