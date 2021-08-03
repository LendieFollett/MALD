rm(list = ls())  ## DON'T FORGET TO SET WD!!
library(ald)
library(ggplot2)
library(dplyr)
library(truncnorm)
library(mvtnorm)
library(tmvtnorm)
library(Rcpp)
library(MCMCpack)

sourceCpp("pgas_2d.cpp") #C++ updates
#source("pgas_mstuart_2d.R") #R updates
#initialize values, create space to save draws
n_chns <- 1
thin <- 5
fix <- FALSE
B <- 10000 #how many burn in draws to throw away
R <- 100000 #how many draws to keep after burn in
use_starting_values <- FALSE

for (s in 1:10){
  exp_jumps <- norm_jumps <- ind <- FALSE
  source("sim_data_alt_param_2d.R") #simulate data based on simulation id s
  Y = y
  yprim = matrix(0,nrow=T,ncol=2)
  #source("starting_values_2d.R") #initialize values
  source("run_mcmc_2d.R") #R+B iterations of pgas.R and pgas.cpp updates
  #save results! not sure how to do this... 
  #maybe save each keeps object as separate RDS file?
  saveRDS(keeps,paste0("keeps/keeps",s,".rds"))
}


