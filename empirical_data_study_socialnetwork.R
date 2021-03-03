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

fix <- FALSE
sourceCpp("pgas.cpp") #C++ updates
B <- 20000 #how many burn in draws to throw away
R <- 10000 #how many draws to keep after burn in
n_chns <- 4 #how many total chains to run
getSymbols("SNAP",from = "2017-05-01",to = "2021-01-31")
SNAP <- as.data.frame(SNAP)
SNAP$Date <- as.Date(rownames(SNAP))

getSymbols("TWTR",from = "2014-05-01",to = "2018-01-31")
TWTR <- as.data.frame(TWTR)
TWTR$Date <- as.Date(rownames(TWTR))

getSymbols("YELP",from = "2013-05-01",to = "2017-01-31")
YELP <- as.data.frame(YELP)
YELP$Date <- as.Date(rownames(YELP))

getSymbols("MTCH",from = "2015-05-01",to = "2019-01-31")
MTCH <- as.data.frame(MTCH)
MTCH$Date <- as.Date(rownames(MTCH))

y <- 100*(log(TWTR[-1,c("TWTR.Close")]) - log(TWTR[-nrow(TWTR),c("TWTR.Close")]))
x <- rep(0,length(y))
source("run_mcmc.R") 
#R+B iterations of pgas.R and pgas.cpp updates
#save results! not sure how to do this... 
#maybe save each keeps object as separate RDS file?
saveRDS(keeps,paste0("keepsTWTR.rds"))

y <- 100*(log(SNAP[-1,c("SNAP.Close")]) - log(SNAP[-nrow(SNAP),c("SNAP.Close")]))
x <- rep(0,length(y))
source("run_mcmc.R") 
#R+B iterations of pgas.R and pgas.cpp updates
#save results! not sure how to do this... 
#maybe save each keeps object as separate RDS file?
saveRDS(keeps,paste0("keepsSNAP.rds"))

y <- 100*(log(YELP[-1,c("YELP.Close")]) - log(YELP[-nrow(YELP),c("YELP.Close")]))
x <- rep(0,length(y))
source("run_mcmc.R") 
#R+B iterations of pgas.R and pgas.cpp updates
#save results! not sure how to do this... 
#maybe save each keeps object as separate RDS file?
saveRDS(keeps,paste0("keepsYELP.rds"))

y <- 100*(log(MTCH[-1,c("MTCH.Close")]) - log(MTCH[-nrow(MTCH),c("MTCH.Close")]))
x <- rep(0,length(y))
source("run_mcmc.R") 
#R+B iterations of pgas.R and pgas.cpp updates
#save results! not sure how to do this... 
#maybe save each keeps object as separate RDS file?
saveRDS(keeps,paste0("keepsMTCH.rds"))

