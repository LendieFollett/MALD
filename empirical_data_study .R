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
#setwd("/Users/000766412/Box Sync/ALD_Codes/mstuart/Empirical Study")
#################################################### 
# SAMPLE ALL PARAMETERS - FIXED AND LATENT ----------
#################################################### 
#set fix = FALSE to let theta parameters be sampled
fix <- FALSE

#sourceCpp("pgas_mstuart_2d.cpp") #C++ updates
#source("pgas_mstuart_2d.R") #R updates
#initialize values, create space to save draws
B <- 10000 #how many burn in draws to throw away
R <- 10000 #how many draws to keep after burn in
n_chns <- 1 #how many chains to run

getSymbols("BTC-USD",from = "2014-09-15",to = "2020-09-30")
BTC <- as.data.frame(`BTC-USD`)
BTC$Date <- seq(as.Date("2014-09-17"),as.Date("2020-09-30"),by="days")
BTC$`BTC-USD.Close`[BTC$Date=="2020-04-17"] <- 7096.18
getSymbols("^GSPC",from = "2014-09-15",to = "2020-09-30")
SP500 <- as.data.frame(`GSPC`)
SP500$Date <- as.Date(rownames(SP500))
S <- merge(BTC,SP500)
T <- nrow(S) - 1
S <- merge(BTC,SP500)

# #2-D MODEL MCMC
# y <- as.matrix(100*(log(S[-1,c("BTC-USD.Close","GSPC.Close")]) - log(S[-nrow(S),c("BTC-USD.Close","GSPC.Close")])))
# x <- array(0,dim=dim(y))
# #source("starting_values_2d.R") #initialize values
# source("run_mcmc_2d.R") #R+B iterations of pgas.R and pgas.cpp updates
# #save results! not sure how to do this... 
# #maybe save each keeps object as separate RDS file?
# saveRDS(keeps,paste0("keeps/keepsBTCSP.rds"))

#-----1 D MODEL MCMC------
sourceCpp("pgas_mstuart.cpp") #C++ updates
y <- 100*(log(S[-1,c("BTC-USD.Close")]) - log(S[-nrow(S),c("BTC-USD.Close")]))
x <- rep(0,length(y))
#source("starting_values.R") #initialize values
source("run_mcmc.R") 
#R+B iterations of pgas.R and pgas.cpp updates
#save results! not sure how to do this... 
#maybe save each keeps object as separate RDS file?
saveRDS(keeps,paste0("../keeps/keepsBTC.rds"))
y <- 100*(log(S[-1,c("GSPC.Close")]) - log(S[-nrow(S),c("GSPC.Close")]))
x <- rep(0,length(y))
#source("starting_values.R") #initialize values
source("run_mcmc.R") 
#R+B iterations of pgas.R and pgas.cpp updates
#save results! not sure how to do this... 
#maybe save each keeps object as separate RDS file?
saveRDS(keeps,paste0("../keeps/keepsSP.rds"))


#################################################### 
# SAMPLE ONLY LATENT - FIX THETA-HAT ----------
#################################################### 
#keepsBTCSP <- readRDS("keepsBTCSP.rds")
keepsBTC <- readRDS("../keeps/keepsBTC.rds")
keepsSP <- readRDS("../keeps/keepsSP.rds")

# #store MAP estimates
# THETA_map2d <- list()
# 
# THETA_map2d$lambda<-  keepsBTCSP$lambda %>%  apply(2, mean) %>% (function(x){x/sum(x)})
# THETA_map2d$sigma_v <- keepsBTCSP$sigma_v %>% apply(2, mean)
# THETA_map2d$sigma_c <- keepsBTCSP$sigma_c %>% apply(2, mean)
# THETA_map2d$rhoc <- keepsBTCSP$rhoc %>%mean
# THETA_map2d$phi <- keepsBTCSP$phi %>%apply(2, mean)
# THETA_map2d$theta<- keepsBTCSP$theta %>%apply(2, mean)
# THETA_map2d$mu <- keepsBTCSP$mu %>%apply(2, mean)
# THETA_map2d$rho<- keepsBTCSP$rho %>%apply(2, mean)
# THETA_map2d$xi_yeta <- cbind(keepsBTCSP$xi_y1eta,keepsBTCSP$xi_y2eta )%>%apply(2, mean)
# THETA_map2d$xi_yw <- cbind(keepsBTCSP$xi_y1w,keepsBTCSP$xi_y2w )%>%apply(2, mean)

#1d model
THETA_map1d <- list()

THETA_map1d$lambda<-  cbind(keepsBTC$lambda[1],keepsSP$lambda[1]) %>%  apply(2, mean) 
THETA_map1d$sigma_v <- cbind(keepsBTC$sigma_v,keepsSP$sigma_v) %>% apply(2, mean)
THETA_map1d$phi <-cbind(keepsBTC$phi,keepsSP$phi) %>%apply(2, mean)
THETA_map1d$theta<- cbind(keepsBTC$theta,keepsSP$theta) %>%apply(2, mean)
THETA_map1d$mu <- cbind(keepsBTC$mu,keepsSP$mu) %>%apply(2, mean)
THETA_map1d$rho<- cbind(keepsBTC$rho,keepsSP$rho) %>%apply(2, mean)
THETA_map1d$xi_yeta <- cbind(keepsBTC$xi_yeta,keepsSP$xi_yeta) %>%apply(2, mean)
THETA_map1d$xi_yw<- cbind(keepsBTC$xi_yw,keepsSP$xi_yw) %>%apply(2, mean)


#########1-d case 
#FIX = TRUE means we're fixing the parameters
fix <- TRUE
MAPs <- sapply(THETA_map1d, "[[", 1) %>%as.list() #where to find MAP estimates

#BTC
sourceCpp("pgas_mstuart.cpp") #C++ updates
y <- 100*(log(S[-1,c("BTC-USD.Close")]) - log(S[-nrow(S),c("BTC-USD.Close")]))
x <- rep(0,length(y))
#source("starting_values.R") #initialize values
source("run_mcmc.R") 
saveRDS(keeps,paste0("../keeps/keepsBTC_MAP.rds"))

#S&P
fix <- TRUE
MAPs <- sapply(THETA_map1d, "[[", 2) %>%as.list() #where to find MAP estimates
y <- 100*(log(S[-1,c("GSPC.Close")]) - log(S[-nrow(S),c("GSPC.Close")]))
x <- rep(0,length(y))
#source("starting_values.R") #initialize values
source("run_mcmc.R") 
#R+B iterations of pgas.R and pgas.cpp updates
#save results! not sure how to do this... 
#maybe save each keeps object as separate RDS file?
saveRDS(keeps,paste0("../keeps/keepsSP_MAP.rds"))

# #########2-d case 
# y <- as.matrix(100*(log(S[-1,c("BTC-USD.Close","GSPC.Close")]) - log(S[-nrow(S),c("BTC-USD.Close","GSPC.Close")])))
# x <- array(0,dim=dim(y))
# fix <- TRUE
# MAPs <- THETA_map2d #where to find MAP estimates
# #source("starting_values_2d.R") #initialize values
# source("run_mcmc_2d.R") #R+B iterations of pgas.R and pgas.cpp updates
# #save results! not sure how to do this... 
# #maybe save each keeps object as separate RDS file?
# saveRDS(keeps,paste0("keeps/keepsBTCSP.rds"))

