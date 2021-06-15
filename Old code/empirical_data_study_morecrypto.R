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
sourceCpp("pgas_mstuart.cpp") #C++ updates
B <- 20000 #how many burn in draws to throw away
R <- 10000 #how many draws to keep after burn in
n_chns <- 4 #how many total chains to run
getSymbols("ETH-USD",from = "2016-01-01",to = "2021-01-31")
ETH <- as.data.frame(`ETH-USD`)
ETH$Date <- as.Date(rownames(ETH))
ETH$`ETH-USD.Close`[ETH$Date=="2020-10-13"] <- 381.19
ETH$`ETH-USD.Close`[ETH$Date=="2020-10-12"] <- 387.73
ETH$`ETH-USD.Close`[ETH$Date=="2020-10-09"] <- 356.59
ETH$`ETH-USD.Close`[ETH$Date=="2020-04-17"] <- 171.64

getSymbols("XRP-USD",from = "2016-01-01",to = "2021-01-31")
XRP <- as.data.frame(`XRP-USD`)
XRP$Date <- as.Date(rownames(XRP))
XRP$`XRP-USD.Close`[XRP$Date=="2020-10-13"] <- 0.2563
XRP$`XRP-USD.Close`[XRP$Date=="2020-10-12"] <- 0.2564
XRP$`XRP-USD.Close`[XRP$Date=="2020-10-09"] <- 0.2543
XRP$`XRP-USD.Close`[XRP$Date=="2020-04-17"] <- 0.1902

getSymbols("BTC-USD",from = "2016-01-01",to = "2021-01-31")
BTC <- as.data.frame(`BTC-USD`)
BTC$Date <- as.Date(rownames(BTC))
BTC$`BTC-USD.Close`[BTC$Date=="2020-10-13"] <- 11425.90
BTC$`BTC-USD.Close`[BTC$Date=="2020-10-12"] <- 11555.36
BTC$`BTC-USD.Close`[BTC$Date=="2020-10-09"] <- 11064.46
BTC$`BTC-USD.Close`[BTC$Date=="2020-04-17"] <- 7096.18

getSymbols("GILD",from = "2016-01-01",to = "2021-01-31")
GILD <- as.data.frame(GILD)
GILD$Date <- as.Date(rownames(GILD))

getSymbols("VRTX",from = "2016-01-01",to = "2021-01-31")
VRTX <- as.data.frame(VRTX)
VRTX$Date <- as.Date(rownames(VRTX))

S <- right_join(ETH,XRP) %>%
  right_join(BTC) %>%
  right_join(GILD) %>%
  right_join(VRTX)
T <- nrow(S) - 1

y <- 100*(log(S[-1,c("ETH-USD.Close")]) - log(S[-nrow(S),c("ETH-USD.Close")]))
x <- rep(0,length(y))
source("run_mcmc.R") 
#R+B iterations of pgas.R and pgas.cpp updates
#save results! not sure how to do this... 
#maybe save each keeps object as separate RDS file?
saveRDS(keeps,paste0("../keeps/keepsETH.rds"))

y <- 100*(log(S[-1,c("XRP-USD.Close")]) - log(S[-nrow(S),c("XRP-USD.Close")]))
x <- rep(0,length(y))
source("run_mcmc.R") 
#R+B iterations of pgas.R and pgas.cpp updates
#save results! not sure how to do this... 
#maybe save each keeps object as separate RDS file?
saveRDS(keeps,paste0("../keeps/keepsXRP.rds"))

y <- 100*(log(S[-1,c("BTC-USD.Close")]) - log(S[-nrow(S),c("BTC-USD.Close")]))
x <- rep(0,length(y))
source("run_mcmc.R") 
#R+B iterations of pgas.R and pgas.cpp updates
#save results! not sure how to do this... 
#maybe save each keeps object as separate RDS file?
saveRDS(keeps,paste0("../keeps/keepsBTC.rds"))

y <- 100*(log(S[-1,c("GILD.Close")]) - log(S[-nrow(S),c("GILD.Close")]))
x <- rep(0,length(y))
source("run_mcmc.R") 
#R+B iterations of pgas.R and pgas.cpp updates
#save results! not sure how to do this... 
#maybe save each keeps object as separate RDS file?
saveRDS(keeps,paste0("../keeps/keepsGILD.rds"))

y <- 100*(log(S[-1,c("VRTX.Close")]) - log(S[-nrow(S),c("VRTX.Close")]))
x <- rep(0,length(y))
source("run_mcmc.R") 
#R+B iterations of pgas.R and pgas.cpp updates
#save results! not sure how to do this... 
#maybe save each keeps object as separate RDS file?
saveRDS(keeps,paste0("../keeps/keepsVRTX.rds"))
