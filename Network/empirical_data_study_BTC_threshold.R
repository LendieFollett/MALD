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
# getSymbols("BTC-USD",from = "2020-10-01",to = "2021-05-31")
# BTC <- as.data.frame(`BTC-USD`)
# BTC$Date <- seq(as.Date("2020-10-01"),as.Date("2021-05-31"),by="days")
# BTC$`BTC-USD.Close`[BTC$Date=="2020-04-17"] <- 7096.18
BTC <- read.csv("BTC_USD_2014-11-03_2022-01-12-CoinDesk.csv")
BTC <- BTC %>% dplyr::select(!Currency)
names(BTC) <- c("Date","BTC-USD.Close","BTC-USD.Open","BTC-USD.High","BTC-USD.Low")
BTC$Date <- as.Date(BTC$Date)
getSymbols("GME",from = "2020-10-01",to = "2021-12-31")
GME <- as.data.frame(GME)
GME$Date <- as.Date(rownames(GME))
getSymbols("AMC",from = "2020-10-01",to = "2021-12-31")
AMC <- as.data.frame(AMC)
AMC$Date <- as.Date(rownames(AMC))
# getSymbols("DOGE-USD",from = "2020-10-01",to = "2021-05-31")
# DOGE <- as.data.frame(`DOGE-USD`)
# DOGE$Date <- as.Date(rownames(DOGE))
DOGE <- read.csv("DOGE_USD_2019-02-27_2022-01-12-CoinDesk.csv")
DOGE <- DOGE %>% dplyr::select(!Currency)
DOGE$Date <- as.Date(DOGE$Date)
names(DOGE) <- c("Date","DOGE-USD.Close","DOGE-USD.Open","DOGE-USD.High","DOGE-USD.Low")
getSymbols("^GSPC",from = "2020-10-01",to = "2021-12-31")
SP500 <- as.data.frame(`GSPC`)
SP500$Date <- as.Date(rownames(SP500))
getSymbols("MRNA",from = "2020-10-01",to = "2021-12-31")
MRNA <- as.data.frame(MRNA)
MRNA$Date <- as.Date(rownames(MRNA))

S <- BTC %>% merge(GME) %>% merge(AMC) %>% merge(DOGE) %>% merge(SP500) %>% merge(MRNA)
T <- nrow(S) - 1


###Data frame of model parameters
models <- data.frame(exp_jumps =  c(FALSE,   TRUE,  FALSE,     FALSE),
                     norm_jumps = c(FALSE,   FALSE, TRUE,      FALSE),
                     ind =        c(FALSE,   FALSE, FALSE,     TRUE),
                     model =      c("SVMALD", "SVLD", "SVMVN", "SVIND"))

#################################################### 
# ALL MODELS ---------- BTC
#################################################### 
for (k in 1:1){
  print(paste0("----- > Starting ", models$model[k], " model < -------"))
keeps <- readRDS(paste0("keeps/keeps_","SVMALD" ,"_","BTC", ".rds"))
thresholds <- sqrt(252*apply(keeps$v[1:20000,,1], 2, mean)) %>% quantile(probs=seq(0.15,0.25,0.01))
print(thresholds)
  use_starting_values <- FALSE
  sourceCpp("pgas_2d_threshold.cpp") #C++ updates
  # #2-D MODEL MCMCb        cfv09
  y <- as.matrix(100*(log(S[-1,c("BTC-USD.Close","GSPC.Close")]) - log(S[-nrow(S),c("BTC-USD.Close","GSPC.Close")])))
  yprim <- array(0,dim=dim(y))
  exp_jumps <- models$exp_jumps[k]
  norm_jumps <- models$norm_jumps[k]
  ind <- models$ind[k]
  for (threshold in thresholds){
  source("run_mcmc_2d_threshold.R") #R+B iterations of pgas.R and pgas.cpp updates
  saveRDS(keeps,paste0("keeps/keeps_",models$model[k] ,"_BTC_",threshold,".rds"))
}
}
#ALL OF THE ABOVE WAS FOR THE 'SHORT TIME PERIOD' ANALYSIS


#################################################### 
# CONVERGENCE CHECKS ----------
#################################################### 
#library(LaplacesDemon)

#total <- 20000 #number of mcmc iterations saved after burn-in, thinning
#doESS <- function(x, total){
#  R <- total
#  if(!is.null(dim(x))){ #if it's a data frame
#    return(apply(x[1:R,], 2, ESS))
#  }else{
#    return(ESS(x[1:R]))
#  }
#}


#SVMALD
#lapply(keeps[c(4,6:17)], doESS, total = 20000) %>% str()


#plot(keeps$sigma_c[,1], type = "l")
#plot(keeps$sigma_c[,2], type = "l")
#plot(keeps$rhoc, type = "l")
#plot(keeps$xi_cw[,1], type = "l")
#plot(keeps$xi_cw[,2], type = "l")
#plot(keeps$xi_y1eta, type = "l")
#plot(keeps$xi_y2eta, type = "l")
#plot(keeps$xi_y1w, type = "l")
#plot(keeps$xi_y2w, type = "l")
#SVALD
#lapply(keepsIND[c(4,6:17)], doESS, total = total) %>% str()
#SVMVN
#lapply(keepsBTCSP_MVN[c(4,6:17)], doESS, total = total) %>% str()
#SVLD
#lapply(keepsBTCSP_LD[c(4,6:17)], doESS, total = total) %>% str()



#plot(keepsBTCSP$sigma_c[1:total, 1]);length(unique(keepsBTCSP$sigma_c[1:total, 1]))/total
#plot(keepsBTCSP$sigma_c[1:total, 1])
#plot(keepsBTCSP$rhoc[1:total])
#plot(keepsBTCSP$xi_y2eta[1:total])

#use for starting values?
#starting_values <- lapply(keeps, domean, total = 20000)
#starting_values %>% str
#starting_values$delta <- apply(starting_values$delta , 2, round)
#starting_values$xi_y1 <- starting_values$xi_y1c <- starting_values$J[,1] + rnorm(length(starting_values$xi_y1), 0, .001)
#starting_values$xi_y2 <- starting_values$xi_y2c <- starting_values$J[,2]+ rnorm(length(starting_values$xi_y1), 0, .001)
#saveRDS(starting_values,"starting_values_MALD.rds")

