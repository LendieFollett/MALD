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
library(xtable)


thin <- 5 #thinning param
B <- 10000 #how many burn in draws to throw away
R <- 100000 #how many draws to keep after burn in
n_chns <- 1 #how many chains to run
#load data

getSymbols("^GSPC",from = "2020-10-01",to = "2021-06-30")
SP500 <- as.data.frame(`GSPC`)
SP500$Date <- as.Date(rownames(SP500))
getSymbols("DIS",from = "2020-10-01",to = "2021-06-30")
DIS <- as.data.frame(DIS)
DIS$Date <- as.Date(rownames(DIS))
getSymbols("BBY",from = "2020-10-01",to = "2021-06-30")
BBY <- as.data.frame(BBY)
BBY$Date <- as.Date(rownames(BBY))
getSymbols("PFE",from = "2020-10-01",to = "2021-06-30")
PFE <- as.data.frame(PFE)
PFE$Date <- as.Date(rownames(PFE))

S <- DIS %>% merge(BBY) %>% merge(PFE) %>% merge(SP500)
T <- nrow(S) - 1


######FIGURE 1
# Date <- S$Date[-1]
# cbind(100*(log(S[-1,c("GSPC.Close", "GME.Close", "BTC-USD.Close", "AMC.Close", "DOGE-USD.Close") ])-
#       log(S[-nrow(S),c("GSPC.Close", "GME.Close", "BTC-USD.Close", "AMC.Close", "DOGE-USD.Close") ])),
#   Date) %>%
#   melt(id.vars = c("Date")) %>%
#   mutate(variable = factor(variable, levels = c("GSPC.Close", "GME.Close", "BTC-USD.Close", "AMC.Close", "DOGE-USD.Close"),
#                            labels = c("S&P", "GME", "BTC", "AMC", "DOGE"))) %>%
#   ggplot() +
#   geom_line(aes(x = Date, y = value)) +
#   facet_grid(variable~., scales = "free_y") +
#   theme_bw() +
#   scale_x_date(breaks = "month")
# 
# ggsave("data_plot_meme.pdf", height = 10, width = 8)

###Data frame of model parameters
models <- data.frame(exp_jumps =  c(FALSE,   TRUE,  FALSE,     FALSE),
                     norm_jumps = c(FALSE,   FALSE, TRUE,      FALSE),
                     ind =        c(FALSE,   FALSE, FALSE,     TRUE),
                     model =      c("SVMALD", "SVLD", "SVMVN", "SVIND"))


#################################################### 
# ALL MODELS ---------- Disney
#################################################### 
for (k in 1:nrow(models)){
  print(paste0("----- > Starting ", models$model[k], " model < -------"))
use_starting_values <- FALSE
sourceCpp("pgas_2d.cpp") #C++ updates
# #2-D MODEL MCMCb        cfv09
y <- as.matrix(100*(log(S[-1,c("DIS.Close","GSPC.Close")]) - log(S[-nrow(S),c("DIS.Close","GSPC.Close")])))
yprim <- array(0,dim=dim(y))
exp_jumps <- models$exp_jumps[k]
norm_jumps <- models$norm_jumps[k]
ind <- models$ind[k]
source("run_mcmc_2d.R") #R+B iterations of pgas.R and pgas.cpp updates
saveRDS(keeps,paste0("keeps/keeps_",models$model[k] ,"_DIS.rds"))
}

#################################################### 
# ALL MODELS ---------- Best Buy
#################################################### 
for (k in 1:nrow(models)){
  print(paste0("----- > Starting ", models$model[k], " model < -------"))
use_starting_values <- FALSE
sourceCpp("pgas_2d.cpp") #C++ updates
# #2-D MODEL MCMCb        cfv09
y <- as.matrix(100*(log(S[-1,c("BBY.Close","GSPC.Close")]) - log(S[-nrow(S),c("BBY.Close","GSPC.Close")])))
yprim <- array(0,dim=dim(y))
exp_jumps <- models$exp_jumps[k]
norm_jumps <- models$norm_jumps[k]
ind <- models$ind[k]
source("run_mcmc_2d.R") #R+B iterations of pgas.R and pgas.cpp updates
saveRDS(keeps,paste0("keeps/keeps_",models$model[k] ,"_BBY.rds"))
}
#################################################### 
# ALL MODELS ---------- Pfizer
#################################################### 
for (k in 1:nrow(models)){
  print(paste0("----- > Starting ", models$model[k], " model < -------"))
use_starting_values <- FALSE
sourceCpp("pgas_2d.cpp") #C++ updates
# #2-D MODEL MCMCb        cfv09
y <- as.matrix(100*(log(S[-1,c("PFE.Close","GSPC.Close")]) - log(S[-nrow(S),c("PFE.Close","GSPC.Close")])))
yprim <- array(0,dim=dim(y))
exp_jumps <- models$exp_jumps[k]
norm_jumps <- models$norm_jumps[k]
ind <- models$ind[k]
source("run_mcmc_2d.R") #R+B iterations of pgas.R and pgas.cpp updates
saveRDS(keeps,paste0("keeps/keeps_",models$model[k] ,"_PFE.rds"))
}

#ALL OF THE ABOVE WAS FOR THE 'SHORT TIME PERIOD' ANALYSIS


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
#starting_values <- lapply(keeps, domean, total = 20000)
#starting_values %>% str
#starting_values$delta <- apply(starting_values$delta , 2, round)
#starting_values$xi_y1 <- starting_values$xi_y1c <- starting_values$J[,1] + rnorm(length(starting_values$xi_y1), 0, .001)
#starting_values$xi_y2 <- starting_values$xi_y2c <- starting_values$J[,2]+ rnorm(length(starting_values$xi_y1), 0, .001)
#saveRDS(starting_values,"starting_values_MALD.rds")

