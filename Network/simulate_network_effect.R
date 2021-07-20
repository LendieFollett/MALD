rm(list = ls())  
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
library(lubridate)

users <- read.csv("Network/my-wallet-n-users.txt")
pmts <- read.csv("Network/n-payments.txt")
trans <- read.csv("Network/n-transactions.txt")
addy <- read.csv("Network/n-unique-addresses.txt")

users$Timestamp <- as.Date(users$Timestamp)
pmts$Timestamp <- as.Date(pmts$Timestamp)
trans$Timestamp <- as.Date(trans$Timestamp)
addy$Timestamp <- as.Date(addy$Timestamp)

users <- users %>% group_by(Timestamp)%>%summarise(users = mean(my.wallet.n.users))
pmts <- pmts %>% group_by(Timestamp)%>%summarise(pmts = mean(n.payments))
trans <- trans %>% group_by(Timestamp)%>%summarise(trans = mean(n.transactions))
addy <- addy %>% group_by(Timestamp)%>%summarise(addy = mean(n.unique.addresses))

net <- users %>%merge(pmts,all=TRUE)%>%merge(trans,all=TRUE)%>%merge(addy,all=TRUE)

delta_net <- net %>% mutate_at(c(2:5), function(x){x - lag(x)})

getSymbols("BTC-USD",from = min(delta_net$Timestamp),to =  max(delta_net$Timestamp))
BTC <- as.data.frame(`BTC-USD`)
BTC$Date <- seq(min(delta_net$Timestamp),max(delta_net$Timestamp),by="days")
BTC$`BTC-USD.Close`[BTC$Date=="2020-04-17"] <- 7096.18
BTC_return_dat <-  data.frame(Date = BTC$Date[-1], 
                              BTC = 100*(log(BTC[-1,c( "BTC-USD.Close") ])-
                                       log(BTC[-nrow(BTC),c("BTC-USD.Close") ])))

BTC_dat <- data.frame(Date = BTC$Date[-1], 
                      BTC = BTC[-1,c( "BTC-USD.Close") ])

all <- merge(delta_net, BTC_return_dat, by.x = "Timestamp", by.y= "Date")

lm_users <- lm(BTC ~ users , data = all)
summary(lm_users)
lm_pmts <- lm(BTC ~ pmts, data = all)
summary(lm_pmts)
lm_trans <- lm(BTC ~ trans , data = all)
summary(lm_trans)
lm_addy <- lm(BTC ~ addy, data = all)
summary(lm_addy)
lm5 <- lm(BTC ~  pmts + trans + addy + users, data = all)
summary(lm5)



ynew <- predict(lm5, all) 
ynew[is.na(ynew)] <- coef(lm5)[1]
ynew <- ynew + rnorm(nrow(all), 0, (sum(lm5$residuals^2)/lm5$df.residual)^0.5)

plot(BTC_return_dat$BTC, type = "l")
plot(ynew, type  ="l")
lines(predict(lm5, all), col = "red")
all$ynew <- ynew


#------fit model

#################################################### 
# SVMALD MODEL ----------
#################################################### 
thin <- 5 #thinning param
B <- 10000 #how many burn in draws to throw away
R <- 100000 #how many draws to keep after burn in
n_chns <- 1 #how many chains to run
use_starting_values <- FALSE
sourceCpp("pgas_2d.cpp") #C++ updates
# #2-D MODEL MCMCb        
getSymbols("^GSPC",from = min(all$Timestamp),to = max(all$Timestamp))
SP500 <- as.data.frame(`GSPC`)
SP500$Date <- as.Date(rownames(SP500))
S <- merge(all[,c("Timestamp", "ynew")],SP500, by.x = "Timestamp", by.y = "Date")
T <- nrow(S) - 1
y <- cbind(S$ynew[-1], 100*(log(S[-1,c("GSPC.Close")]) - log(S[-nrow(S),c("GSPC.Close")])))

yprim <- array(0,dim=dim(y))
#source("starting_values_2d.R") #initialize values (performed within run_mcmc_2d.R)
exp_jumps <- norm_jumps <- ind <- FALSE
source("run_mcmc_2d.R") #R+B iterations of pgas.R and pgas.cpp updates
saveRDS(keeps,paste0("Network/keepsBTCSPSIM_MALD.rds"))

apply(keeps$rho, 2, mean)
#[1]  0.003716513 -0.020089649 -0.005576129 -0.6242593
mean(keeps$rho[,3]>0)
#[1] 0.4986

