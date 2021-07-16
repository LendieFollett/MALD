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

delta_net <- net %>% mutate_at(c(2:4), function(x){x - lag(x)})

getSymbols("BTC-USD",from = min(delta_net$Timestamp),to =  max(delta_net$Timestamp))
BTC <- as.data.frame(`BTC-USD`)
BTC$Date <- seq(min(delta_net$Timestamp),max(delta_net$Timestamp),by="days")
BTC$`BTC-USD.Close`[BTC$Date=="2020-04-17"] <- 7096.18
BTC_return_dat <- data.frame(Date = BTC$Date[-1], 
                 BTC = 100*(log(BTC[-1,c( "BTC-USD.Close") ])-
       log(BTC[-nrow(BTC),c("BTC-USD.Close") ])))

BTC_dat <- data.frame(Date = BTC$Date[-1], 
                      BTC = BTC[-1,c( "BTC-USD.Close") ])

all <- merge(delta_net, BTC_return_dat, by.x = "Timestamp", by.y= "Date")


lm_pmts <- lm(BTC ~ pmts, data = all)
summary(lm_pmts)
lm_trans <- lm(BTC ~ trans , data = all)
summary(lm_trans)
lm_addy <- lm(BTC ~ addy, data = all)
summary(lm_addy)
lm5 <- lm(BTC ~  pmts + trans + addy, data = all)
summary(lm5)
lm_users <- lm(BTC ~ users , data = all)
##summary(lm_users)


ynew <- predict(lm5, all) + rnorm(nrow(all), 0, (sum(lm5$residuals^2)/lm5$df.residual)^0.5)
plot(ynew, type = "l")


