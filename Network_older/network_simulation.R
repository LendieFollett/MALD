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

users <- read.csv("Network_older/my-wallet-n-users")
trans <- read.csv("Network_older/n-transactions")
addy <- read.csv("Network_older/n-unique-addresses")

users$Timestamp <- as.Date(users$Timestamp)
trans$Timestamp <- as.Date(trans$Timestamp)
addy$Timestamp <- as.Date(addy$Timestamp)

users <- users %>% group_by(Timestamp)%>%summarise(users = mean(my.wallet.n.users)) %>%subset(Timestamp > as.Date("2011-11-01")) %>% mutate_at(c(2), function(x){x - lag(x)})
trans <- trans %>% group_by(Timestamp)%>%summarise(trans = mean(n.transactions)) %>%subset(Timestamp > as.Date("2011-11-01")) %>% mutate_at(c(2), function(x){x - lag(x)})
addy <- addy %>% group_by(Timestamp)%>%summarise(addy = mean(n.unique.addresses)) %>%subset(Timestamp > as.Date("2011-11-01")) %>% mutate_at(c(2), function(x){x - lag(x)})

delta_net <- users%>%merge(trans, all.x=TRUE, all.y=TRUE)%>%merge(addy, all.x=TRUE, all.y=TRUE) %>%
  subset(Timestamp < as.Date("2018-01-01"))

#delta_net <- net %>% mutate_at(c(2:4), function(x){x - lag(x)})

getSymbols("BTC-USD",from = min(delta_net$Timestamp),to =  max(delta_net$Timestamp))
BTC <- as.data.frame(`BTC-USD`)
BTC$Date <- as.Date(rownames(BTC))#seq(min(delta_net$Timestamp),max(delta_net$Timestamp),by="days")
#BTC$`BTC-USD.Close`[BTC$Date=="2020-04-17"] <- 7096.18
BTC_return_dat <- data.frame(Date = BTC$Date[-1], 
                             BTC = 100*(log(BTC[-1,c( "BTC-USD.Close") ])-
                                          log(BTC[-nrow(BTC),c("BTC-USD.Close") ])))

BTC_dat <- data.frame(Date = BTC$Date[-1], 
                      BTC = BTC[-1,c( "BTC-USD.Close") ])

all <- merge(delta_net, BTC_return_dat, by.x = "Timestamp", by.y= "Date", all=TRUE) %>%subset(Timestamp < as.Date("2015-01-01"))

cor(delta_net[,-1])

lm_trans <- lm(BTC ~ trans , data = all)
summary(lm_trans)
lm_addy <- lm(BTC ~ addy, data = all)
summary(lm_addy)
lm_users <- lm(BTC ~ users , data = all)
summary(lm_users)
lm5 <- lm(BTC ~  trans + addy + users, data = all)
summary(lm5)


