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

users2 <- users %>% group_by(Timestamp)%>%summarise(users = mean(my.wallet.n.users)) %>%subset(Timestamp > as.Date("2011-11-01")) %>%mutate(month = format(Timestamp, "%m%Y"))%>%group_by(month)%>%summarise(users = mean(users))%>%mutate_at(c(2), function(x){x - lag(x)}) %>%mutate(month = as.Date(paste0("01", month), format = "%d%m%Y"))
trans2 <- trans %>% group_by(Timestamp)%>%summarise(trans = mean(n.transactions)) %>%subset(Timestamp > as.Date("2011-11-01")) %>%mutate(month = format(Timestamp, "%m%Y"))%>%group_by(month)%>%summarise(trans = mean(trans)) %>% mutate_at(c(2), function(x){x - lag(x)})%>%mutate(month = as.Date(paste0("01", month), format = "%d%m%Y"))
addy2 <- addy %>% group_by(Timestamp)%>%summarise(addy = mean(n.unique.addresses)) %>%subset(Timestamp > as.Date("2011-11-01")) %>%mutate(month = format(Timestamp, "%m%Y"))%>%group_by(month)%>%summarise(addy = mean(addy))%>% mutate_at(c(2), function(x){x - lag(x)})%>%mutate(month = as.Date(paste0("01", month), format = "%d%m%Y"))

delta_net <- users2%>%merge(trans2, all.x=TRUE, all.y=TRUE)%>%merge(addy2, all.x=TRUE, all.y=TRUE) %>%
  arrange(month)

#delta_net <- net %>% mutate_at(c(2:4), function(x){x - lag(x)})

getSymbols("BTC-USD",from = min(users$Timestamp),to =  max(users$Timestamp))
BTC <- as.data.frame(`BTC-USD`)
BTC$Date <- as.Date(rownames(BTC))#seq(min(delta_net$Timestamp),max(delta_net$Timestamp),by="days")
#BTC$`BTC-USD.Close`[BTC$Date=="2020-04-17"] <- 7096.18
BTC <- BTC%>%mutate(month = format(Date, "%m%Y"))%>%group_by(month)%>%summarise(`BTC-USD.Close` = mean(`BTC-USD.Close`))%>%
  arrange(month)
BTC_return_dat <- data.frame(Date = BTC$month[-1], 
                             BTC = 100*((BTC[-1,c( "BTC-USD.Close") ])-
                                          (BTC[-nrow(BTC),c("BTC-USD.Close") ]))/BTC[-nrow(BTC),c("BTC-USD.Close") ])%>%
  mutate(month = as.Date(paste0("01", Date), format = "%d%m%Y"))

BTC_dat <- data.frame(Date = BTC$month[-1], 
                      BTC = BTC[-1,c( "BTC-USD.Close") ])%>%
  arrange(Date)

all <- merge(delta_net, BTC_return_dat, by.x = "month", by.y= "month", all=TRUE) #%>%subset(Timestamp > as.Date("2015-01-01"))

cor(delta_net[,-1])

lm_trans <- lm(BTC.USD.Close ~ trans , data = all)
summary(lm_trans)
lm_addy <- lm(BTC.USD.Close ~ addy, data = all)
summary(lm_addy)
lm_users <- lm(BTC.USD.Close ~ users , data = all)
summary(lm_users)
lm5 <- lm(BTC.USD.Close ~  trans + addy+users , data = all)
summary(lm5)


