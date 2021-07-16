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

net <- users %>%merge(pmts)%>%merge(trans)%>%merge(addy)




