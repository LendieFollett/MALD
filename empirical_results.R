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
library(xtable)

#load data
getSymbols("BTC-USD",from = "2020-12-01",to = "2021-05-31")
BTC <- as.data.frame(`BTC-USD`)
BTC$Date <- seq(as.Date("2020-12-01"),as.Date("2021-05-31"),by="days")
BTC$`BTC-USD.Close`[BTC$Date=="2020-04-17"] <- 7096.18
getSymbols("GME",from = "2020-12-01",to = "2021-05-31")
GME <- as.data.frame(GME)
GME$Date <- as.Date(rownames(GME))
getSymbols("AMC",from = "2020-12-01",to = "2021-05-31")
AMC <- as.data.frame(AMC)
AMC$Date <- as.Date(rownames(AMC))
getSymbols("DOGE-USD",from = "2020-12-01",to = "2021-05-31")
DOGE <- as.data.frame(`DOGE-USD`)
DOGE$Date <- as.Date(rownames(DOGE))
getSymbols("^GSPC",from = "2020-12-01",to = "2021-05-31")
SP500 <- as.data.frame(`GSPC`)
SP500$Date <- as.Date(rownames(SP500))

S <- BTC %>%merge(GME) %>% merge(AMC) %>% merge(DOGE) %>% merge(SP500)
T <- nrow(S) - 1
Date <- S$Date[-1]

# MEANS  ----------


domean<- function(x, total){
  R <- total
  if(length(dim(x)) == 2){ #if it's a data frame
    return(apply(x[1:R,], 2, median))
  }else if (length(dim(x)) > 2){
    return(apply(x[1:R,,], 2:3, function(x){(median(x))}))
  }else{
    return(median(x[1:R]))
  }
}


###########################
#FOR ANALYSIS OF SHORT TERM TIME SERIES CRYPTO, MEME STOCKS
###########################

keeps <- readRDS(paste0("keeps_short/keeps_",SVMALD ,"_","BTC", ".rds"))
names <- lapply(keeps[c(4,6:17)], domean, total = 10) %>%unlist %>%names

keeps_summary <- array(dim = c(16, length(names)))
keeps_v1 <- array(dim = c(16,dim(keeps$v)[2]))
keeps_v2 <- array(dim = c(16,dim(keeps$v)[2]))
model <- rep(NA, 16)
data <- rep(NA, 16)
colnames(keeps_summary) <- names

j = 0
for (m in c("SVMALD", "SVMVN", "SVLD", "SVIND")){
  for (i in c("BTC", "DOGE", "AMC", "GMC")){
    j = j + 1
    keeps <- readRDS(paste0("keeps_short/keeps_",m ,"_",i, ".rds"))
    keeps_summary[j,] <- lapply(keeps[c(4,6:17)], domean, total = 20000) %>%unlist
    keeps_v1[j,] <- apply(keeps$v[1:R,,1], 2:3, function(x){(mean(x))}) #alternative sv
    keeps_v2[j,] <- apply(keeps$v[1:R,,2], 2:3, function(x){(mean(x))}) #sp sv
    model[j] <-m
    data[j] <- i
  }
}

keeps_summary <- keeps_summary %>%as.data.frame()%>%
  mutate_all(round, digits = 2) %>%
  mutate(lambda = paste0("(",lambda1,", ", lambda2,", ", lambda3,", ", lambda4,")"))%>%
  dplyr::select(-c(lambda1, lambda2, lambda3, lambda4))
keeps_summary$model <- model; keeps_summary$series <- data

#TABLE XXX POSTERIOR MEANS OF PARAMETERS------
keeps_summary[,c(ncol(keeps_summary), ncol(keeps_summary)-1, 1:(ncol(keeps_summary)-2))] %>%
  arrange(series) %>%
  xtable() %>%
  print()

#FIGURE XXX STOCHASTIC VOLATILITY--------
keeps_v1 %>%
  mutate(model = model, series = data) %>%
  melt(id.vars = c(model, series)) %>%
  mutate(Date = Date)%>% #CORRECT FOR ORDERING!
  ggplot() +
  geom_line(aes(x = Var2, y = value, colour = model)) +
  facet_grid(series~., scales ="free_y")


###########################
#FOR ANALYSIS OF LONG TERM TIME SERIES (BTC only)
###########################

keeps <- readRDS(paste0("keeps_long/keeps_",SVMALD ,"_","BTC", ".rds"))
names <- lapply(keeps[c(4,6:17)], domean, total = 10) %>%unlist %>%names

keeps_summary <- array(dim = c(16, length(names)))
keeps_v1 <- array(dim = c(16,dim(keeps$v)[2]))
keeps_v2 <- array(dim = c(16,dim(keeps$v)[2]))
model <- rep(NA, 16)
colnames(keeps_summary) <- names

j = 0
for (m in c("SVMALD", "SVMVN", "SVLD", "SVIND")){
    j = j + 1
    keeps <- readRDS(paste0("keeps_long/keeps_",m ,"_","BTC", ".rds"))
    keeps_summary[j,] <- lapply(keeps[c(4,6:17)], domean, total = 20000) %>%unlist
    keeps_v1[j,] <- apply(keeps$v[1:R,,1], 2:3, function(x){(mean(x))}) #alternative sv
    keeps_v2[j,] <- apply(keeps$v[1:R,,2], 2:3, function(x){(mean(x))}) #sp sv
    model[j] <-m
}


keeps_summary <- keeps_summary %>%as.data.frame()%>%
  mutate_all(round, digits = 2) %>%
  mutate(lambda = paste0("(",lambda1,", ", lambda2,", ", lambda3,", ", lambda4,")"))%>%
  dplyr::select(-c(lambda1, lambda2, lambda3, lambda4))
keeps_summary$model <- model

#TABLE XXX POSTERIOR MEANS OF PARAMETERS------
keeps_summary[,c(ncol(keeps_summary), ncol(keeps_summary)-1, 1:(ncol(keeps_summary)-2))] %>%
  xtable() %>%
  print()

#FIGURE XXX STOCHASTIC VOLATILITY--------
keeps_v1_long <- keeps_v1 %>%
  mutate(model = model) %>%
  melt(id.vars = c(model)) %>%
  mutate(Date = Date,
         series = "BTC")

keeps_v2_long <- keeps_v2 %>%
  mutate(model = model) %>%
  melt(id.vars = c(model)) %>%
  mutate(Date = Date,
         series = "S&P")
  
  
  rbind(keeps_v1_long,keeps_v2_long) %>%
  ggplot() +
  geom_line(aes(x = Var2, y = value, colour = model)) +
    facet_grid(series~model, scales = "free_y")

