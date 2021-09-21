rm(list=ls())
library(tidyverse)
library(quantmod)
##### RHO PLOTS FOR LONG TIME SERIES #####
getSymbols("BTC-USD",from = "2014-09-15",to = "2020-09-30")
BTC <- as.data.frame(`BTC-USD`)
BTC$Date <- seq(as.Date("2014-09-17"),as.Date("2020-09-30"),by="days")
BTC$`BTC-USD.Close`[BTC$Date=="2020-04-17"] <- 7096.18
getSymbols("^GSPC",from = "2014-09-15",to = "2020-09-30")
SP500 <- as.data.frame(`GSPC`)
SP500$Date <- as.Date(rownames(SP500))
S <- merge(BTC,SP500)
T <- nrow(S) - 1
Date <- S$Date
tmp <- data.frame(Date=S$Date[-1])
for (m in c("IND","LD","MVN","MALD")){
  keeps <- readRDS(paste0("keeps_long/keepsBTCSP_",m,".rds"))
  mu <- keeps$mu[,1]
  theta <- keeps$theta[,1]
  kappa <- 1 - keeps$phi[,1]
  sigmav <- keeps$sigma_v[,1]
  J <- keeps$J[,,1]
  V <- keeps$v[,,1]
  eps_y <- NULL
  eps_v <- NULL
  for (k in 1:20000){
    print(k)
    eps_y <- cbind(eps_y,(y - mu[k] - J[k,])/sqrt(V[k,-(T+1)]))
    eps_v <- cbind(eps_v,(V[k,-1] - V[k,-(T+1)] - kappa[k]*(theta[k] - V[k,-(T+1)]))/(sigmav[k]*sqrt(V[k,-(T+1)])))
  }
  rho <- (apply(eps_y*eps_v,1,mean) - apply(eps_y,1,mean)*apply(eps_v,1,mean))/sqrt(apply(eps_y,1,var)*apply(eps_v,1,var))
  tmp$tm <- rho
  names(tmp)[names(tmp)=="tm"] = m
}
p <- tmp %>% gather("Model","Rho",-Date) %>%
  ggplot() +
  geom_line(aes(x=Date,y=Rho,colour=Model)) + theme_bw() + theme(text = element_text(size = 20))
ggsave("rho_t.pdf",p, height = 10, width = 14) 


##### RHO PLOTS FOR SHORT TIME SERIES #####
getSymbols("^GSPC",from = "2020-10-01",to = "2021-06-30")
SP500 <- as.data.frame(`GSPC`)
SP500$Date <- as.Date(rownames(SP500))
BTC <- read.csv("BTC_USD_2020-07-01_2021-06-30-CoinDesk.csv")
BTC <- BTC %>% dplyr::select(!Currency)
names(BTC) <- c("Date","BTC-USD.Close","BTC-USD.Open","BTC-USD.High","BTC-USD.Low")
BTC$Date <- as.Date(BTC$Date)
DOGE <- read.csv("DOGE_USD_2020-07-01_2021-06-30-CoinDesk.csv")
DOGE <- DOGE %>% dplyr::select(!Currency)
DOGE$Date <- as.Date(DOGE$Date)
names(DOGE) <- c("Date","DOGE-USD.Close","DOGE-USD.Open","DOGE-USD.High","DOGE-USD.Low")
getSymbols("AMC",from = "2020-10-01",to = "2021-06-30")
AMC <- as.data.frame(AMC)
AMC$Date <- as.Date(rownames(AMC))
getSymbols("GME",from = "2020-10-01",to = "2021-06-30")
GME <- as.data.frame(GME)
GME$Date <- as.Date(rownames(GME))
getSymbols("MRNA",from = "2020-10-01",to = "2021-06-30")
MRNA <- as.data.frame(MRNA)
MRNA$Date <- as.Date(rownames(MRNA))
getSymbols("DIS",from = "2020-10-01",to = "2021-06-30")
DIS <- as.data.frame(DIS)
DIS$Date <- as.Date(rownames(DIS))
getSymbols("BBY",from = "2020-10-01",to = "2021-06-30")
BBY <- as.data.frame(BBY)
BBY$Date <- as.Date(rownames(BBY))
getSymbols("BMY",from = "2020-10-01",to = "2021-06-30")
BMY <- as.data.frame(BMY)
BMY$Date <- as.Date(rownames(BMY))

S <- BTC %>% merge(DOGE) %>% merge(AMC) %>% merge(GME) %>% merge(MRNA) %>% merge(DIS) %>%
  merge(BBY) %>% merge(BMY)
T <- nrow(S) - 1

for (i in c("BTC","DOGE","AMC","GME","MRNA","DIS","BBY","BMY")){
  if (i == "BTC"){
    y = log(S$`BTC-USD.Close`)[-1] - log(S$`BTC-USD.Close`)[-(T+1)]
  } else if (i == "DOGE"){
    y = log(S$`DOGE-USD.Close`)[-1] - log(S$`DOGE-USD.Close`)[-(T+1)]
  }else if (i == "AMC"){
    y = log(S$AMC.Close)[-1] - log(S$AMC.Close)[-(T+1)]
  } else if (i == "GME"){
    y = log(S$GME.Close)[-1] - log(S$GME.Close)[-(T+1)]
  } else if (i == "MRNA"){
    y = log(S$MRNA.Close)[-1] - log(S$MRNA.Close)[-(T+1)]
  } else if (i == "DIS"){
    y = log(S$DIS.Close)[-1] - log(S$DIS.Close)[-(T+1)]
  } else if (i == "BBY"){
    y = log(S$BBY.Close)[-1] - log(S$BBY.Close)[-(T+1)]
  } else {
    y = log(S$BMY.Close)[-1] - log(S$BMY.Close)[-(T+1)]
  }
  tmp <- data.frame(Date=S$Date[-1])
  for (m in c("SVIND","SVLD","SVMVN","SVMALD")){
    keeps <- readRDS(paste0("keeps/keeps_",m,"_",i,".rds"))
    mu <- keeps$mu[,1]
    theta <- keeps$theta[,1]
    kappa <- 1 - keeps$phi[,1]
    sigmav <- keeps$sigma_v[,1]
    J <- keeps$J[,,1]
    V <- keeps$v[,,1]
    eps_y <- NULL
    eps_v <- NULL
    for (k in 1:20000){
      print(k)
      eps_y <- cbind(eps_y,(y - mu[k] - J[k,])/sqrt(V[k,-(T+1)]))
      eps_v <- cbind(eps_v,(V[k,-1] - V[k,-(T+1)] - kappa[k]*(theta[k] - V[k,-(T+1)]))/(sigmav[k]*sqrt(V[k,-(T+1)])))
    }
    rho <- (apply(eps_y*eps_v,1,mean) - apply(eps_y,1,mean)*apply(eps_v,1,mean))/sqrt(apply(eps_y,1,var)*apply(eps_v,1,var))
    tmp$tm <- rho
    names(tmp)[names(tmp)=="tm"] = m
  }
  p <- tmp %>% gather("Model","Rho",-Date) %>%
    ggplot() +
    geom_line(aes(x=Date,y=Rho,colour=Model)) + theme_bw() + theme(text = element_text(size = 20))
  ggsave(paste0(i,"_rho_t.pdf"),p, height = 10, width = 14) 
}

