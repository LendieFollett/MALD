library(tidyverse)
library(quantmod)

getSymbols("^GSPC",from = "2020-10-01",to = "2021-06-30")
SP500 <- as.data.frame(`GSPC`)
SP500$Date <- as.Date(rownames(SP500))
getSymbols("BTC-USD",from = "2020-10-01",to = "2021-06-30")
BTC <- as.data.frame(BTC)
BTC$Date <- as.Date(rownames(BTC))
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

S <- BTC %>% merge(AMC) %>% merge(GME) %>% merge(MRNA) %>% merge(DIS) %>%
  merge(BBY) %>% merge(BMY)
T <- nrow(S) - 1

for (i in c("BTC","AMC","GME","MRNA","DIS","BBY","BMY")){
  if (i == "BTC"){
    y = log(S$`BTC-USD.close`)[-1] - log(S$`BTC-USD.close`)[-(T+1)]
  } else if (i == "AMC"){
    y = log(S$AMC.close)[-1] - log(S$AMC.close)[-(T+1)]
  } else if (i == "GME"){
    y = log(S$GME.close)[-1] - log(S$GME.close)[-(T+1)]
  } else if (i == "MRNA"){
    y = log(S$MRNA.close)[-1] - log(S$MRNA.close)[-(T+1)]
  } else if (i == "DIS"){
    y = log(S$DIS.close)[-1] - log(S$DIS.close)[-(T+1)]
  } else if (i == "BBY"){
    y = log(S$BBY.close)[-1] - log(S$BBY.close)[-(T+1)]
  } else {
    y = log(S$BMY.close)[-1] - log(S$BMY.close)[-(T+1)]
  }
  for (m in c("SVIND","SVLD","SVMVN","SVMALD")){
    keeps <- readRDS(paste0("keeps/keeps_",i,"_",m,".rds"))
    mu <- apply(keeps$mu,2,mean)
    theta <- apply(keeps$theta,2,mean)
    kappa <- 1 - apply(keeps$phi,2,mean)
    sigmav <- apply(keeps$sigma_v,2,mean)
    J <- apply(keeps$J,c(1,3),mean)
    V <- apply(keeps$V,c(1,3),mean)
    eps_y <- (y - mu - J)/sqrt(V[-(T+1)])
    eps_v <- (V[-1] - V[-(T+1)] - kappa*(theta - V[-(T+1)]))/sqrt(V[-(T+1)])
    rho <- eps_y*eps_v
    ggplot() + geom_point(aes(x=1:T,y=rho)) + xlab("Time") + ylab("rho")
  }
}

