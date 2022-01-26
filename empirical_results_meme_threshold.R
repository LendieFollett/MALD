rm(list = ls()) 
library(ald)
library(ggplot2)
library(grid)
library(gridExtra)
library(dplyr)
library(truncnorm)
library(mvtnorm)
library(tmvtnorm)
library(Rcpp)
library(MCMCpack)
library(quantmod)
library(RcppTN)
library(xtable)
library(tibble)#rownames_to_column()
library(tidyr)
library(reshape2)

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

############################################################
#----SHORT TIME SERIES ONLY---------------------------------
#ALL CRYPTO, MEME STOCKS
############################################################
#load data
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
Date <- S$Date

prob_rho_1 <- NULL
prob_rho_2 <- NULL
prob_rho_12 <- NULL

for (i in c("BTC", "DOGE", "AMC", "GME","MRNA")){
  summary <- NULL
  if (i == "BTC"){
    y <- 100*(log(S[-1,c("BTC-USD.Close","GSPC.Close")]) - 
                     log(S[-nrow(S),c("BTC-USD.Close","GSPC.Close")]))
    thresholds <- c(60:65)
  } else if (i == "DOGE"){
    y <- 100*(log(S[-1,c("DOGE-USD.Close","GSPC.Close")]) - 
                     log(S[-nrow(S),c("DOGE-USD.Close","GSPC.Close")]))
    thresholds <- c(60:65)
  } else if (i == "AMC"){
    y <- 100*(log(S[-1,c("AMC.Close","GSPC.Close")]) - 
                     log(S[-nrow(S),c("AMC.Close","GSPC.Close")]))
    thresholds <- c(75:80)
  } else if (i == "GME") {
    y <- 100*(log(S[-1,c("GME.Close","GSPC.Close")]) - 
                     log(S[-nrow(S),c("GME.Close","GSPC.Close")]))
    thresholds <- c(55:60)
  } else {
    y <- 100*(log(S[-1,c("MRNA.Close","GSPC.Close")]) - 
                     log(S[-nrow(S),c("MRNA.Close","GSPC.Close")]))
    thresholds <- c(60:65)
  }
  j = 0
  thr <- rep(NA, 6)
  data <- rep(i, 1)
  keeps_v1_long <- NULL
  for (m in c("SVMALD")){
    for (thresh in thresholds){
      j = j + 1
      keeps <- readRDS(paste0("keeps/keeps_",m ,"_",i, "_",thresh,".rds"))
      keeps_v1 <- apply(keeps$v[1:20000,,1], 2, mean) #alternative sv
      keeps_v2 <- apply(keeps$v[1:20000,,2], 2, mean) #sp sv
      keeps_j1 <- apply(keeps$J[1:20000,,1], 2, mean) #alternative sv
      keeps_j2 <- apply(keeps$J[1:20000,,2], 2, mean) #sp sv
      thr[j] <- thresh
      summary <- c(summary,
                   paste0(round(median(keeps$rho[,3]),2),
                   "(",round(quantile(keeps$rho[,3],0.025),3),
                   ",",round(quantile(keeps$rho[,3],0.975),3),
                   ")"),
                   paste0(round(median(keeps$rho[,4]),2),
                   "(",round(quantile(keeps$rho[,4],0.025),3),
                   ",",round(quantile(keeps$rho[,4],0.975),3),
                   ")"),
      round(length(which(keeps$rho[,3] > 0))/20000,3),
      round(length(which(keeps$rho[,4] > 0))/20000,3),
      round(length(which(keeps$rho[,3] > keeps$rho[,4]))/20000,3),
      round(ks.test((y[,1]-mean(keeps$mu[,1])-keeps_j1)/sqrt(keeps_v1[-(T+1)]),"pnorm")$p.value,3))
      keeps_v1_long <- rbind(keeps_v1_long,
        keeps_v1 %>% as.data.frame()%>%
        mutate(series = data,
               threshold = as.character(j)) %>%
        melt(id.vars = c("series","threshold"))%>%
        mutate(Date = Date)) #CHECK THIS STRUCTURE
    }
  }
  summary %>% matrix(nrow=6,byrow=T) %>% xtable() %>% print()
  # 
  # keeps_summary <- keeps_summary %>% as.data.frame()%>%
  #   #mutate_all(round, digits = 2) %>%
  #   #mutate(lambda = paste0("(",lambda1,", ", lambda2,", ", lambda3,", ", lambda4,")"))%>%
  #   #dplyr::select(-c(lambda1, lambda2, lambda3, lambda4)) %>%
  #   mutate(model = model) %>%t()
  # 
  # #TABLE XXX POSTERIOR MEANS OF PARAMETERS------
  # keeps_summary%>%
  #   xtable() %>%
  #   print()
  # # 
  #FIGURE XXX STOCHASTIC VOLATILITY--------
  #ALTERNATIVE CURRENCY
  assign(paste0("keeps_",i,"_long"),keeps_v1_long)
}

V <- rbind(keeps_AMC_long,keeps_BTC_long,keeps_DOGE_long,keeps_GME_long,keeps_MRNA_long)
p <- ggplot(V) +
  geom_line(aes(x = Date, y = sqrt(252*value), linetype = threshold)) +
  facet_grid(series~., scales = "free_y") +
  theme_bw() +
  scale_colour_grey() +
  labs(x = "Date", y = "Volatility")+ 
  theme(text = element_text(size = 20))
ggsave(paste0("Volatility_", "all_threshold", ".pdf"),p, width = 10, height = 14)
