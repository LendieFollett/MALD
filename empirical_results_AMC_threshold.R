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
library(latex2exp)
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

get_qq_short <- function(keeps,data,i,threshold){## Function to make the QQ Plots
  #browser()
  # delta <- rep(0,T)
  # V <- array(0, dim = c(T+1, 2))
  # J <- y <- x <- array(0, dim = c(T, 2))
  # V[1,] <- apply(keeps$v[,1,], 2, mean)
  y <- x <- array(0, dim = c(T, 2))
  
  sim <- 0
  for (t in 1:T){
    # print(t)
    # set.seed(t + 4922035+ sim)
    # delta <- sample(c(0:3),prob=apply(keeps$lambda, 2, mean), 1)
    # 
    # set.seed(15866245 + t + sim)
    # if (m == "MVN"){
    #   B <- 1
    # } else {
    #   B <- rexp(1)
    # }
    # xi_y1 <- mean(keeps$xi_y1w)*B + sqrt(B)*rnorm(1,0,mean(keeps$xi_y1eta)) #SHOULD THIS BE SQRT(ETA)? No
    # if (m == "MVN"){
    #   B <- 1
    # } else {
    #   B <- rexp(1)
    # }
    # xi_y2 <- mean(keeps$xi_y2w)*B + sqrt(B)*rnorm(1,0,mean(keeps$xi_y2eta)) #SHOULD THIS BE SQRT(ETA)? No
    # if (m == "IND"){
    #   xi_c = cbind(xi_y1,xi_y2)
    # } else {
    #   if (m == "MVN"){
    #     B <- 1
    #   } else {
    #     B <- rexp(1)
    #   }
    #   Sigma <- matrix(c(mean(keeps$sigma_c[,1])^2,
    #                     mean(keeps$rhoc)*mean(keeps$sigma_c[,1])*mean(keeps$sigma_c[,2]),
    #                     mean(keeps$rhoc)*mean(keeps$sigma_c[,1])*mean(keeps$sigma_c[,2]),
    #                     mean(keeps$sigma_c[,2]^2)),
    #                   nrow=2)
    #   xi_c <- apply(keeps$xi_cw, 2, mean)*B+
    #     sqrt(B)*rtmvnorm(n = 1, mean = c(0,0), sigma = Sigma)
    # }
    # 
    # J[t,] = xi_c*(delta==2) +cbind(xi_y1,0)*(delta==0) + cbind(0,xi_y2)*(delta==1)
    J = apply(keeps$J[,t,], 2, mean)
    V = apply(keeps$v[,t,], 2, mean)
    
    # Sigma <- matrix(c(V[t,1],
    #                   mean(keeps$rho[,1])*sqrt(prod(V[t,])),
    #                   mean(keeps$rho[,3])*mean(keeps$sigma_v[,1])*V[t,1],0,
    #                   mean(keeps$rho[,1])*sqrt(prod(V[t,])),V[t,2],0,mean(keeps$rho[,4])*mean(keeps$sigma_v[,2])*V[t,2],
    #                   mean(keeps$rho[,3])*mean(keeps$sigma_v[,1])*V[t,1],0,mean(keeps$sigma_v[,1])^2*V[t,1],mean(keeps$rho[,2])*prod(apply(keeps$sigma_v, 2, mean))*sqrt(prod(V[t,])),
    #                   0,mean(keeps$rho[,4])*mean(keeps$sigma_v[,2])*V[t,2],mean(keeps$rho[,2])*prod(apply(keeps$sigma_v, 2, mean))*sqrt(prod(V[t,])),mean(keeps$sigma_v[,2])^2*V[t,2]),nrow=4)
    
    Sigma <- matrix(c(V[1],
                      mean(keeps$rho[,1])*sqrt(prod(V)),
                      mean(keeps$rho[,3])*mean(keeps$sigma_v[,1])*V[1],0,
                      mean(keeps$rho[,1])*sqrt(prod(V)),V[2],0,mean(keeps$rho[,4])*mean(keeps$sigma_v[,2])*V[2],
                      mean(keeps$rho[,3])*mean(keeps$sigma_v[,1])*V[1],0,mean(keeps$sigma_v[,1])^2*V[1],mean(keeps$rho[,2])*prod(apply(keeps$sigma_v, 2, mean))*sqrt(prod(V)),
                      0,mean(keeps$rho[,4])*mean(keeps$sigma_v[,2])*V[2],mean(keeps$rho[,2])*prod(apply(keeps$sigma_v, 2, mean))*sqrt(prod(V)),mean(keeps$sigma_v[,2])^2*V[2]),nrow=4)
    set.seed(463468+t)
    temp <- rtmvnorm(n = 1,
                     mean = c(apply(keeps$mu, 2 ,mean) + J + Sigma[1:2,3:4] %*% solve(Sigma[3:4,3:4]) %*% 
                                (apply(keeps$v[,t+1,],2,mean) - (apply(keeps$theta,2, mean) + apply(keeps$phi, 2, mean)*(V - apply(keeps$theta, 2, mean))))),
                     sigma = Sigma[1:2,1:2] - Sigma[1:2,3:4] %*% solve(Sigma[3:4,3:4]) %*% Sigma[3:4,1:2])
    
    y[t,] <- temp
    if( t+1 <= T){ x[t+1] <- 0 }
  }
  
  QQdat = cbind(data,y)
  names(QQdat) = c("V1","V2","V3","V4")
  # return(ks.test(QQdat$V1,QQdat$V3)$p.value)
  
  p1 <- ggplot() +
    geom_point(aes(x=quantile(QQdat$V1,seq(0.01,0.99,0.01)),y=quantile(QQdat$V3,seq(0.01,0.99,0.01)))) +
    geom_abline(slope=1,intercept=0) +
    #xlim(c(-15,15)) + ylim(c(-15,15)) +
    xlab("Actual Quantiles") + ylab("Simulated Quantiles") + theme_bw() + ggtitle(paste0("Threshold = ",round(threshold,2)))
  #p1

  return(p1 + theme(plot.title = element_text(hjust = 0.5,size = 20)))
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
getSymbols("AMC",from = "2021-05-15",to = "2021-12-31")
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

for (i in c("AMC")){
  m <- "SVMALD"
  keeps <- readRDS(paste0("keeps/keeps_AMC/keeps_",m ,"_",i,".rds"))
  summary <- NULL
  thresholds <- sqrt(252*apply(keeps$v[1:20000,,1], 2, mean)) %>% quantile(probs=seq(0.15,0.25,0.01))
  if (i == "BTC"){
    y <- 100*(log(S[-1,c("BTC-USD.Close","GSPC.Close")]) - 
                     log(S[-nrow(S),c("BTC-USD.Close","GSPC.Close")]))
  } else if (i == "DOGE"){
    y <- 100*(log(S[-1,c("DOGE-USD.Close","GSPC.Close")]) - 
                     log(S[-nrow(S),c("DOGE-USD.Close","GSPC.Close")]))
  } else if (i == "AMC"){
    y <- 100*(log(S[-1,c("AMC.Close","GSPC.Close")]) - 
                     log(S[-nrow(S),c("AMC.Close","GSPC.Close")]))
  } else if (i == "GME") {
    y <- 100*(log(S[-1,c("GME.Close","GSPC.Close")]) - 
                     log(S[-nrow(S),c("GME.Close","GSPC.Close")]))
  } else {
    y <- 100*(log(S[-1,c("MRNA.Close","GSPC.Close")]) - 
                     log(S[-nrow(S),c("MRNA.Close","GSPC.Close")]))
  }
  j = 0
  prob_rho_1_mn <- rep(NA,11)
  prob_rho_1_pl_1_mn <- rep(NA,11)
  ks_p_value <- rep(NA,11)
  data <- rep(i, 1)
  keeps_v1_long <- NULL
  for (thresh in thresholds){
    j = j + 1
    keeps <- readRDS(paste0("keeps/keeps_AMC/keeps_",m ,"_",i, "_",thresh,".rds"))
    keeps_v1 <- apply(keeps$v[1:20000,,1], 2, mean) #alternative sv
    prob_rho_1_mn[j] <- length(which(keeps$rho[,4] > 0))/20000
    prob_rho_1_pl_1_mn[j] <- length(which(keeps$rho[,3] > keeps$rho[,4]))/20000
    # keeps_v1_long <- rbind(keeps_v1_long,
    #   keeps_v1 %>% as.data.frame()%>%
    #   mutate(series = data,
    #          threshold = as.character(j)) %>%
    #   melt(id.vars = c("series","threshold"))%>%
    #   mutate(Date = Date)) #CHECK THIS STRUCTURE
    # ks_p_value[j] <- get_qq_short(keeps,y,i,m)
    if ((j == 1) | (j == 11)){
      plot <- get_qq_short(keeps,y,i,thresh)
      assign(paste0("plot",j),plot)
    }
  }
  # assign(paste0(i,"_p_values"),ks_p_value)
  p1 <- data.frame(Threshold=thresholds,Probability=prob_rho_1_mn) %>%
    ggplot(aes(x=Threshold,y=Probability)) +
    geom_point() + 
    xlim(floor(thresholds[1]),ceiling(thresholds[11])) +
    ylab(TeX("$P(\\rho_{1-} > 0 | {\\textbf{y}})$")) + 
    stat_smooth(method="lm",se=FALSE,colour="black",fullrange=TRUE) +
    theme_bw()
  p2 <- data.frame(Threshold=thresholds,Probability=prob_rho_1_pl_1_mn) %>%
    ggplot(aes(x=Threshold,y=Probability)) +
    geom_point() + 
    xlim(floor(thresholds[1]),ceiling(thresholds[11])) +
    ylab(TeX("$P(\\rho_{1+} > \\rho_{1-} | {\\textbf{y}})$")) + 
    stat_smooth(method="lm",se=FALSE,colour="black",fullrange=TRUE) +
    theme_bw()
  p <- grid.arrange(p1,p2, nrow = 1)
  ggsave(paste0("Probability_Plots_",i, ".pdf"),
         p, width = 10, height = 7)
  p <- grid.arrange(plot1,plot11,nrow = 2)
  ggsave(paste0("QQ_Plots_",i,"_threshold", ".pdf"),
         p, width = 12, height = 12)
  # mod <- summary(lm(prob_rho_1_mn ~ thresholds))
  # print(ifelse(mod$coef[2,1] < 0,
  #              1 - mod$coef[2,4] / 2,
  #              mod$coef[2,4] / 2))
  # print(mod$coef[2,1] + 1.96*c(-1,1)*mod$coef[2,2])
  # mod <- summary(lm(prob_rho_1_pl_1_mn ~ thresholds))
  # print(ifelse(mod$coef[2,1] > 0,
  #              1 - mod$coef[2,4] / 2,
  #              mod$coef[2,4] / 2))
  # print(mod$coef[2,1] + 1.96*c(-1,1)*mod$coef[2,2])
}
