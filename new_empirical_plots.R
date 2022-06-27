rm(list = ls()) 
library(ald)
library(ggplot2)
library(grid)
library(gridExtra)
library(dplyr)
library(quantmod)
library(tmvtnorm)

get_qq <- function(keeps, data){## Function to make the QQ Plots
  #browser()
  # delta <- rep(0,T)
  # V <- array(0, dim = c(T+1, 2))
  # J <- y <- x <- array(0, dim = c(T, 2))
  # V[1,] <- apply(keeps$v[,1,], 2, mean)
  y <- x <- array(0, dim = c(T, 2))
  
  sim <- 0
  for (t in 1:T){
    print(t)
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
  return(ks.test(QQdat[,1],QQdat[,3])$p.value)
  
  # p1 <- ggplot() +
  #   geom_point(aes(x=quantile(QQdat$V1,seq(0.01,0.99,0.01)),y=quantile(QQdat$V3,seq(0.01,0.99,0.01)))) +
  #   geom_abline(slope=1,intercept=0) +
  #   #xlim(c(-15,15)) + ylim(c(-15,15)) +
  #   xlab("Actual Quantiles") + ylab("Simulated Quantiles") + theme_bw() + ggtitle(paste0(substr(m,1,2),"-",substr(m,3,6)))
  # #p1
  # 
  # return(p1 + theme(plot.title = element_text(hjust = 0.5,size = 20)))
}

BTC <- read.csv("BTC-USD.csv")
BTC <- BTC %>% dplyr::select(Date,Close) %>% rename(`BTC-USD.Close`=Close)
SP <- read.csv("SPY.csv")
SP <- SP %>% dplyr::select(Date,Close) %>% rename(`GSPC.Close`=Close)
S <- merge(BTC,SP)
S$Date <- as.Date(S$Date)
T <- nrow(S) - 1

QQdat <- cbind(100*log(S$`BTC-USD.Close`[-1]) - 100*log(S$`BTC-USD.Close`[-(T+1)]),
               100*log(S$`GSPC.Close`[-1]) - 100*log(S$`GSPC.Close`[-(T+1)]))

data.frame(Date=S$Date[-1],
                BTC=100*log(S$`BTC-USD.Close`[-1]) - 100*log(S$`BTC-USD.Close`[-(T+1)]),
                SP=100*log(S$`GSPC.Close`[-1]) - 100*log(S$`GSPC.Close`[-(T+1)])) %>%
  gather("Asset","Value",-Date) %>%
  mutate(Asset = replace(Asset,Asset=="SP","S&P 500")) %>%
  ggplot() + 
  geom_line(aes(x=Date,y=Value)) +
  facet_grid(Asset~.,scales="free_y") +
  ylab("% Return") +
  theme_bw()
ggsave("data_plot.pdf", height = 10, width = 14)

keepsBTCSP_MALD <- readRDS("keeps/keepsBTCSP_MALD.rds")
get_qq(keepsBTCSP_MALD,QQdat)
MALD_v1 <- sqrt(252*apply(keepsBTCSP_MALD$v[1:20000,,1], 2, mean))
MALD_v2 <- sqrt(252*apply(keepsBTCSP_MALD$v[1:20000,,2], 2, mean))
MALD_J1 <- apply(keepsBTCSP_MALD$J[1:20000,,1], 2, mean)
MALD_J2 <- apply(keepsBTCSP_MALD$J[1:20000,,2], 2, mean)
rm(keepsBTCSP_MALD)
keepsBTCSP_IND <- readRDS("keeps/keepsBTCSP_IND.rds")
get_qq(keepsBTCSP_IND,QQdat)
IND_v1 <- sqrt(252*apply(keepsBTCSP_IND$v[1:20000,,1], 2, mean))
IND_v2 <- sqrt(252*apply(keepsBTCSP_IND$v[1:20000,,2], 2, mean))
IND_J1 <- apply(keepsBTCSP_IND$J[1:20000,,1], 2, mean)
IND_J2 <- apply(keepsBTCSP_IND$J[1:20000,,2], 2, mean)
rm(keepsBTCSP_IND)
keepsBTCSP_LD <- readRDS("keeps/keepsBTCSP_LD.rds")
get_qq(keepsBTCSP_LD,QQdat)
LD_v1 <- sqrt(252*apply(keepsBTCSP_LD$v[1:20000,,1], 2, mean))
LD_v2 <- sqrt(252*apply(keepsBTCSP_LD$v[1:20000,,2], 2, mean))
LD_J1 <- apply(keepsBTCSP_LD$J[1:20000,,1], 2, mean)
LD_J2 <- apply(keepsBTCSP_LD$J[1:20000,,2], 2, mean)
rm(keepsBTCSP_LD)
keepsBTCSP_MVN <- readRDS("keeps/keepsBTCSP_MVN.rds")
get_qq(keepsBTCSP_MVN,QQdat)
MVN_v1 <- sqrt(252*apply(keepsBTCSP_MVN$v[1:20000,,1], 2, mean))
MVN_v2 <- sqrt(252*apply(keepsBTCSP_MVN$v[1:20000,,2], 2, mean))
MVN_J1 <- apply(keepsBTCSP_MVN$J[1:20000,,1], 2, mean)
MVN_J2 <- apply(keepsBTCSP_MVN$J[1:20000,,2], 2, mean)
rm(keepsBTCSP_MVN)

dt <- data.frame(Date=S$Date,
           S=S$`BTC-USD.Close`,
           `MALD`=MALD_v1,
           `IND`=IND_v1,
           `LD`=LD_v1,
           `MVN`=MVN_v1) %>%
  gather("Model","V",-c(Date,S))
dt$Model <- factor(dt$Model,
                   levels=c("IND","LD","MVN","MALD"),
                   labels=c("SV-IND","SV-LD","SV-MVN","SV-MALD"))
p <- dt %>%
  gather("Var","Value",-c(Date,Model)) %>%
  ggplot() +
  geom_line(aes(x = Date, y = Value, linetype = Model)) +
  facet_grid(Var~., scales = "free_y") +
  theme_bw() +
  scale_colour_grey() +   
  labs(x = "Date")+ theme(text = element_text(size = 20))
ggsave(paste0("BTC_Vol.pdf"),p, width = 10, height = 14)

dt <- data.frame(Date=S$Date,
                 S=100*S$`GSPC.Close`,
                 `MALD`=MALD_v2,
                 `IND`=IND_v2,
                 `LD`=LD_v2,
                 `MVN`=MVN_v2) %>%
  gather("Model","V",-c(Date,S))
dt$Model <- factor(dt$Model,
                   levels=c("IND","LD","MVN","MALD"),
                   labels=c("SV-IND","SV-LD","SV-MVN","SV-MALD"))
p <- dt %>%
  gather("Var","Value",-c(Date,Model)) %>%
  ggplot() +
  geom_line(aes(x = Date, y = Value, linetype = Model)) +
  facet_grid(Var~., scales = "free_y") +
  theme_bw() +
  scale_colour_grey() +   
  labs(x = "Date")+ theme(text = element_text(size = 20))
ggsave(paste0("SP_Vol.pdf"),p, width = 10, height = 14)

BTC_J <- data.frame(Date=S$Date[-1],
                 MALD = MALD_J1,
                 IND = IND_J1,
                 LD = LD_J1,
                 MVN = MVN_J1,
                 series = "BTC") %>%
  gather("Model","Value",-c(Date,series))

SP_J <- data.frame(Date=S$Date[-1],
                  MALD = MALD_J2,
                  IND = IND_J2,
                  LD = LD_J2,
                  MVN = MVN_J2,
                  series = "SP") %>%
  gather("Model","Value",-c(Date,series))

p <- rbind(BTC_J,SP_J) %>%
  mutate(Model = factor(Model,
                        levels=c("IND","LD","MVN","MALD"),
                        labels=c("SV-IND","SV-LD","SV-MVN","SV-MALD"))) %>%
  ggplot() +
  geom_line(aes(x = Date, y = Value, linetype = Model)) +
  facet_grid(series~Model, scales = "free_y") +
  theme_bw() +
  scale_colour_grey() +
  labs(x = "Date", y = "Jump size")+ theme(text = element_text(size = 20))
ggsave("jump_sizes.pdf",p, height = 10, width = 14) 


