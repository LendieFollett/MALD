library(ggplot2)
library(grid)
library(gridExtra)
library(MASS)
library(tmvtnorm)
library(ald)
library(ggplot2)
library(dplyr)
library(truncnorm)
library(mvtnorm)
library(tmvtnorm)
library(Rcpp)
library(MCMCpack)
library(quantmod)


keepsBTCSP <- readRDS("/Users/000766412/Box Sync/ALD_Codes/keepsBTCSP.rds")
keepsBTC <- readRDS("/Users/000766412/Box Sync/ALD_Codes/keepsBTC.rds")
keepsSP <- readRDS("/Users/000766412/Box Sync/ALD_Codes/keepsSP.rds")

## Generate SP, BTC Data for 1-d, 2-d model
getSymbols("BTC-USD",from = "2014-09-15",to = "2020-09-30")
BTC <- as.data.frame(`BTC-USD`)
BTC$Date <- seq(as.Date("2014-09-17"),as.Date("2020-09-30"),by="days")
BTC$`BTC-USD.Close`[BTC$Date=="2020-04-17"] <- 7096.18
getSymbols("^GSPC",from = "2014-09-15",to = "2020-09-30")
SP500 <- as.data.frame(`GSPC`)
SP500$Date <- as.Date(rownames(SP500))
Rsequence <- seq(1,nrow(keepsBTC$v), by = 100)
S <- merge(BTC,SP500)
T <- nrow(S) - 1
S <- merge(BTC,SP500)
y <- as.matrix(100*(log(S[-1,c("BTC-USD.Close","GSPC.Close")]) - log(S[-nrow(S),c("BTC-USD.Close","GSPC.Close")])))
x <- array(0,dim=dim(y))
r2 <- 0
y <-array(0, dim = c(T, 2, length(Rsequence)))

#2d model
for(r in Rsequence){ #loop over posterior draws -->posterior predictive distributions
  print(r)
  r2 <- r2 + 1
  delta <- rep(0,T)
  xiy1 <- xiy2 <- xiy1s <- xiy2s <- xics <- rep(0,T)
  V <- array(0,dim=c(T+1,2))
  xic <- J <-  x <- array(0,dim=c(T,2))
  V[1,] <- keepsBTCSP$v[r,1,]

  sim <- 0 
  for ( i in 1:T){
  
  # Sigma_c = matrix(c((keepsBTCSP$sigma_c[r,1])^2,
  #                    (keepsBTCSP$rhoc[r])*prod(keepsBTCSP$sigma_c[r,]),
  #                    (keepsBTCSP$rhoc[r])*prod(keepsBTCSP$sigma_c[r,]),
  #                    (keepsBTCSP$sigma_c[r,2])^2),nrow=2)
  # whichone <- sample(c(0:3),prob=apply(keepsBTCSP$lambda,2,mean), 1)
  # delta[i] <- whichone
  # xics[i] <- rexp(1)
  # X<- mvrnorm(n = 1, mu = c(0,0), Sigma = Sigma_c)
  # xic[i,] <- (keepsBTCSP$xi_cw[r,]*xics[i]  + sqrt(xics[i])*X)
  # 
  # 
  # xiy1s[i] <- rexp(1)
  # xiy1[i] <- keepsBTCSP$xi_y1w[r]*xiy1s[i]+ 
  #               sqrt(xiy1s[i])*keepsBTCSP$xi_y1eta[r]*rnorm(1,0,1)
  # xiy2s[i] <- rexp(1)
  # xiy2[i] <- keepsBTCSP$xi_y2w[r]*xiy2s[i]+ 
  #               sqrt(xiy2s[i])*keepsBTCSP$xi_y2eta[r]*rnorm(1,0,1)
  # 
  # J[i,] = (delta[i] == 2)*xic[i,] + (delta[i] == 0)*c(xiy1[i],0) + (delta[i] == 1)*c(0,xiy2[i])
  # 
  J[i,] = keepsBTCSP$J[r,i,]  
    
  Sigma <- matrix(c(V[i,1],(keepsBTCSP$rho[r,1])*sqrt(prod(V[i,])),(keepsBTCSP$rho[r,3])*(keepsBTCSP$sigma_v[r,1])*V[i,1],0,
                    (keepsBTCSP$rho[r,1])*sqrt(prod(V[i,])),V[i,2],0,(keepsBTCSP$rho[r,4])*(keepsBTCSP$sigma_v[r,2])*V[i,2],
                    (keepsBTCSP$rho[r,3])*(keepsBTCSP$sigma_v[r,1])*V[i,1],0,(keepsBTCSP$sigma_v[r,1])^2*V[i,1],(keepsBTCSP$rho[r,2])*prod(keepsBTCSP$sigma_v[r,])*sqrt(prod(V[i,])),
                    0,(keepsBTCSP$rho[r,4])*(keepsBTCSP$sigma_v[r,2])*V[i,2],(keepsBTCSP$rho[r,2])*prod(keepsBTCSP$sigma_v[r,])*sqrt(prod(V[i,])),(keepsBTCSP$sigma_v[r,2])^2*V[i,2]),nrow=4)

  temp <- rtmvnorm(n = 1, 
                   mean = c(x[i,] + keepsBTCSP$mu[r,]+ J[i,],
                            keepsBTCSP$theta[r,] + keepsBTCSP$phi[r,]*(V[i,] - keepsBTCSP$theta[r,])),
                   sigma = Sigma, lower=c(-Inf,-Inf,0,0))
  
  V[i+1,] <- temp[3:4]
  y[i,,r2] <- temp[1:2]
  if( i+1 <= T){ x[i+1,] <- c(0,0)}  
  }
  print(summary(y[,,r2]))
}

y_2d <- y


## SP Only Model
delta <- rep(0,T)
xiy1 <- xiy1s <- rep(0,T)
V <- rep(0,T+1)
J <- x <- rep(0,T)

y <-array(0, dim = c(T, 2, length(Rsequence)))
r2 = 0
for(r in Rsequence){ #loop over posterior draws -->posterior predictive distributions
 r2 = r2 + 1
 print(r)
   V[1] <- c(keepsSP$v[r,1])
  for ( i in 1:T){
  # whichone <- sample(c(0:1),prob=c(keepsSP$lambda[r],1-keepsSP$lambda[r]), 1)
  # delta[i] <- whichone
  # 
  # xiy1s[i] <- rexp(1)
  # xiy1[i] <- (keepsSP$xi_yw[r]*xiy1s[i]+ 
  #               sqrt(xiy1s[i])*keepsSP$xi_yeta[r]*rnorm(1,0,1))
  # 
  J[i] = keepsSP$J[r,i]#(delta[i] == 0)*xiy1[i]
  
  Sigma <- matrix(c(V[i],keepsSP$rho[r]*keepsSP$sigma_v[r]*V[i],keepsSP$rho[r]*keepsSP$sigma_v[r]*V[i],keepsSP$sigma_v[r]^2*V[i]),
                  nrow=2)

  temp <- rtmvnorm(n = 1, 
                   mean = c(x[i] + keepsSP$mu[r] + J[i],
                            keepsSP$theta[r] + keepsSP$phi[r]*(V[i] - keepsSP$theta[r])),
                   sigma = Sigma, lower=c(-Inf,0))
  
  V[i+1] <- temp[2]
  y[i,2,r2] <- temp[1]
  if( i+1 <= T){ x[i+1] <- 0 }  
}
  
}


## BTC Only Model
delta <- rep(0,T)
xiy1 <- xiy1s <- rep(0,T)
V <- rep(0,T+1)
J <- x <- rep(0,T)
r2 = 0
for(r in Rsequence){ #loop over posterior draws -->posterior predictive distributions
  r2 = r2 + 1
  print(r)
  V[1] <- c(keepsBTC$v[r,1])
  for ( i in 1:T){
    # whichone <- sample(c(0:1),prob=c(keepsBTC$lambda[r],1-keepsBTC$lambda[r]), 1)
    # delta[i] <- whichone
    # 
    # xiy1s[i] <- rexp(1)
    # xiy1[i] <- (keepsBTC$xi_yw[r]*xiy1s[i]+ 
    #               sqrt(xiy1s[i])*keepsBTC$xi_yeta[r]*rnorm(1,0,1))
    
    #J[i] = (delta[i] == 0)*xiy1[i]
    J[i] = keepsBTC$J[r,i]
    Sigma <- matrix(c(V[i],keepsBTC$rho[r]*keepsBTC$sigma_v[r]*V[i],keepsBTC$rho[r]*keepsBTC$sigma_v[r]*V[i],keepsBTC$sigma_v[r]^2*V[i]),
                    nrow=2)
    
    temp <- rtmvnorm(n = 1, 
                     mean = c(x[i] + keepsBTC$mu[r] + J[i],
                              keepsBTC$theta[r] + keepsBTC$phi[r]*(V[i] - keepsBTC$theta[r])),
                     sigma = Sigma, lower=c(-Inf,0))
    
    V[i+1] <- temp[2]
    y[i,1,r2] <- temp[1]
    if( i+1 <= T){ x[i+1] <- 0 }  
  }
  
}

y_1d <- y

###------COMPARE 2D VS 2 1D MODELS

QQdatBTCSP = data.frame(S1 = 100*(log(S$`BTC-USD.Close`[-1]) - log(S$`BTC-USD.Close`[-(T+1)])),
                        S2 = 100*(log(S$`GSPC.Close`[-1]) - log(S$`GSPC.Close`[-(T+1)])))

#----basic residual plots

str(y_2d); str(y_1d)

r_2d <- do.call(rbind,apply(y_2d, 3, function(x){x - QQdatBTCSP}) )%>%
  mutate(t = rep(c(1:T), length(Rsequence)),
         v1 = as.vector(keepsBTCSP$v[Rsequence,-1 , 1]),
         v2 = as.vector(keepsBTCSP$v[Rsequence, -1, 2])) %>%
  mutate(S1 = S1/sqrt(v1),
         S2 = S2/sqrt(v2)) %>%
  group_by(t) %>%
  summarise(meanS1 = mean(S1),#posterior summaries of residuals
            meanS2 = mean(S2),
            lowerS1 = quantile(S1, .05),
            upperS1 = quantile(S1, .95),
            lowerS2 = quantile(S2, .05),
            upperS2 = quantile(S2, .95))

r_1d <- do.call(rbind,apply(y_1d, 3, function(x){x - QQdatBTCSP}) )%>%
  mutate(t = rep(c(1:T), length(Rsequence)),
         v1 = as.vector(keepsBTC$v[Rsequence,-1 ]),
         v2 = as.vector(keepsSP$v[Rsequence, -1])) %>%
  mutate(S1 = S1/sqrt(v1),
         S2 = S2/sqrt(v2)) %>%
  group_by(t) %>%
  summarise(meanS1 = mean(S1),
            meanS2 = mean(S2),
            lowerS1 = quantile(S1, .05),
            upperS1 = quantile(S1, .95),
            lowerS2 = quantile(S2, .05),
            upperS2 = quantile(S2, .95))

p1_2d <- r_2d %>%
  ggplot(aes(x =apply(y_2d[,1,], 1, mean), y = meanS1)) +
  geom_point()+
  geom_pointrange(aes(ymin = lowerS1, ymax = upperS1),alpha = I(.2)) +
  theme_bw() +labs(x = "Y-hat", y = "Residual")

p2_2d <-r_2d %>%
  ggplot(aes(x =apply(y_2d[,2,], 1, mean), y = meanS2)) +
  geom_pointrange(aes(ymin = lowerS2, ymax = upperS2),alpha = I(.2))+
  theme_bw()+labs(x = "Y-hat", y = "Residual")


p1_1d <-r_1d %>%
  ggplot(aes(x =apply(y_1d[,1,], 1, mean), y = meanS1))+
  geom_pointrange(aes(ymin = lowerS1, ymax = upperS1),alpha = I(.2))+
  theme_bw()+labs(x = "Y-hat", y = "Residual")


p2_1d <-r_1d %>%
  ggplot(aes(x =apply(y_1d[,2,], 1, mean), y = meanS2)) +
  geom_pointrange(aes(ymin = lowerS2, ymax = upperS2),alpha = I(.2))+
  theme_bw()+labs(x = "Y-hat", y = "Residual")



grid.arrange(p1_2d, p2_2d,p1_1d, p2_1d, nrow = 2)


#count number of jumps larger (in absolute magnitude) than percent%
#y1 is the first time point (true) to get back to original scale
njumps <- function(y, y1, percent=.25){
 jump <-  y/(cumsum(c(0,y[-(T)])) + y1)
  sum(abs(jump) > percent)
}
njumps2 <- function(y, y1, percent=.25){
  jump <-  y/(cumsum(c(0,y[-(T)])) + y1)
  abs(jump) > percent
}



keeps <- list()
keeps3 <- list()

r2 <- 0
for(r in Rsequence){ 
  r2 <- r2 + 1
#H0: that x and y were drawn from the same continuous distribution
#grab the p-values

keeps[[r2]] <- data.frame(y1 = QQdatBTCSP$S1,
                          y2 = QQdatBTCSP$S2,
                          y1_1d = y_1d[,1,r2],
                          y2_1d = y_1d[,2,r2],
                          y1_2d = y_2d[,1,r2],
                          y2_2d = y_2d[,2,r2],
                          iteration = r2)
keeps3[[r2]] <- data.frame(n_y1 = njumps(QQdatBTCSP$S1,log(S$`BTC-USD.Close`[1])),
                           n_y2 = njumps(QQdatBTCSP$S2,log(S$`GSPC.Close`[1]) ),
                           n_both = sum(njumps2(QQdatBTCSP$S1,log(S$`BTC-USD.Close`[1])) &njumps2(QQdatBTCSP$S2,log(S$`GSPC.Close`[1]) ) ),
                          
                           n_1d_1 = njumps(y_1d[,1,r2],log(S$`BTC-USD.Close`[1])),
                           n_1d_2 = njumps(y_1d[,2,r2],log(S$`GSPC.Close`[1])),
                           n_1d_both <- sum(njumps2(y_1d[,1,r2],log(S$`BTC-USD.Close`[1])) & njumps2(y_1d[,2,r2],log(S$`GSPC.Close`[1]) ) ),
                            
                           n_2d_1 = njumps(y_2d[,1,r2],log(S$`BTC-USD.Close`[1])),
                           n_2d_2 = njumps(y_2d[,2,r2],log(S$`GSPC.Close`[1])),
                           n_2d_both = sum(njumps2(y_2d[,1,r2],log(S$`BTC-USD.Close`[1])) & njumps2(y_2d[,2,r2],log(S$`GSPC.Close`[1]) ) ),
                           
                           iteration = r2)
  

}

library(reshape2)


do.call(rbind,keeps3) %>%
  select(c( "n_1d_1", "n_2d_1", "iteration","n_y1", "n_y2"))%>%
  melt(id.var = c("iteration", "n_y1", "n_y2")) %>%
  ggplot() + geom_density(aes(x = value, fill = variable),  alpha = I(.3)) +
  geom_vline(aes(xintercept = n_y1))

(do.call(rbind,keeps3)$n_1d_1 > do.call(rbind,keeps3)$n_y1) %>%mean
(do.call(rbind,keeps3)$n_2d_1 > do.call(rbind,keeps3)$n_y1) %>%mean

do.call(rbind,keeps3) %>%
  select(c( "n_1d_2", "n_2d_2", "iteration","n_y1", "n_y2"))%>%
  melt(id.var = c("iteration", "n_y1", "n_y2")) %>%
  ggplot() + geom_density(aes(x = value, fill = variable),  alpha = I(.3)) +
  geom_vline(aes(xintercept = n_y2))

(do.call(rbind,keeps3)$n_1d_2 > do.call(rbind,keeps3)$n_y2) %>%mean
(do.call(rbind,keeps3)$n_2d_2 > do.call(rbind,keeps3)$n_y2) %>%mean

(do.call(rbind,keeps3)$n_1d_both > do.call(rbind,keeps3)$n_both) %>%mean
(do.call(rbind,keeps3)$n_2d_both > do.call(rbind,keeps3)$n_both) %>%mean


