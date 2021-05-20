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


keepsIND <- readRDS(paste0(filepath,"keepsBTCSP_IND.rds")) #independence
keepsBTCSP_MALD <- readRDS(paste0(filepath,"keepsBTCSP.rds")) #MALD jump;s
keepsBTCSP_MVN <- readRDS(paste0(filepath,"keepsBTCSP_MVN.rds")) #multivariate normal jumps
keepsBTCSP_LD <- readRDS(paste0(filepath,"keepsBTCSP_LD.rds")) #laplacian jumps

## Generate SP, BTC Data for 1-d, 2-d model
getSymbols("BTC-USD",from = "2014-09-15",to = "2020-09-30")
BTC <- as.data.frame(`BTC-USD`)
BTC$Date <- seq(as.Date("2014-09-17"),as.Date("2020-09-30"),by="days")
BTC$`BTC-USD.Close`[BTC$Date=="2020-04-17"] <- 7096.18
getSymbols("^GSPC",from = "2014-09-15",to = "2020-09-30")
SP500 <- as.data.frame(`GSPC`)
SP500$Date <- as.Date(rownames(SP500))
Rsequence <- seq(1,nrow(keepsBTCSP_MALD$v), by = 100)
S <- merge(BTC,SP500)
T <- nrow(S) - 1
S <- merge(BTC,SP500)
y <- as.matrix(100*(log(S[-1,c("BTC-USD.Close","GSPC.Close")]) - log(S[-nrow(S),c("BTC-USD.Close","GSPC.Close")])))
x <- array(0,dim=dim(y))

#to store posterior predictive draws
y_IND <-array(0, dim = c(T, 2, length(Rsequence)))
y_MALD <-array(0, dim = c(T, 2, length(Rsequence)))
y_MVN <-array(0, dim = c(T, 2, length(Rsequence)))
y_LD <-array(0, dim = c(T, 2, length(Rsequence)))

#2d model
r2 <- 0
for(r in Rsequence){ #loop over posterior draws -->posterior predictive distributions
  print(r)
  r2 <- r2 + 1
  V_MVN <- array(0,dim=c(T+1,2))
  V_MALD <- array(0,dim=c(T+1,2))
  V_LD <-array(0,dim=c(T+1,2))
  V_IND <- array(0,dim=c(T+1,2))
  
  V_MVN <- keepsBTCSP_MVN$v[r,,]
  V_MALD <- keepsBTCSP_MALD$v[r,,]
  V_LD <- keepsBTCSP_LD$v[r,,]
  V_IND <- keepsIND$v[r,,]

  sim <- 0 
  for ( i in 1:T){
  
    #MALD JUMPS
  Sigma <- matrix(c(V_MALD[i,1],(keepsBTCSP_MALD$rho[r,1])*sqrt(prod(V_MALD[i,])),(keepsBTCSP_MALD$rho[r,3])*(keepsBTCSP_MALD$sigma_v[r,1])*V_MALD[i,1],0,
                    (keepsBTCSP_MALD$rho[r,1])*sqrt(prod(V_MALD[i,])),V_MALD[i,2],0,(keepsBTCSP_MALD$rho[r,4])*(keepsBTCSP_MALD$sigma_v[r,2])*V_MALD[i,2],
                    (keepsBTCSP_MALD$rho[r,3])*(keepsBTCSP_MALD$sigma_v[r,1])*V_MALD[i,1],0,(keepsBTCSP_MALD$sigma_v[r,1])^2*V_MALD[i,1],(keepsBTCSP_MALD$rho[r,2])*prod(keepsBTCSP_MALD$sigma_v[r,])*sqrt(prod(V_MALD[i,])),
                    0,(keepsBTCSP_MALD$rho[r,4])*(keepsBTCSP_MALD$sigma_v[r,2])*V_MALD[i,2],(keepsBTCSP_MALD$rho[r,2])*prod(keepsBTCSP_MALD$sigma_v[r,])*sqrt(prod(V_MALD[i,])),(keepsBTCSP_MALD$sigma_v[r,2])^2*V_MALD[i,2]),nrow=4)

  temp <- mvrnorm(n = 1, 
                   mu = c(x[i,] + keepsBTCSP_MALD$mu[r,]+ keepsBTCSP_MALD$J[r,i,]) ,
                            #keepsBTCSP_MALD$theta[r,] + keepsBTCSP_MALD$phi[r,]*(V_MALD[i,] - keepsBTCSP_MALD$theta[r,])),
                   Sigma = Sigma[1:2,1:2])#, lower=c(-Inf,-Inf,0,0))
  
  #V_MALD[i+1,] <- keepsBTCSP_MALD$v[r,i,]#pmax(temp[3:4], c(0.001,0.001))
  y_MALD[i,,r2] <- temp[1:2]
  
  #INDEPENDENCE
  
  Sigma <- matrix(c(V_IND[i,1],(keepsIND$rho[r,1])*sqrt(prod(V_IND[i,])),(keepsIND$rho[r,3])*(keepsIND$sigma_v[r,1])*V_IND[i,1],0,
                    (keepsIND$rho[r,1])*sqrt(prod(V_IND[i,])),V_IND[i,2],0,(keepsIND$rho[r,4])*(keepsIND$sigma_v[r,2])*V_IND[i,2],
                    (keepsIND$rho[r,3])*(keepsIND$sigma_v[r,1])*V_IND[i,1],0,(keepsIND$sigma_v[r,1])^2*V_IND[i,1],(keepsIND$rho[r,2])*prod(keepsIND$sigma_v[r,])*sqrt(prod(V_IND[i,])),
                    0,(keepsIND$rho[r,4])*(keepsIND$sigma_v[r,2])*V_IND[i,2],(keepsIND$rho[r,2])*prod(keepsIND$sigma_v[r,])*sqrt(prod(V_IND[i,])),(keepsIND$sigma_v[r,2])^2*V_IND[i,2]),nrow=4)
  
  temp <- mvrnorm(n = 1, 
                   mu = c(x[i,] + keepsIND$mu[r,]+ keepsIND$J[r,i,]) ,
                            #keepsIND$theta[r,] + keepsIND$phi[r,]*(V_IND[i,] - keepsIND$theta[r,])),
                   Sigma = Sigma[1:2,1:2])#, lower=c(-Inf,-Inf,0,0))
  #V_IND[i+1,] <-  pmax(temp[3:4], c(0.001,0.001))
  y_IND[i,,r2] <- temp[1:2]
  
  #MULTIVARIATE NORMAL JUMPS
  Sigma <- matrix(c(V_MVN[i,1],(keepsBTCSP_MVN$rho[r,1])*sqrt(prod(V_MVN[i,])),(keepsBTCSP_MVN$rho[r,3])*(keepsBTCSP_MVN$sigma_v[r,1])*V_MVN[i,1],0,
                    (keepsBTCSP_MVN$rho[r,1])*sqrt(prod(V_MVN[i,])),V_MVN[i,2],0,(keepsBTCSP_MVN$rho[r,4])*(keepsBTCSP_MVN$sigma_v[r,2])*V_MVN[i,2],
                    (keepsBTCSP_MVN$rho[r,3])*(keepsBTCSP_MVN$sigma_v[r,1])*V_MVN[i,1],0,(keepsBTCSP_MVN$sigma_v[r,1])^2*V_MVN[i,1],(keepsBTCSP_MVN$rho[r,2])*prod(keepsBTCSP_MVN$sigma_v[r,])*sqrt(prod(V_MVN[i,])),
                    0,(keepsBTCSP_MVN$rho[r,4])*(keepsBTCSP_MVN$sigma_v[r,2])*V_MVN[i,2],(keepsBTCSP_MVN$rho[r,2])*prod(keepsBTCSP_MVN$sigma_v[r,])*sqrt(prod(V_MVN[i,])),(keepsBTCSP_MVN$sigma_v[r,2])^2*V_MVN[i,2]),nrow=4)
  
  temp <- mvrnorm(n = 1, 
                   mu = c(x[i,] + keepsBTCSP_MVN$mu[r,]+ keepsBTCSP_MVN$J[r,i,] ) ,
                          #  keepsBTCSP_MVN$theta[r,] + keepsBTCSP_MVN$phi[r,]*(V_MVN[i,] - keepsBTCSP_MVN$theta[r,])),
                   Sigma = Sigma[1:2,1:2])#, lower=c(-Inf,-Inf,0,0))
  #V_MVN[i+1,] <-  pmax(temp[3:4], c(0.001,0.001))
  y_MVN[i,,r2] <- temp[1:2]
  
  #LAPLACIAN JUMPS
  Sigma <- matrix(c(V_LD[i,1],(keepsBTCSP_LD$rho[r,1])*sqrt(prod(V_LD[i,])),(keepsBTCSP_LD$rho[r,3])*(keepsBTCSP_LD$sigma_v[r,1])*V_LD[i,1],0,
                    (keepsBTCSP_LD$rho[r,1])*sqrt(prod(V_LD[i,])),V_LD[i,2],0,(keepsBTCSP_LD$rho[r,4])*(keepsBTCSP_LD$sigma_v[r,2])*V_LD[i,2],
                    (keepsBTCSP_LD$rho[r,3])*(keepsBTCSP_LD$sigma_v[r,1])*V_LD[i,1],0,(keepsBTCSP_LD$sigma_v[r,1])^2*V_LD[i,1],(keepsBTCSP_LD$rho[r,2])*prod(keepsBTCSP_LD$sigma_v[r,])*sqrt(prod(V_LD[i,])),
                    0,(keepsBTCSP_LD$rho[r,4])*(keepsBTCSP_LD$sigma_v[r,2])*V_LD[i,2],(keepsBTCSP_LD$rho[r,2])*prod(keepsBTCSP_LD$sigma_v[r,])*sqrt(prod(V_LD[i,])),(keepsBTCSP_LD$sigma_v[r,2])^2*V_LD[i,2]),nrow=4)
  
  temp <- mvrnorm(n = 1, 
                   mu = c(x[i,] + keepsBTCSP_LD$mu[r,]+ keepsBTCSP_LD$J[r,i,]),
                          #  keepsBTCSP_LD$theta[r,] + keepsBTCSP_LD$phi[r,]*(V_LD[i,] - keepsBTCSP_LD$theta[r,])),
                   Sigma = Sigma[1:2,1:2])#, lower=c(-Inf,-Inf,0,0))
  
  #V_LD[i+1,] <-  pmax(temp[3:4], c(0.001,0.001))
  y_LD[i,,r2] <- temp[1:2]
  
  if( i+1 <= T){ x[i+1,] <- c(0,0)}  
  }
}




###------COMPARE 4 COMPETING MODELS

QQdatBTCSP = data.frame(S1 = 100*(log(S$`BTC-USD.Close`[-1]) - log(S$`BTC-USD.Close`[-(T+1)])),
                        S2 = 100*(log(S$`GSPC.Close`[-1]) - log(S$`GSPC.Close`[-(T+1)])))
library(reshape2)
QQdatBTCSP %>%
  melt() %>%
  mutate(variable = factor(variable, levels = c("S1", "S2"),
                           labels = c("BTC", "SP")))%>%
  ggplot() +
  geom_line(aes(x = as.Date(rep(S$Date[-1], 2)), y = value )) +
  facet_wrap(~variable, nrow = 2)

qplot(QQdatBTCSP[,1], apply(y_LD, 1:2,mean )[,1]) + geom_abline(aes(slope = 1, intercept = 0))
qplot(QQdatBTCSP[,2], apply(y_LD, 1:2,mean )[,2]) + geom_abline(aes(slope = 1, intercept = 0))


y_MALD[,,1] %>%as.data.frame()%>%
  melt() %>%
  mutate(variable = factor(variable, levels = c("V1", "V2"),
                           labels = c("BTC", "SP")))%>%
  ggplot() +
  geom_line(aes(x = as.Date(rep(S$Date[-1], 2)), y = value )) +
  facet_wrap(~variable, nrow = 2)

#make 'lineup'
lu_fun <- function(y, index,title){
as.data.frame.table(y) %>%
  mutate(Date = rep(as.Date(S$Date[-1]), 200)) %>%
  rename(variable = Var2,
           iteration= Var3,
           value = Freq) %>%
  select(-c(Var1))  %>%
  mutate(variable = factor(variable, levels = c("A", "B"),
                                             labels = c("BTC", "SP"))) %>%
  rbind((QQdatBTCSP %>%
           melt() %>%
           mutate(variable = factor(variable, levels = c("S1", "S2"),
                                    labels = c("BTC", "SP")),
                  iteration = "truth",
                  Date = as.Date(rep(S$Date[-1], 2))))) %>%
  subset(variable == index& #only look at one index
           iteration %in% c("truth", toupper(letters)[16:26])&
         Date > "2020-01-01")%>%
  ggplot() +
  geom_line(aes(x = Date, y = value )) +
  facet_wrap(~iteration) +
    scale_x_date(date_labels = "%b-%y",date_breaks = "4 month")+
  ggtitle(title)
  }

grid.arrange(lu_fun(y=y_MALD, index="BTC",title= "SVMALD, BTC"),
             lu_fun(y=y_IND, index="BTC",title= "SVIND, BTC"),
             lu_fun(y=y_MVN, index="BTC",title= "SVMVN, BTC"),
             lu_fun(y=y_LD, index="BTC",title= "SVLD, BTC")
)

grid.arrange(lu_fun(y=y_MALD, index="SP",title= "SVMALD, SP"),
             lu_fun(y=y_IND, index="SP",title= "SVIND, SP"),
             lu_fun(y=y_MVN, index="SP",title= "SVMVN, SP"),
             lu_fun(y=y_LD, index="SP",title= "SVLD, SP")
)


#notes:
#SVIND clearly misses the large negative jump in BTC in early 2020

apply(keepsIND$v, 2:3, mean) %>%
  as.data.frame() %>%
  melt() %>%
  mutate(Date = as.Date(rep(S$Date, 2))) %>%
  ggplot() + geom_line(aes(x = Date, y = value)) +
  facet_wrap(~variable, scales = "free_y", nrow = 2)

#----basic residual plots


r_summary <- function(y_2d){
  do.call(rbind,apply(y_2d, 3, function(x){x - QQdatBTCSP}) )%>%
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
            upperS2 = quantile(S2, .95))}


r_MALD <- r_summary(y_MALD)
r_IND <- r_summary(y_IND)
r_MVN <- r_summary(y_MVN)
r_LD <- r_summary(y_LD)

r_MALD %>%
  ggplot(aes(x =apply(y_MALD[,1,], 1, mean), y = meanS1)) +
  geom_point()+
  geom_pointrange(aes(ymin = lowerS1, ymax = upperS1),alpha = I(.2)) +
  theme_bw() +labs(x = "Y-hat", y = "Residual")

r_MALD %>%
  ggplot(aes(x =apply(y_MALD[,2,], 1, mean), y = meanS2)) +
  geom_pointrange(aes(ymin = lowerS2, ymax = upperS2),alpha = I(.2))+
  theme_bw()+labs(x = "Y-hat", y = "Residual")


r_IND %>%
  ggplot(aes(x =apply(y_IND[,1,], 1, mean), y = meanS1))+
  geom_pointrange(aes(ymin = lowerS1, ymax = upperS1),alpha = I(.2))+
  theme_bw()+labs(x = "Y-hat", y = "Residual")

r_IND %>%
  ggplot(aes(x =apply(y_IND[,1,], 1, mean), y = meanS1))+
  geom_pointrange(aes(ymin = lowerS1, ymax = upperS1),alpha = I(.2))+
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


