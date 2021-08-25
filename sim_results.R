library(ald)
library(MASS)
library(tmvtnorm)
library(tidyverse)
library(grid)
library(gridExtra)

# mu <- NULL
# theta <- NULL
# kappa <- NULL
# sigmav <- NULL
# rho <- NULL
# 
# muy1 <- NULL
# nuy1 <- NULL
# muy2 <- NULL
# nuy2 <- NULL
# 
# muyc <- NULL
# nuyc <- NULL
# rhoc <- NULL
# 
# lambda <- NULL
# 
# for (i in 1:100){
#   keeps <- readRDS(paste0("keeps/keeps",i,".rds"))
#   mu <- cbind(mu,keeps$mu)
#   theta <- cbind(theta,keeps$theta)
#   kappa <- cbind(kappa,1-keeps$phi)
#   sigmav <- cbind(sigmav,keeps$sigma_v)
#   rho <- cbind(rho,keeps$rho)
#   
#   muy1 <- c(muy1,keeps$xi_y1w)
#   nuy1 <- c(nuy1,keeps$xi_y1eta)
#   muy2 <- c(muy2,keeps$xi_y2w)
#   nuy2 <- c(nuy2,keeps$xi_y2eta)
#   
#   muyc <- cbind(muyc,keeps$xi_cw)
#   nuyc <- cbind(nuyc,keeps$sigma_c)
#   rhoc <- c(rhoc,keeps$rhoc)
#   
#   lambda <- cbind(lambda,keeps$lambda)
# }

s <- 100
keeps <- readRDS(paste0("keeps/keeps",s,".rds")) 
set.seed(s)
sim = runif(1,1, 50000)
print(sim)
idx =235235235
idx2  = 3452
T <- 1500
t <- T + 1
true_omega <- array(0, dim=c(T+1,2))
true_sigma_v <- c(.2,.3) #standard deviation on volatility
true_rho <- c(.1,-.2,-.4,.4) #correlation between y's, volatility's and leverage effects for assets 1 and 2

true_mu <- c(.05,.05)

y <- array(0, dim=c(T,2))
#x <- array(0, dim=c(T,2))
#x[1,] <- 0

#volatility model parameters, sparting point
true_theta <- c(1,1.5)
true_phi <- c(.98,.96)
true_omega[1,] <- c(0.9,1.6)

#contemporaneous jump parameters
true_sigma_c <- c(3.5,3.5)
true_w_c <- c(-3,3)
true_rhoc <- -0.5

#asset 1 jump parameters
true_eta_y1<- 3.5
true_w_y1<- -3

#asset 2 jump parameters
true_eta_y2 <- 3.5
true_w_y2 <- 3


true_Sigma_c <- matrix(c(true_sigma_c[1]*true_sigma_c[1], 
                         true_rhoc*(true_sigma_c[2])*sqrt(true_sigma_c[1]), 
                         true_rhoc*(true_sigma_c[2])*sqrt(true_sigma_c[1]), 
                         true_sigma_c[2]*true_sigma_c[2]), ncol = 2)



true_xic <- matrix(0,nrow = T, ncol = 2)
true_xiy1 <- rep(0, T)
true_xiy2 <- rep(0, T)

true_xics <- rep(0, T)
true_xiy1s <- rep(0, T)
true_xiy2s <- rep(0, T)

true_J<- matrix(0,nrow = T, ncol = 2)

true_delta <- rep(0, T)


for ( i in 1:T){ 
  print(i)
  #true jump distribution
  set.seed(i + 4922035+ sim)
  whichone <- sample(c(0:3),prob = c(0.05,0.05, 0.05,.85), 1)
  true_delta[i] <- whichone
  
  set.seed(i + 2352350+ sim)
  true_xics[i] <- rexp(1)
  set.seed(i + 52352350+ sim)
  X<- mvrnorm(n = 1, mu = c(0,0), Sigma = true_Sigma_c)
  true_xic[i,] <- (true_w_c*true_xics[i]  + sqrt(true_xics[i] )*X)
  
  
  set.seed(15866245 + i+ sim)
  true_xiy1s[i] <- rexp(1)
  true_xiy1[i] <- (true_w_y1*true_xiy1s[i]+ sqrt(true_xiy1s[i])*true_eta_y1*rnorm(1,0,1))
  set.seed(57759235+i)
  true_xiy2s[i] <- rexp(1)
  true_xiy2[i] <- (true_w_y2*true_xiy2s[i]+ sqrt(true_xiy2s[i])*true_eta_y2*rnorm(1,0,1))
  
  true_J[i,] = (true_delta[i] == 2)*true_xic[i,] + (true_delta[i] == 0)*c(true_xiy1[i],0) + (true_delta[i] == 1)*c(0,true_xiy2[i])
  
  Sigma <- matrix(c(true_omega[i,1],true_rho[1]*sqrt(prod(true_omega[i,])),true_rho[3]*true_sigma_v[1]*true_omega[i,1],0,
                    true_rho[1]*sqrt(prod(true_omega[i,])),true_omega[i,2],0,true_rho[4]*true_sigma_v[2]*true_omega[i,2],
                    true_rho[3]*true_sigma_v[1]*true_omega[i,1],0,true_sigma_v[1]^2*true_omega[i,1],true_rho[2]*prod(true_sigma_v)*sqrt(prod(true_omega[i,])),
                    0,true_rho[4]*true_sigma_v[2]*true_omega[i,2],true_rho[2]*prod(true_sigma_v)*sqrt(prod(true_omega[i,])),true_sigma_v[2]^2*true_omega[i,2]),nrow=4)
  idx <- 0
  set.seed(463468+i + idx + sim)
  temp <- rtmvnorm(n = 1, 
                   mean = c(true_mu + true_J[i,],
                            true_theta + true_phi*(true_omega[i,] - true_theta)),
                   sigma = Sigma, lower=c(-Inf,-Inf,0,0))
  
  true_omega[i+1,] <- temp[3:4]
  y[i,] <- temp[1:2]
}

Vol1 <- data.frame(Time=c(1:1501),
                   True=true_omega[,1],
                   Est=keeps$v_mean[,1])
Vol2 <- data.frame(Time=c(1:1501),
                   True=true_omega[,2],
                   Est=keeps$v_mean[,2])
J1 <- data.frame(Time=c(1:1500),
                   True=true_J[,1],
                   Est=keeps$J_mean[,1])
J2 <- data.frame(Time=c(1:1500),
                   True=true_J[,2],
                   Est=keeps$J_mean[,2])
g1 <- ggplot(Vol1) +
  geom_line(aes(x=Time,y=True),color="red") +
  geom_line(aes(x=Time,y=Est),color="blue") +
  ylab("Volatility -- Asset 1")
g2 <- ggplot(Vol2) +
  geom_line(aes(x=Time,y=True),color="red") +
  geom_line(aes(x=Time,y=Est),color="blue") +
  ylab("Volatility -- Asset 2")
g3 <- ggplot(J1) +
  geom_line(aes(x=Time,y=True)) +
  geom_point(aes(x=Time,y=Est),color='purple') +
  ylab("Jump -- Asset 1")
g4 <- ggplot(J2) +
  geom_line(aes(x=Time,y=True)) +
  geom_point(aes(x=Time,y=Est),color='purple') +
  ylab("Jump -- Asset 2")

p1 <- grid.arrange(g1,g3,g2,g4,nrow=2)
library(reshape2)
vol_plot <- rbind(data.frame(Vol1, series = "Asset 1"), 
      data.frame(Vol2, series = "Asset 2")) %>%
  melt(id.var = c("Time", "series")) %>%
  ggplot() + 
  geom_line(aes(x = Time, y =value, colour = variable, linetype =variable )) +
  facet_wrap(series~., nrow = 2, scales = "free") +
  theme_bw() +
  scale_colour_grey(start = .2, end = .6) +
  theme(legend.position = "none") +
  labs(x = "Time", y = "Volatility")
vol_plot
ggsave("sim_volatility.pdf",vol_plot)



jump_plot <- rbind(data.frame(J1, series = "Asset 1"), 
      data.frame(J2, series = "Asset 2")) %>%
  ggplot() + 
  geom_point(aes(x = Time, y =True)) +
  geom_line(aes(x = Time, y =Est), colour = "grey60") +
  facet_wrap(series~., nrow = 2, scales = "free") +
  theme_bw() +
  scale_colour_grey(start = .2, end = .6) +
  theme(legend.position = "none") +
  labs(x = "Time", y = "Volatility")
jump_plot
ggsave("sim_jump.pdf",vol_plot)



get_qq <- function(keeps, simdat){## Function to make the QQ Plots
  delta <- rep(0,T)
  V <- array(0, dim = c(T+1, 2))
  J <- y <- x <- array(0, dim = c(T, 2))
  #V[1,] <- apply(keeps$v[,1,], 2, mean)
  V <- keeps$v_mean
  #overwrite J
  J <- keeps$J_mean
  sim <- 0
  for (t in 1:T){
    print(t)
    set.seed(t + 4922035+ sim)
    #delta <- sample(c(0:3),prob=apply(keeps$lambda, 2, mean), 1)
    
    set.seed(15866245 + t + sim)
    B <- rexp(1)
    Sigma <- matrix(c(mean(keeps$sigma_c[1])^2,
                      mean(keeps$rhoc)*mean(keeps$sigma_c[1])*mean(keeps$sigma_c[2]),
                      mean(keeps$rhoc)*mean(keeps$sigma_c[1])*mean(keeps$sigma_c[2]),
                      mean(keeps$sigma_c[2]^2)),
                    nrow=2)
    # xi_c <- apply(keeps$xi_cw, 2, mean)*B+
    #   sqrt(B)*rtmvnorm(n = 1, mean = c(0,0), sigma = Sigma)
    # 
    # if(model == "MVN"){
    #   B <- 1
    # }else{
    #   B <- rexp(1)
    # }
    # xi_y1 <- mean(keeps$xi_y1w)*B + sqrt(B)*rnorm(1,0,mean(keeps$xi_y1eta)) #SHOULD THIS BE SQRT(ETA)?
    # if(model == "MVN"){
    #   B <- 1
    # }else{
    #   B <- rexp(1)
    # }
    # xi_y2 <- mean(keeps$xi_y2w)*B + sqrt(B)*rnorm(1,0,mean(keeps$xi_y2eta)) #SHOULD THIS BE SQRT(ETA)?
    # 
    # 
    # J = xi_c*(delta==2) +cbind(xi_y1,0)*(delta==0) + cbind(0,xi_y2)*(delta==1)
    # 
    Sigma <- matrix(c(V[t,1],
                      mean(keeps$rho[1])*sqrt(prod(V[t,])),
                      mean(keeps$rho[3])*mean(keeps$sigma_v[1])*V[t,1],0,
                      mean(keeps$rho[1])*sqrt(prod(V[t,])),V[t,2],0,mean(keeps$rho[4])*mean(keeps$sigma_v[2])*V[t,2],
                      mean(keeps$rho[3])*mean(keeps$sigma_v[1])*V[t,1],0,mean(keeps$sigma_v[1])^2*V[t,1],mean(keeps$rho[2])*prod(keeps$sigma_v)*sqrt(prod(V[t,])),
                      0,mean(keeps$rho[4])*mean(keeps$sigma_v[2])*V[t,2],mean(keeps$rho[2])*prod(keeps$sigma_v)*sqrt(prod(V[t,])),mean(keeps$sigma_v[2])^2*V[t,2]),nrow=4)
    
    
    
    set.seed(463468+t)
    temp <- rtmvnorm(n = 1,
                     mean = c(keeps$mu + J[t,],
                              keeps$theta + keeps$phi*(V[t,] - keeps$theta)),
                     sigma = Sigma, lower=c(-Inf,-Inf, 0, 0))
    
    #V[t+1,] <- temp[c(3:4)]
    y[t,] <- temp[c(1,2)]
    if( t+1 <= T){ x[t+1] <- 0 }
  }
  
  QQdat = data.frame(simdat,
                     y)
  names(QQdat) = c("V1","V2","V3","V4")
  
  p1 <- ggplot() +
    geom_point(aes(x=quantile(QQdat$V1,seq(0.01,0.99,0.01)),y=quantile(QQdat$V3,seq(0.01,0.99,0.01)))) +
    geom_abline(slope=1,intercept=0) +
    #xlim(c(-15,15)) + ylim(c(-15,15)) +
    xlab("Actual Quantiles") + ylab("Simulated Quantiles") + theme_bw() + ggtitle("Simulation -- Asset 1")
  
  p2 <- ggplot() +
    geom_point(aes(x=quantile(QQdat$V2,seq(0.01,0.99,0.01)),y=quantile(QQdat$V4,seq(0.01,0.99,0.01)))) +
    geom_abline(slope=1,intercept=0) +
    #xlim(c(-15,15)) + ylim(c(-15,15)) +
    xlab("Actual Quantiles") + ylab("Simulated Quantiles") + theme_bw() + ggtitle("Simulation -- Asset 2")
  p3 <- grid.arrange(p1,p2, nrow = 1)
  return(p3 + theme(plot.title = element_text(hjust = 0.5,size = 20)))
}

get_qq(keeps,y)
