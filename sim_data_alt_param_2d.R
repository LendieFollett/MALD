library(ald)
library(MASS)

#s represents the simulation index - makes sure each simulated data set is distinct and reproducible.
set.seed(s)
sim = runif(1,1, 50000)
print(sim)
idx =235235235
idx2  = 3452
T <- 5000
t <- T + 1
true_omega <- array(0, dim=c(T+1,2))
true_sigma_v <- c(.2,.3) #standard deviation on volatility
true_rho <- c(.4,.5,-.4,.4) #correlation between y's, volatility's and leverage effects for assets 1 and 2

true_mu <- c(.05,.05)

y <- array(0, dim=c(T,2))
x <- array(0, dim=c(T,2))
x[1,] <- 0

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
                  mean = c(x[i,] + true_mu + true_J[i,],
                         true_theta + true_phi*(true_omega[i,] - true_theta)),
                  sigma = Sigma, lower=c(-Inf,-Inf,0,0))
  
  true_omega[i+1,] <- temp[3:4]
  y[i,] <- temp[1:2]
  if( i+1 <= T){ x[i+1,] <- y[i,] }  
}



p <-  qplot(2:length(y), true_xic[-1,1]) + 
  geom_segment(aes(x = 2:length(y),xend = 2:length(y),
                   y = true_xic[-1,1], 
                   yend = true_xic[-1,1] + y[-1] - ( y[-length(y)] + true_xic[-1,1]) ))

plot(true_xic[,2]);lines(rep(true_sigma_v, length = length(true_omega)),col="red")
plot(true_xic[,1]); lines(sqrt(true_omega),col="red")

plot(y-x); lines(sqrt(true_omega),col="red")
plot( y, type = "l")
plot(true_omega,type = "l")


