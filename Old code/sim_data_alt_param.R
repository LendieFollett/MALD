<<<<<<< HEAD
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
true_omega <- rep(NA, t)
true_sigma_v <- .2 #standard deviation on volatility
true_rho <- -.4 #correlation between y and (future) volatility

true_mu <- .05

y <- rep(NA, T)
x <- rep(NA, T)
x[1]<-0

#volatility model parameters, sparting point
true_theta <- 1
true_phi <- .98
true_omega[1] <- 0.9

#contemporaneous jump parameters
true_sigma_c <- c(3.5,2)
true_w_c <- c(-3,1)
true_rhoc <- -0.5

#return jump parameters
true_eta_y<- 3.5
true_w_y<- -3

#volatility jump parameters
true_eta_v <- 2
true_w_v <- 1
  
  
true_Sigma_c <- matrix(c(true_sigma_c[1]*true_sigma_c[1], 
                         true_rhoc*(true_sigma_c[2])*sqrt(true_sigma_c[1]), 
                         true_rhoc*(true_sigma_c[2])*sqrt(true_sigma_c[1]), 
                       true_sigma_c[2]*true_sigma_c[2]), ncol = 2)



true_xic <- matrix(0,nrow = T, ncol = 2)
true_xiy <- rep(0, T)
true_xiw <- rep(0, T)

true_J<- matrix(0,nrow = T, ncol = 2)

true_delta <- rep(0, T)


for ( i in 1:T){ 
  print(i)
  #true jump distribution
  set.seed(i + 4922035+ sim)
  whichone <- sample(c(0:3),prob = c(0.05,0.05, 0.05,.85), 1)
  true_delta[i] <- whichone

  # set.seed(i + 2352350+ sim)
  # z <- rexp(1)
  # set.seed(i + 52352350+ sim)
  # X<- mvrnorm(n = 1, mu = c(0,0), Sigma = true_Sigma_c)
  # true_xic[i,] <- (whichone == 2)*(true_w_c*z + sqrt(z)*X)


  set.seed(15866245 + i+ sim)
  z <- rexp(1)
  true_xiy[i] <- (whichone == 0)*(true_w_y*z+ sqrt(z)*true_eta_y*rnorm(1,0,1))
  # set.seed(57759235+i)
  # z <- rexp(1)
  # true_xiw[i] <- (whichone == 1)*(true_w_v*z+ sqrt(z)*true_eta_v*rnorm(1,0,1))
  
  Sigma <- diag(c(sqrt(true_omega[i]),true_sigma_v*sqrt(true_omega[i])))%*%(diag(2)+(1-diag(2))*true_rho)%*%diag(c(sqrt(true_omega[i]),true_sigma_v*sqrt(true_omega[i])))
  idx <- 0
  set.seed(463468+i + idx + sim)
  temp <- mvrnorm(n = 1, 
                  mu = c(x[i] +true_mu + true_xic[i,1] + true_xiy[i] ,
                         true_theta + true_phi*(true_omega[i] - true_theta)+true_xic[i,2] + true_xiw[i]),
                  Sigma = Sigma)
  while(temp[2] < 0){
    set.seed(i + 4922035+ sim+ idx)
    whichone <- sample(c(0:3),prob = c(0.05,0.05, 0.05,.85 ), 1)
    true_delta[i] <- whichone

    # set.seed(i + 2352350+ sim+ idx)
    # z <- rexp(1)
    # set.seed(i + 52352350+ sim+ idx)
    # X<- mvrnorm(n = 1, mu = c(0,0), Sigma = true_Sigma_c)
    # true_xic[i,] <- (whichone == 2)*(true_w_c*z + sqrt(z)*X)


    set.seed(15866245 + i+ sim+ idx)
    z <- rexp(1)
    true_xiy[i] <- (whichone == 0)*(true_w_y*z+ sqrt(z)*true_eta_y*rnorm(1,0,1))
    # set.seed(57759235+i+idx)
    # z <- rexp(1)
    # true_xiw[i] <- (whichone == 1)*(true_w_v*z+ sqrt(z)*true_eta_v*rnorm(1,0,1))
    Sigma <- diag(c(sqrt(true_omega[i]),true_sigma_v*sqrt(true_omega[i])))%*%(diag(2)+(1-diag(2))*true_rho)%*%diag(c(sqrt(true_omega[i]),true_sigma_v*sqrt(true_omega[i])))
    idx = idx + 1
    set.seed(463468+i + idx + sim)
  temp <- mvrnorm(n = 1, 
                  mu = c(x[i] +true_mu + true_xic[i,1] + true_xiy[i] ,
                         true_theta + true_phi*(true_omega[i] - true_theta)+true_xic[i,2] + true_xiw[i]),
                  Sigma = Sigma
  )
  }
  
  true_omega[i+1] <- temp[2]
  y[i] <-temp[1]
  if( i+1 <= T){ x[i+1] <- y[i] }  
}

true_J = true_xic + cbind(true_xiy,0) + cbind(0,true_xiw)


p <-  qplot(2:length(y), true_xic[-1,1]) + 
  geom_segment(aes(x = 2:length(y),xend = 2:length(y),
                   y = true_xic[-1,1], 
                   yend = true_xic[-1,1] + y[-1] - ( y[-length(y)] + true_xic[-1,1]) ))

plot(true_xic[,2]);lines(rep(true_sigma_v, length = length(true_omega)),col="red")
plot(true_xic[,1]); lines(sqrt(true_omega),col="red")

plot(y-x); lines(sqrt(true_omega),col="red")
plot( y, type = "l")
plot(true_omega,type = "l")


=======
library(ald)
library(MASS)

#s represents the simulation index - makes sure each simulated data set is distinct and reproducible.
set.seed(s)
sim = runif(1,1, 50000)
print(sim)
idx =235235235
idx2  = 3452
T <- 1000
t <- T + 1
true_omega <- rep(NA, t)
true_sigma_v <- .1 #standard deviation on volatility
true_rho <- -.4 #correlation between y and (future) volatility

true_mu <- .05

y <- rep(NA, T)
x <- rep(NA, T)
x[1]<-0

#volatility model parameters, sparting point
true_theta <- 1
true_phi <- .98
true_omega[1] <- 0.94

#contemporaneous jump parameters
true_sigma_c <- c(5,2)
true_w_c <- c(-3,1)
true_rhoc <- -0.5

#return jump parameters
true_eta_y<- 5
true_w_y<- -3

#volatility jump parameters
true_eta_v <- 2
true_w_v <- 1
  
  
true_Sigma_c <- matrix(c(true_sigma_c[1]*true_sigma_c[1], 
                         true_rhoc*(true_sigma_c[2])*sqrt(true_sigma_c[1]), 
                         true_rhoc*(true_sigma_c[2])*sqrt(true_sigma_c[1]), 
                       true_sigma_c[2]*true_sigma_c[2]), ncol = 2)



true_xic <- matrix(0,nrow = T, ncol = 2)
true_xiy <- rep(0, T)
true_xiw <- rep(0, T)

true_J<- matrix(0,nrow = T, ncol = 2)

true_delta <- rep(0, T)


for ( i in 1:T){ 
  print(i)
  #true jump distribution
  # set.seed(i + 4922035+ sim)
  # whichone <- sample(c(0:3),prob = c(0.05,0.05, 0.05,.85), 1)
  # true_delta[i] <- whichone
  # 
  # set.seed(i + 2352350+ sim)
  # z <- rexp(1)
  # set.seed(i + 52352350+ sim)
  # X<- mvrnorm(n = 1, mu = c(0,0), Sigma = true_Sigma_c)
  # true_xic[i,] <- (whichone == 2)*(true_w_c*z + sqrt(z)*X)
  # 
  # 
  # set.seed(15866245 + i+ sim)  
  # z <- rexp(1)
  # true_xiy[i] <- (whichone == 0)*(true_w_y*z+ sqrt(z)*true_eta_y*rnorm(1,0,1))
  # set.seed(57759235+i) 
  # z <- rexp(1)
  # true_xiw[i] <- (whichone == 1)*(true_w_v*z+ sqrt(z)*true_eta_v*rnorm(1,0,1))
  
  Sigma <- diag(c(sqrt(true_omega[i]),true_sigma_v*sqrt(true_omega[i])))%*%(diag(2)+(1-diag(2))*true_rho)%*%diag(c(sqrt(true_omega[i]),true_sigma_v*sqrt(true_omega[i])))
  idx <- 0
  set.seed(463468+i + idx + sim)
  temp <- mvrnorm(n = 1, 
                  mu = c(x[i] +true_mu + true_xic[i,1] + true_xiy[i] ,
                         true_theta + true_phi*(true_omega[i] - true_theta)+true_xic[i,2] + true_xiw[i]),
                  Sigma = Sigma)
  while(temp[2] < 0){
    # set.seed(i + 4922035+ sim+ idx)
    # whichone <- sample(c(0:3),prob = c(0.05,0.05, 0.05,.85 ), 1)
    # true_delta[i] <- whichone
    # 
    # set.seed(i + 2352350+ sim+ idx)
    # z <- rexp(1)
    # set.seed(i + 52352350+ sim+ idx)
    # X<- mvrnorm(n = 1, mu = c(0,0), Sigma = true_Sigma_c)
    # true_xic[i,] <- (whichone == 2)*(true_w_c*z + sqrt(z)*X)
    # 
    # 
    # set.seed(15866245 + i+ sim+ idx)  
    # z <- rexp(1)
    # true_xiy[i] <- (whichone == 0)*(true_w_y*z+ sqrt(z)*true_eta_y*rnorm(1,0,1))
    # set.seed(57759235+i+idx) 
    # z <- rexp(1)
    # true_xiw[i] <- (whichone == 1)*(true_w_v*z+ sqrt(z)*true_eta_v*rnorm(1,0,1))  
    Sigma <- diag(c(sqrt(true_omega[i]),true_sigma_v*sqrt(true_omega[i])))%*%(diag(2)+(1-diag(2))*true_rho)%*%diag(c(sqrt(true_omega[i]),true_sigma_v*sqrt(true_omega[i])))
    idx = idx + 1
    set.seed(463468+i + idx + sim)
  temp <- mvrnorm(n = 1, 
                  mu = c(x[i] +true_mu + true_xic[i,1] + true_xiy[i] ,
                         true_theta + true_phi*(true_omega[i] - true_theta)+true_xic[i,2] + true_xiw[i]),
                  Sigma = Sigma
  )
  }
  
  true_omega[i+1] <- temp[2]
  y[i] <-temp[1]
  if( i+1 <= T){ x[i+1] <- y[i] }  
}

true_J = true_xic + cbind(true_xiy,0) + cbind(0,true_xiw)


p <-  qplot(2:length(y), true_xic[-1,1]) + 
  geom_segment(aes(x = 2:length(y),xend = 2:length(y),
                   y = true_xic[-1,1], 
                   yend = true_xic[-1,1] + y[-1] - ( y[-length(y)] + true_xic[-1,1]) ))

plot(true_xic[,2]);lines(rep(true_sigma_v, length = length(true_omega)),col="red")
plot(true_xic[,1]); lines(sqrt(true_omega),col="red")

plot(y-x); lines(sqrt(true_omega),col="red")
plot( y, type = "l")
plot(true_omega,type = "l")


>>>>>>> 25b38294694a191ddeb463ab100d608b8263a58c
