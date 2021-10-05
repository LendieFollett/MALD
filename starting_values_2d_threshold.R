#starting values
#for PGAS mixture density
p <- .5

w_start <- array(0,dim = c(T,2))
w_start[,1] <- lowess((y[,1]-mean(y[,1]))^2, f = 0.1)$y #informed starting values
w_start[,2] <- lowess((y[,2]-mean(y[,2]))^2, f = 0.1)$y #informed starting values
#w_start[T+1,] <- c(w_start[T,1],w_start[T,2])
#w_start <- c(w_start, tail(w_start,1))
v <- rbind(w_start,c(w_start[T,1],w_start[T,2]))
#plot(v[,1])
#plot(v[,2])

theta <- apply(w_start,2,mean)
phi <- c(0.95,0.95)
rho <- c(0,0,0,0,0)
sigma_v <- c(0.2,0.2)
mu <- c(0,0)

#if (norm_jumps == TRUE){
xi_y1s <- xi_y2s <- rep(1, T)
xi_cs <- rep(1, T)
# } else {
#   xi_y1s <- .1*sapply(abs(xi_y1),function(x){ifelse(x==0,1,x)})
#   xi_y2s <- .1*sapply(abs(xi_y2),function(x){ifelse(x==0,1,x)})
#   xi_cs <- xi_y1s
# }
# if(exp_jumps == TRUE){
xi_y1w <- 0
xi_y2w <- 0
xi_cw <- c(0,0)
# } else {
#   xi_y1w <- mean(xi_y1)
#   xi_y2w <- mean(xi_y2)
#   xi_cw <- c(mean(xi_y1),mean(xi_y2))
# }

xi_y1eta <- 5; xi_y2eta <- 5;

rhoc <- .3
sigma_c <- c(5,5)
Sigma_c <- matrix(c(sigma_c[1]^2, 
                    rhoc*prod(sigma_c), 
                    rhoc*prod(sigma_c), 
                    sigma_c[2]^2), ncol = 2)

s_test <- rexp(10000*T)

xi_y1 <- matrix(((y[,1]-yprim[,1]-mu[1]) / w_start[,1] + xi_y1w / xi_y1eta^2) / (1 / w_start[,1] + 1 / (rexp(1000*T)*xi_y1eta^2)),ncol=1000) %>% apply(1,mean)
xi_y2 <- matrix(((y[,2]-yprim[,2]-mu[2]) / w_start[,2] + xi_y2w / xi_y2eta^2) / (1 / w_start[,2] + 1 / (rexp(1000*T)*xi_y2eta^2)),ncol=1000) %>% apply(1,mean)
xi_c <- cbind(xi_y1,xi_y2)
#make sure the s values don't get too big! 

lambda <- c(0.05,0.05,0.05,0.85)

probs <- cbind(dnorm(y[,1],yprim[,1] + mu[1] + xi_y1,sqrt(w_start[,1])) * dnorm(y[,2],yprim[,2] + mu[2],sqrt(w_start[,2])) * lambda[1],
               dnorm(y[,1],yprim[,1] + mu[1],sqrt(w_start[,1])) * dnorm(y[,2],yprim[,2] + mu[2] + xi_y2,sqrt(w_start[,2])) * lambda[2],
               dnorm(y[,1],yprim[,1] + mu[1] + xi_c[,1],sqrt(w_start[,1])) * dnorm(y[,2],yprim[,2] + mu[2] + xi_c[,2],sqrt(w_start[,2])) * lambda[3],
               dnorm(y[,1],yprim[,1] + mu[1],sqrt(w_start[,1])) * dnorm(y[,2],yprim[,2] + mu[2],sqrt(w_start[,2])) * lambda[4])

delta <- apply(probs,1,function(p){sample(c(0:3),1,prob=p)})

J = xi_c*(delta==2) +cbind(xi_y1,0)*(delta==0) + cbind(0,xi_y2)*(delta==1) 
N_y1 <- as.numeric(delta == 0)
N_y2 <- as.numeric(delta == 1)
N_c <- as.numeric(delta == 2)


J_mean <- array(0, dim = c(T,2))
v_mean <- array(0, dim = c(T+1,2))

if (use_starting_values == TRUE){
  sv <- readRDS("starting_values_MALD.rds")
  v <- sv$v
  theta <- sv$theta
  phi <- sv$phi
  rho <- sv$rho + rnorm(4, 0, .0001)
  sigma_v <- sv$sigma_v+ rnorm(2, 0, .0001)
  mu <- sv$mu
  xi_y1w <- sv$xi_y1w
  xi_y2w <- sv$xi_y2w
  xi_cw <- sv$xi_cw
  xi_y1eta <- sv$xi_y1eta; xi_y2eta <- sv$xi_y2eta;
  
  rhoc <- sv$rhoc +rnorm(1, 0, .0001)
  sigma_c <- sv$sigma_c
  Sigma_c <- matrix(c(sigma_c[1]^2, 
                      rhoc*prod(sigma_c), 
                      rhoc*prod(sigma_c), 
                      sigma_c[2]^2), ncol = 2)
  
  xi_y1 <- sv$xi_y1
  xi_y2 <- sv$xi_y2
  xi_c <- cbind(sv$xi_y1c,sv$xi_y2c)
  
  delta <- sv$delta
  
  J = xi_c*(delta==2) +cbind(xi_y1,0)*(delta==0) + cbind(0,xi_y2)*(delta==1) 
  N_y1 <- as.numeric(delta == 0)
  N_y2 <- as.numeric(delta == 1)
  N_c <- as.numeric(delta == 2)
}

####### KEEPS IS NOW DEFINED IN run_mcmc_2d.R
#keeps is a list that stores kept draws (everything after burnin B)
# keeps <- list(
#   v = array(NA, dim = c(R, T+1, 2)),
#   delta =  array(NA, dim = c(R, T)),
#   lambda = array(0, dim = c(R,4)),
#   J = array(NA, dim = c(R, T, 2)),
#   sigma_v = array(NA, dim=c(R,2)),
#   
#   sigma_c = array(NA, dim = c(R,2)),
#   rhoc = rep(NA, R),
#   xi_cw = array(NA, dim = c(R, 2)),
#   
#   xi_y1eta = rep(NA, R),
#   xi_y1w = rep(NA, R),
#   
#   xi_y2eta = rep(NA, R),
#   xi_y2w = rep(NA, R),
#   
#   phi = array(NA, dim = c(R,2)),
#   theta = array(NA, dim = c(R,2)),
#   
#   mu = array(NA, dim = c(R,2)),
#   rho = array(NA, dim = c(R,4))
#   
#   #xi_c =array(NA, dim = c(R, T,2)),
#   #xi_y=array(NA, dim = c(R, T,2))
#   
#   #xi_cv =array(NA, dim = c(R, T)),
#   #delta =array(NA, dim = c(R, T))
#   
# )