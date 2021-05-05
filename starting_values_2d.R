#starting values
delta <- sample(c(0:3), T, replace=T, prob=c(0.1,0.1,0.1,0.7))
N_y1 <- as.numeric(delta == 0)
N_y2 <- as.numeric(delta == 1)
N_c <- as.numeric(delta == 2)

#for PGAS mixture density
p <- .5

xi_y1 <- (y[,1]-x[,1])*.8
xi_y2 <- (y[,2]-x[,2])*.8
xi_c <- cbind(xi_y1,xi_y2)
#make sure the s values don't get too big! 
xi_y1s <- .1*sapply(abs(xi_y1),function(x){ifelse(x==0,1,x)}); xi_y2s <- .1*sapply(abs(xi_y2),function(x){ifelse(x==0,1,x)}); xi_cs <- xi_y1s
xi_y1w <- mean(xi_y1); xi_y2w <- mean(xi_y2); xi_cw <- c(mean(xi_y1),mean(xi_y2))
xi_y1eta <- 5; xi_y2eta <- 5;

J = xi_c*(delta==2) +cbind(xi_y1,0)*(delta==0) + cbind(0,xi_y2)*(delta==1) 

rhoc <- .3
sigma_c <- c(5,5)
Sigma_c <- matrix(c(sigma_c[1]^2, 
                    rhoc*prod(sigma_c), 
                    rhoc*prod(sigma_c), 
                    sigma_c[2]^2), ncol = 2)



w_start <- array(0,dim = c(T+1,2))
w_start[1:T,1] <- lowess((y[,1]-mean(y[,1]))^2, f = 0.2)$y #informed starting values
w_start[1:T,2] <- lowess((y[,2]-mean(y[,2]))^2, f = 0.2)$y #informed starting values
w_start[T+1,] <- c(w_start[T,1],w_start[T,2])
#w_start <- c(w_start, tail(w_start,1))
v <- w_start
#plot(v[,1])
#plot(v[,2])

theta <- c(1,1)
phi <- c(0.95,0.95)
rho <- c(0,0,0,0)
sigma_v <- c(0.2,0.2)
mu <- c(0,0)



lambda <- c(0.1,0.1,0.1,0.7)


J_mean <- array(0, dim = c(T,2))
v_mean <- array(0, dim = c(T+1,2))



#keeps is a list that stores kept draws (everything after burnin B)
keeps <- list(
  v = array(NA, dim = c(R, T+1, 2)),
  delta =  array(NA, dim = c(R, T)),
  lambda = array(0, dim = c(R,4)),
  J = array(NA, dim = c(R, T, 2)),
  sigma_v = array(NA, dim=c(R,2)),
  
  sigma_c = array(NA, dim = c(R,2)),
  rhoc = rep(NA, R),
  xi_cw = array(NA, dim = c(R, 2)),
  
  xi_y1eta = rep(NA, R),
  xi_y1w = rep(NA, R),
  
  xi_y2eta = rep(NA, R),
  xi_y2w = rep(NA, R),
  
  phi = array(NA, dim = c(R,2)),
  theta = array(NA, dim = c(R,2)),
  
  mu = array(NA, dim = c(R,2)),
  rho = array(NA, dim = c(R,4))
  
  #xi_c =array(NA, dim = c(R, T,2)),
  #xi_y=array(NA, dim = c(R, T,2))
  
  #xi_cv =array(NA, dim = c(R, T)),
  #delta =array(NA, dim = c(R, T))
  
)