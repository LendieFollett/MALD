#starting values
delta <- sample(c(0:3), T, replace=T)
N_y <- as.numeric(delta == 0)
#N_v <- as.numeric(delta == 1)
#N_c <- as.numeric(delta == 2)

#for PGAS mixture density
p <- .5

xi_y <- (y-x)*.8
#xi_v <- rnorm(T)
#xi_c <- cbind(xi_y,xi_v)
#make sure the s values don't get too big! 
xi_ys <- .1*sapply(abs(xi_y),function(x){ifelse(x==0,1,x)})#; xi_vs <- .1*abs(xi_v); xi_cs <- .1*abs(xi_y)
xi_yw <- mean(xi_ys)#;xi_vw <- 0; xi_cw <- c(mean(xi_ys),0)
xi_yeta <- 5#;xi_veta <- 1;

J = xi_y*N_y

rhoc <- 0
sigma_c <- c(5,1)
Sigma_c <- matrix(c(sigma_c[1]^2, 
                    rhoc*prod(sigma_c), 
                    rhoc*prod(sigma_c), 
                    sigma_c[2]^2), ncol = 2)


#w_start <- rep(var(y-x), T+1) #naive starting values
#or
w_start <- lowess((y-x)^2, f = 0.1)$y #informed starting values
w_start <- c(w_start, tail(w_start,1))
v <- w_start

theta <- mean(w_start)
phi <- .99
rho <- 0
sigma_v <- 0.5
mu <- 0



lambda <- .1


J_mean <- array(0, dim = c(length(xi_y),1))
v_mean <- rep(0, T+1)



#keeps is a list that stores kept draws (everything after burnin B)
keeps <- list(
  v = array(NA, dim = c(R, T+1)),
  delta =  array(NA, dim = c(R, T)),
  lambda = array(0, dim = c(R,2)),
  J = array(NA, dim = c(R, T)),
  sigma_v = rep(NA, R),
  
#  sigma_c = array(0, dim = c(R,2)),
#  rhoc = rep(NA, R),
#  xi_cw = array(NA, dim = c(R, 2)),
  
  xi_yeta = rep(NA, R),
  xi_yw = rep(NA, R),
  
#  xi_veta = rep(NA, R),
#  xi_vw = rep(NA, R),
  
  phi = rep(NA, R),
  theta = rep(NA, R),
  
  mu = rep(NA, R),
  rho = rep(NA, R),
  
 # xi_cy =array(NA, dim = c(R, T)),
#  xi_cv =array(NA, dim = c(R, T)),
  delta =array(NA, dim = c(R, T)),
  
  xi_ys = array(NA, dim = c(R,T)),
  xi_y = array(NA, dim = c(R,T))
  #TO ESTIMATE MARGINAL LIKELKIHOOD p(y_{1:T} | hteta-hat)
 # py = array(0, dim = c(R, T)) #initialize at 0!
  
  
)