#PARAMETERS
#SVMVN 
  #norm_jumps: =TRUE if jumps should be normally distributed (fix all s params = 1)
#SVLD
  #exp_jumps: =TRUE if jumps should be (symmetric) exponentially distributed (fix all w params = 0)

#keeps is a list that stores kept draws (everything after burnin B)
keeps <- list(
  chain = rep(0,n_chns*R/thin),
  v = array(NA, dim = c(R*n_chns/thin, T+1, 2)),
  delta =  array(NA, dim = c(R*n_chns/thin, T)),
  lambda = array(0, dim = c(R*n_chns/thin,4)),
  J = array(NA, dim = c(R*n_chns/thin, T, 2)),
  sigma_v = array(NA, dim=c(R*n_chns/thin,2)),
  
  sigma_c = array(NA, dim = c(R*n_chns/thin,2)),
  rhoc = rep(NA, R*n_chns/thin),
  xi_cw = array(NA, dim = c(R*n_chns/thin, 2)),
  
  xi_y1eta = rep(NA, R*n_chns/thin),
  xi_y1w = rep(NA, R*n_chns/thin),
  
  xi_y2eta = rep(NA, R*n_chns/thin),
  xi_y2w = rep(NA, R*n_chns/thin),
  
  phi = array(NA, dim = c(R*n_chns/thin,2)),
  theta = array(NA, dim = c(R*n_chns/thin,2)),
  
  mu = array(NA, dim = c(R*n_chns/thin,2)),
  rho = array(NA, dim = c(R*n_chns/thin,4)),
  
  xi_y1 =array(NA, dim = c(R*n_chns/thin, T)),
  xi_y2 =array(NA, dim = c(R*n_chns/thin, T)),
  xi_y1c =array(NA, dim = c(R*n_chns/thin, T)),
  xi_y2c =array(NA, dim = c(R*n_chns/thin, T))
)
h <- 1e-4

for (chn in 1:n_chns){
source("starting_values_2d.R") #initialize values
####### ONLY FOR PFE #########
  #w_start[1,1] = 0.005
#(in total, we're running R + B iterations)
  
  if (ind == TRUE){
   rho <- rep(0,4)
   rhoc <- 0
  }
  
  
for (i in 1:(R + B)){
  print(i)
  #update stochastic volatility using pgas
  v <- pgas_v_cpp(y,yprim,omega=v,J,mu,theta,phi,sigma_v,rho,N=10)
  
  xi_y1<- pgas_xiy1_cpp(y,yprim, omega=v, mu, theta, phi, sigma_v, rho, xi_y1, xi_y2, xi_c, N_y1, N_y2, N_c, xi_y1w, xi_y1eta, xi_y1s, N=10) %>% as.vector()
  J = xi_c*(delta==2) +cbind(xi_y1,0)*(delta==0) + cbind(0,xi_y2)*(delta==1)
  if (norm_jumps == FALSE){
    xi_y1s <- pgas_s_cpp(xi_y1, xi_y1w, xi_y1eta, xi_y1s, N=10) %>%as.vector()
  }
  xi_y2 <- pgas_xiy2_cpp(y,yprim, omega=v, mu, theta, phi, sigma_v, rho, xi_y1, xi_y2, xi_c, N_y1, N_y2, N_c, xi_y2w, xi_y2eta, xi_y2s, N=10) %>% as.vector()
  J = xi_c*(delta==2) +cbind(xi_y1,0)*(delta==0) + cbind(0,xi_y2)*(delta==1)
  if (norm_jumps == FALSE){ 
    xi_y2s <- pgas_s_cpp(xi_y2, xi_y2w, xi_y2eta, xi_y2s, N=10) %>%as.vector()
  }
  if (ind == TRUE){
    xi_c = cbind(xi_y1,xi_y2)
    J = xi_c*(delta==2) +cbind(xi_y1,0)*(delta==0) + cbind(0,xi_y2)*(delta==1)
  } else {
    xi_c <- pgas_xic_cpp(y,yprim, omega=v, mu, theta, phi, sigma_v, rho, xi_y1, xi_y2, xi_c, N_y1, N_y2, N_c, xi_cw, sigma_c, rhoc, xi_cs, N=10)
    J = xi_c*(delta==2) +cbind(xi_y1,0)*(delta==0) + cbind(0,xi_y2)*(delta==1)
    if (norm_jumps == FALSE){
      xi_cs <- pgas_sc_cpp(xi_c, xi_cw, Sigma_c, xi_cs, N=10) %>% as.vector()
    }
  }
  delta <- update_delta(y,yprim,omega=v,xiy1=xi_y1, xiy2=xi_y2, xic=xi_c,mu,theta,phi,sigma_v,rho,lambda)

  N_y1 <- as.numeric(delta == 0)
  N_y2 <- as.numeric(delta == 1)
  N_c <- as.numeric(delta == 2)
  J = xi_c*(delta==2) +cbind(xi_y1,0)*(delta==0) + cbind(0,xi_y2)*(delta==1)
  lambda <- update_lambda(c(sum(N_y1),sum(N_y2),sum(N_c),T-sum(c(N_y1,N_y2,N_c))),c(2,2,2,34))
  
  f <- function(s){(log_pxi(xi_y1,xi_y1s,xi_y1w,s+h/2, w_prior_param = c(0,2.5), eta_prior_param = c(.5, 1)) - 
                      log_pxi(xi_y1,xi_y1s,xi_y1w,s-h/2, w_prior_param = c(0,2.5), eta_prior_param = c(.5, 1))) / h}
  hat <- uniroot(f,c(0.001,100))$root %>% tryCatch(error = function(e) {xi_y1eta})
  sd <- sqrt(-h^2 / (log_pxi(xi_y1,xi_y1s,xi_y1w,hat+h, w_prior_param = c(0,2.5), eta_prior_param = c(.5, 1)) - 
                       2*log_pxi(xi_y1,xi_y1s,xi_y1w,hat, w_prior_param = c(0,2.5), eta_prior_param = c(.5, 1)) + 
                       log_pxi(xi_y1,xi_y1s,xi_y1w,hat-h, w_prior_param = c(0,2.5), eta_prior_param = c(.5, 1))))
  sd <- ifelse(is.nan(sd),0.01,sd)
  xi_y1eta <- update_eta(xi_y1,xi_y1s,xi_y1w,xi_y1eta,hat,sd,w_prior_param = c(0,2.5), eta_prior_param = c(.5, 1))
  
  f <- function(s){(log_pxi(xi_y2,xi_y2s,xi_y2w,s+h/2, w_prior_param = c(0,.25), eta_prior_param = c(.5, 1)) - 
                      log_pxi(xi_y2,xi_y2s,xi_y2w,s-h/2, w_prior_param = c(0,.25), eta_prior_param = c(.5, 1))) / h}
  hat <- uniroot(f,c(0.001,100))$root %>% tryCatch(error = function(e) {xi_y2eta})
  sd <- sqrt(-h^2 / (log_pxi(xi_y2,xi_y2s,xi_y2w,hat+h, w_prior_param = c(0,.25), eta_prior_param = c(.5, 1)) - 
                       2*log_pxi(xi_y2,xi_y2s,xi_y2w,hat, w_prior_param = c(0,.25), eta_prior_param = c(.5, 1)) + 
                       log_pxi(xi_y2,xi_y2s,xi_y2w,hat-h, w_prior_param = c(0,.25), eta_prior_param = c(.5, 1))))
sd <- ifelse(is.nan(sd),0.01,sd)
  xi_y2eta <- update_eta(xi_y2,xi_y2s,xi_y2w,xi_y2eta,hat,sd,w_prior_param = c(0,0.25), eta_prior_param = c(.5, 1))
  
  if(exp_jumps == FALSE){
    f <- function(s){(log_pxi(xi_y1,xi_y1s,s+h/2,xi_y1eta, w_prior_param = c(0,2.5), eta_prior_param = c(.5, 1)) - 
                        log_pxi(xi_y1,xi_y1s,s-h/2,xi_y1eta, w_prior_param = c(0,2.5), eta_prior_param = c(.5, 1))) / h}
    hat <- uniroot(f,c(-100,100))$root %>% tryCatch(error = function(e) {xi_y1w})
    sd <- sqrt(-h^2 / (log_pxi(xi_y1,xi_y1s,hat+h,xi_y1eta, w_prior_param = c(0,2.5), eta_prior_param = c(.5, 1)) - 
                         2*log_pxi(xi_y1,xi_y1s,hat,xi_y1eta, w_prior_param = c(0,2.5), eta_prior_param = c(.5, 1)) + 
                         log_pxi(xi_y1,xi_y1s,hat-h,xi_y1eta, w_prior_param = c(0,2.5), eta_prior_param = c(.5, 1))))
sd <- ifelse(is.nan(sd),0.01,sd)
    xi_y1w <- update_w(xi_y1,xi_y1s,xi_y1w,xi_y1eta,hat,sd, w_prior_param = c(0,2.5), eta_prior_param = c(.5, 1))
    
    f <- function(s){(log_pxi(xi_y2,xi_y2s,s+h/2,xi_y2eta, w_prior_param = c(0,.25), eta_prior_param = c(.5, 1)) - 
                        log_pxi(xi_y2,xi_y2s,s-h/2,xi_y2eta, w_prior_param = c(0,.25), eta_prior_param = c(.5, 1))) / h}
    hat <- uniroot(f,c(-100,100))$root %>% tryCatch(error = function(e) {xi_y2w})
    sd <- sqrt(-h^2 / (log_pxi(xi_y2,xi_y2s,hat+h,xi_y2eta, w_prior_param = c(0,.25), eta_prior_param = c(.5, 1)) - 
                         2*log_pxi(xi_y2,xi_y2s,hat,xi_y2eta, w_prior_param = c(0,.25), eta_prior_param = c(.5, 1)) + 
                         log_pxi(xi_y2,xi_y2s,hat-h,xi_y2eta, w_prior_param = c(0,.25), eta_prior_param = c(.5, 1))))
sd <- ifelse(is.nan(sd),0.01,sd)
    xi_y2w <- update_w(xi_y2,xi_y2s,xi_y2w,xi_y2eta,hat,sd,w_prior_param = c(0,.25), eta_prior_param = c(.5, 1))
    
    if (ind == FALSE){
      f <- function(s){(log_pxi_c(xi_c,xi_cs,c(s+h/2,xi_cw[2]),sigma_c,rhoc) - log_pxi_c(xi_c,xi_cs,c(s-h/2,xi_cw[2]),sigma_c,rhoc)) / h}
      hat <- uniroot(f,c(-100,100))$root %>% tryCatch(error = function(e) {xi_cw[1]})
      sd <- sqrt(-h^2 / (log_pxi_c(xi_c,xi_cs,c(hat+h,xi_cw[2]),sigma_c,rhoc) - 2*log_pxi_c(xi_c,xi_cs,c(hat,xi_cw[2]),sigma_c,rhoc) + log_pxi_c(xi_c,xi_cs,c(hat-h,xi_cw[2]),sigma_c,rhoc)))
      sd <- ifelse(is.nan(sd),0.01,sd)
      xi_cw[1] <- update_w_c(xi_c,xi_cs,xi_cw,sigma_c,rhoc,hat,sd,0)
      
      f <- function(s){(log_pxi_c(xi_c,xi_cs,c(xi_cw[1],s+h/2),sigma_c,rhoc) - log_pxi_c(xi_c,xi_cs,c(xi_cw[1],s-h/2),sigma_c,rhoc)) / h}
      hat <- uniroot(f,c(-100,100))$root %>% tryCatch(error = function(e) {xi_cw[2]})
      sd <- sqrt(-h^2 / (log_pxi_c(xi_c,xi_cs,c(xi_cw[1],hat+h),sigma_c,rhoc) - 2*log_pxi_c(xi_c,xi_cs,c(xi_cw[1],hat),sigma_c,rhoc) + log_pxi_c(xi_c,xi_cs,c(xi_cw[1],hat-h),sigma_c,rhoc)))
      sd <- ifelse(is.nan(sd),0.01,sd)
      xi_cw[2] <- update_w_c(xi_c,xi_cs,xi_cw,sigma_c,rhoc,hat,sd,1)
    }
  }
  
  if (ind == FALSE){
    f <- function(s){(log_pxi_c(xi_c,xi_cs,xi_cw,c(s+h/2,sigma_c[2]),rhoc) - log_pxi_c(xi_c,xi_cs,xi_cw,c(s-h/2,sigma_c[2]),rhoc)) / h}
    hat <- uniroot(f,c(0.001,100))$root %>% tryCatch(error = function(e) {sigma_c[1]})
    sd <- sqrt(-h^2 / (log_pxi_c(xi_c,xi_cs,xi_cw,c(hat+h,sigma_c[2]),rhoc) - 2*log_pxi_c(xi_c,xi_cs,xi_cw,c(hat,sigma_c[2]),rhoc) + log_pxi_c(xi_c,xi_cs,xi_cw,c(hat-h,sigma_c[2]),rhoc)))
    sd <- ifelse(is.nan(sd),0.01,sd)
    sigma_c[1] <- update_sigma_c(xi_c,xi_cs,xi_cw,sigma_c,rhoc,hat,sd,0)
    
    f <- function(s){(log_pxi_c(xi_c,xi_cs,xi_cw,c(sigma_c[1],s+h/2),rhoc) - log_pxi_c(xi_c,xi_cs,xi_cw,c(sigma_c[1],s-h/2),rhoc)) / h}
    hat <- uniroot(f,c(0.001,100))$root %>% tryCatch(error = function(e) {sigma_c[2]})
    sd <- sqrt(-h^2 / (log_pxi_c(xi_c,xi_cs,xi_cw,c(sigma_c[1],hat+h),rhoc) - 2*log_pxi_c(xi_c,xi_cs,xi_cw,c(sigma_c[1],hat),rhoc) + log_pxi_c(xi_c,xi_cs,xi_cw,c(sigma_c[1],hat-h),rhoc)))
    sd <- ifelse(is.nan(sd),0.01,sd)
    sigma_c[2] <- update_sigma_c(xi_c,xi_cs,xi_cw,sigma_c,rhoc,hat,sd,1)
    
    f <- function(s){(log_pxi_c(xi_c,xi_cs,xi_cw,sigma_c,s+h/2) - log_pxi_c(xi_c,xi_cs,xi_cw,sigma_c,s-h/2)) / h}
    hat <- uniroot(f,c(-0.9999,0.9999))$root %>% tryCatch(error = function(e) {rhoc})
    sd <- sqrt(-h^2 / (log_pxi_c(xi_c,xi_cs,xi_cw,sigma_c,hat+h) - 2*log_pxi_c(xi_c,xi_cs,xi_cw,sigma_c,hat) + log_pxi_c(xi_c,xi_cs,xi_cw,sigma_c,hat-h)))
    sd <- ifelse(is.nan(sd),0.01,sd)
    rhoc <- update_rhoc(xi_c,xi_cs,xi_cw,sigma_c,rhoc,hat,sd)
  }
  
  Sigma_c <- matrix(c(sigma_c[1]^2,rhoc*prod(sigma_c),rhoc*prod(sigma_c),sigma_c[2]^2),nrow=2)
  mu <- update_mu(y,yprim,omega=v,J,theta,phi,sigma_v,rho) %>% as.vector
  theta <- update_theta(y,yprim,omega=v,J,mu,theta,phi,sigma_v,rho) %>% as.vector
  phi <- update_phi(y,yprim,omega=v,J,mu,theta,phi,sigma_v,rho) %>% as.vector
  
  f <- function(s){(log_pyv(y,yprim,v,J,mu,theta,phi,c(s+h/2,sigma_v[2]),rho) - log_pyv(y,yprim,v,J,mu,theta,phi,c(s-h/2,sigma_v[2]),rho)) / h}
  hat <- uniroot(f,c(0.001,100))$root %>% tryCatch(error = function(e) {sigma_v[1]})
  sd <- sqrt(-h^2 / (log_pyv(y,yprim,v,J,mu,theta,phi,c(hat+h,sigma_v[2]),rho) - 2*log_pyv(y,yprim,v,J,mu,theta,phi,c(hat,sigma_v[2]),rho) + log_pyv(y,yprim,v,J,mu,theta,phi,c(hat-h,sigma_v[2]),rho)))
  sd <- ifelse(is.nan(sd),0.01,sd)
  sigma_v[1] <- update_sigma_v(y,yprim,v,J,mu,theta,phi,sigma_v,rho,hat,sd,0)
  
  f <- function(s){(log_pyv(y,yprim,v,J,mu,theta,phi,c(sigma_v[1],s+h/2),rho) - log_pyv(y,yprim,v,J,mu,theta,phi,c(sigma_v[1],s-h/2),rho)) / h}
  hat <- uniroot(f,c(0.001,100))$root %>% tryCatch(error = function(e) {sigma_v[2]})
  sd <- sqrt(-h^2 / (log_pyv(y,yprim,v,J,mu,theta,phi,c(sigma_v[1],hat+h),rho) - 2*log_pyv(y,yprim,v,J,mu,theta,phi,c(sigma_v[1],hat),rho) + log_pyv(y,yprim,v,J,mu,theta,phi,c(sigma_v[1],hat-h),rho)))
  sd <- ifelse(is.nan(sd),0.01,sd)
  sigma_v[2] <- update_sigma_v(y,yprim,v,J,mu,theta,phi,sigma_v,rho,hat,sd,1)
  
  if (ind == FALSE){
    f <- function(s){(log_pyv(y,yprim,v,J,mu,theta,phi,sigma_v,c(s+h/2,rho[2:4])) - log_pyv(y,yprim,v,J,mu,theta,phi,sigma_v,c(s-h/2,rho[2:4]))) / h}
    a = rho[2]^2 - 1
    b = -2*prod(rho[-1])
    c = 1 - sum(rho[-1]^2) + prod(rho[3:4])^2
    end = (-b + c(-1,1)*sqrt(b^2 - 4*a*c)) / (2*a)
    hat <- uniroot(f,c(min(end)+0.0001,max(end)-0.0001))$root %>% tryCatch(error = function(e) {rho[1]})
    sd <- sqrt(-h^2 / (log_pyv(y,yprim,v,J,mu,theta,phi,sigma_v,c(hat+h,rho[2:4])) - 2*log_pyv(y,yprim,v,J,mu,theta,phi,sigma_v,c(hat,rho[2:4])) + log_pyv(y,yprim,v,J,mu,theta,phi,sigma_v,c(hat-h,rho[2:4]))))
    sd <- ifelse(is.nan(sd),0.01,sd)
    rho[1] <- update_rho(y,yprim,v,J,mu,theta,phi,sigma_v,rho,hat,sd,0)
    
    f <- function(s){(log_pyv(y,yprim,v,J,mu,theta,phi,sigma_v,c(rho[1],s+h/2,rho[3:4])) - log_pyv(y,yprim,v,J,mu,theta,phi,sigma_v,c(rho[1],s-h/2,rho[3:4]))) / h}
    a = rho[1]^2 - 1
    b = -2*prod(rho[-2])
    c = 1 - sum(rho[-2]^2) + prod(rho[3:4])^2
    end = (-b + c(-1,1)*sqrt(b^2 - 4*a*c)) / (2*a)
    hat <- uniroot(f,c(min(end)+0.0001,max(end)-0.0001))$root %>% tryCatch(error = function(e) {rho[2]})
    sd <- sqrt(-h^2 / (log_pyv(y,yprim,v,J,mu,theta,phi,sigma_v,c(rho[1],hat+h,rho[3:4])) - 2*log_pyv(y,yprim,v,J,mu,theta,phi,sigma_v,c(rho[1],hat,rho[3:4])) + log_pyv(y,yprim,v,J,mu,theta,phi,sigma_v,c(rho[1],hat-h,rho[3:4]))))
    sd <- ifelse(is.nan(sd),0.01,sd)
    rho[2] <- update_rho(y,yprim,v,J,mu,theta,phi,sigma_v,rho,hat,sd,1)
  }
  
  f <- function(s){(log_pyv(y,yprim,v,J,mu,theta,phi,sigma_v,c(rho[1:2],s+h/2,rho[4])) - log_pyv(y,yprim,v,J,mu,theta,phi,sigma_v,c(rho[1:2],s-h/2,rho[4]))) / h}
  a = rho[4]^2 - 1
  b = -2*prod(rho[-3])
  c = 1 - sum(rho[-3]^2) + prod(rho[1:2])^2
  end = (-b + c(-1,1)*sqrt(b^2 - 4*a*c)) / (2*a)
  hat <- uniroot(f,c(min(end)+0.0001,max(end)-0.0001))$root %>% tryCatch(error = function(e) {rho[3]})
  sd <- sqrt(-h^2 / (log_pyv(y,yprim,v,J,mu,theta,phi,sigma_v,c(rho[1:2],hat+h,rho[4])) - 2*log_pyv(y,yprim,v,J,mu,theta,phi,sigma_v,c(rho[1:2],hat,rho[4])) + log_pyv(y,yprim,v,J,mu,theta,phi,sigma_v,c(rho[1:2],hat-h,rho[4]))))
  sd <- ifelse(is.nan(sd),0.01,sd)
  rho[3] <- update_rho(y,yprim,v,J,mu,theta,phi,sigma_v,rho,hat,sd,2)
  
  f <- function(s){(log_pyv(y,yprim,v,J,mu,theta,phi,sigma_v,c(rho[1:3],s+h/2)) - log_pyv(y,yprim,v,J,mu,theta,phi,sigma_v,c(rho[1:3],s-h/2))) / h}
  a = rho[3]^2 - 1
  b = -2*prod(rho[-4])
  c = 1 - sum(rho[-4]^2) + prod(rho[1:2])^2
  end = (-b + c(-1,1)*sqrt(b^2 - 4*a*c)) / (2*a)
  hat <- uniroot(f,c(min(end)+0.0001,max(end)-0.0001))$root %>% tryCatch(error = function(e) {rho[4]})
  sd <- sqrt(-h^2 / (log_pyv(y,yprim,v,J,mu,theta,phi,sigma_v,c(rho[1:3],hat+h)) - 2*log_pyv(y,yprim,v,J,mu,theta,phi,sigma_v,c(rho[1:3],hat)) + log_pyv(y,yprim,v,J,mu,theta,phi,sigma_v,c(rho[1:3],hat-h))))
  sd <- ifelse(is.nan(sd),0.01,sd)
  rho[4] <- update_rho(y,yprim,v,J,mu,theta,phi,sigma_v,rho,hat,sd,3)
  
  #store after burn in
  if (i > B) {
    if (i %% thin == 0){
      #mean jump sizes (incorporates binary indicator delta draws)
      J_mean <- J_mean + J/R
      #mean stochastic volatility
      v_mean <- v_mean + v/R
      
      j = (R*(chn - 1) + i - B)/thin
      print(paste0("Saving Iteration: ",j))
      keeps$sigma_v[j,] <- sigma_v
      keeps$v[j,,] <- v
      keeps$J[j,,] <- J
      keeps$lambda[j,] <- lambda
      keeps$mu[j,] <- mu
      keeps$phi[j,] <- phi
      keeps$theta[j,] <- theta
      keeps$rho[j,] <- rho
      keeps$xi_y1eta[j] <- xi_y1eta
      keeps$xi_y1w[j] <- xi_y1w
      keeps$xi_y2eta[j] <- xi_y2eta
      keeps$xi_y2w[j] <- xi_y2w
      keeps$xi_cw[j,] <- xi_cw
      keeps$sigma_c[j,] <- sigma_c
      keeps$rhoc[j] <- rhoc
      keeps$delta[j,] <- delta
      #keeps$xi_y1[j,] <- xi_y1
      #keeps$xi_y2[j,] <- xi_y2
      #keeps$xi_y1c[j,] <- xi_c[,1]
      #keeps$xi_y2c[j,] <- xi_c[,2]
    }
  }
}
}

# library(reshape2)

# ggplot()+
#   #geom_line(aes(x = Var2, y = value, group = Var1), alpha = I(.2), data = melt(keeps$v[seq(1,R, by = 50),-(T-1)]))+
#   geom_line(aes(1:(T+1),v_mean[,1]), colour = "blue") +
#   geom_line(aes(1:(T+1),true_omega[,1]), colour = "red") 
# ggsave(paste0("Plots2d/vol1Plots/volatility_estimates_",s,".pdf"))
# 
# ggplot()+
#   #geom_line(aes(x = Var2, y = value, group = Var1), alpha = I(.2), data = melt(keeps$v[seq(1,R, by = 50),-(T-1)]))+
#   geom_line(aes(1:(T+1),v_mean[,2]), colour = "blue") +
#   geom_line(aes(1:(T+1),true_omega[,2]), colour = "red")
# ggsave(paste0("Plots2d/vol2Plots/volatility_estimates_",s,".pdf"))
# 
# 
# qplot(1:T,true_J[,1], geom = "line") +
#   geom_point(aes(x = 1:T, y = J_mean[,1]), colour = "purple")
# ggsave(paste0("Plots2d/J1Plots/J1_estimates_",s,".pdf"))
# 
# qplot(1:T,true_J[,2], geom = "line") +
#   geom_point(aes(x = 1:T, y = J_mean[,2]), colour = "purple")
# ggsave(paste0("Plots2d/J2Plots/J2_estimates_",s,".pdf"))
# 
# qplot(true_xic[,1], J_mean[,1]) + geom_abline(aes(slope = 1, intercept = 0))
# 
# qplot((y-(x + mean(keeps$mu) + J_mean[,1]))/sqrt(v_mean[-(T+1)])) +
#   scale_x_continuous(limits = c(-3,3))
# 
# mean(keeps$rhoc)
# mean(keeps$sigma_c[,1])
# mean(keeps$sigma_c[,2])
# mean(keeps$xi_cw[,1])
# mean(keeps$xi_cw[,2])
# 
# mean(keeps$rho)
# 
# mean(keeps$phi)
# 
# mean(keeps$theta)
# 
# plot(apply(keeps$delta, 1, function(x){sum(x == 3)}))
# 
# 
# #tuning: rw parameters should update about 20-40% of the time
# 
# length(unique(keeps$sigma_c[,1]))/R
# length(unique(keeps$rho))/R #needs smaller jump size
# length(unique(keeps$rhoc))/R
