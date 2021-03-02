
#keeps is a list that stores kept draws (everything after burnin B)
keeps <- list(
  chain = rep(0,n_chns*R),
  v = array(NA, dim = c(R*n_chns, T+1, 2)),
  delta =  array(NA, dim = c(R*n_chns, T)),
  lambda = array(0, dim = c(R*n_chns,4)),
  J = array(NA, dim = c(R*n_chns, T, 2)),
  sigma_v = array(NA, dim=c(R*n_chns,2)),
  
  sigma_c = array(NA, dim = c(R*n_chns,2)),
  rhoc = rep(NA, R*n_chns),
  xi_cw = array(NA, dim = c(R*n_chns, 2)),
  
  xi_y1eta = rep(NA, R*n_chns),
  xi_y1w = rep(NA, R*n_chns),
  
  xi_y2eta = rep(NA, R*n_chns),
  xi_y2w = rep(NA, R*n_chns),
  
  phi = array(NA, dim = c(R*n_chns,2)),
  theta = array(NA, dim = c(R*n_chns,2)),
  
  mu = array(NA, dim = c(R*n_chns,2)),
  rho = array(NA, dim = c(R*n_chns,4)),
  
  xi_y1 =array(NA, dim = c(R*n_chns, T)),
  xi_y2 =array(NA, dim = c(R*n_chns, T)),
  xi_y1c =array(NA, dim = c(R*n_chns, T)),
  xi_y2c =array(NA, dim = c(R*n_chns, T))
  
)



for (chn in 1:n_chns){
  source("starting_values_2d.R") #initialize values
  if(fix == TRUE){
    sigma_c <- MAPs$sigma_c
    rhoc <- MAPs$rhoc
    Sigma_c <- matrix(c(sigma_c[1]^2,rhoc*prod(sigma_c),rhoc*prod(sigma_c),sigma_c[2]^2),nrow=2)
    # #update lambda (R)
    lambda <- MAPs$lambda
    mu <- MAPs$mu
    #update theta
    theta <- MAPs$theta
    phi <- MAPs$phi
    #update sigma_v and rho
    sigma_v <- MAPs$sigma_v
    rho <- MAPs$rho
    
    xi_y1eta <- MAPs$xi_yeta[1]
    xi_y1w <- MAPs$xi_yw[1]
    xi_y2eta <- MAPs$xi_yeta[2]
    xi_y2w <- MAPs$xi_yw[2]
    
  }
#(in total, we're running R + B iterations)
for (i in 1:(R + B)){
  print(i)
  #update stochastic volatility using pgas
  v <- pgas_v_cpp(y,x,omega=v,J,mu,theta,phi,sigma_v,rho,N=10)
  #v <- true_omega
  #sample xi^y1, associated parameters s, eta, w
  xi_y1<- pgas_xiy1_cpp(y, x, omega=v, mu, theta, phi, sigma_v, rho, xi_y1, xi_y2, xi_c, N_y1, N_y2, N_c, xi_y1w, xi_y1eta, xi_y1s, N=10) %>% as.vector()
  J = xi_c*(delta==2) +cbind(xi_y1,0)*(delta==0) + cbind(0,xi_y2)*(delta==1)
  #xi_y1 = true_xiy1
  xi_y1s <- pgas_s_cpp(xi_y1, xi_y1w, xi_y1eta, xi_y1s, N=10) %>%as.vector()
  #xi_y1s <- true_xiy1s
  #xi_y1w <- true_w_y1
  
  #sample xi^y2, associated parameters s, eta, w
  xi_y2 <- pgas_xiy2_cpp(y, x, omega=v, mu, theta, phi, sigma_v, rho, xi_y1, xi_y2, xi_c, N_y1, N_y2, N_c, xi_y2w, xi_y2eta, xi_y2s, N=10) %>% as.vector()
  J = xi_c*(delta==2) +cbind(xi_y1,0)*(delta==0) + cbind(0,xi_y2)*(delta==1)
  #xi_y2 <- true_xiy2
  xi_y2s <- pgas_s_cpp(xi_y2, xi_y2w, xi_y2eta, xi_y2s, N=10) %>%as.vector()
  #xi_y2s <- true_xiy2s
  
  #sample xi^c, associated parameters s, w
  xi_c <- pgas_xic_cpp(y, x, omega=v, mu, theta, phi, sigma_v, rho, xi_y1, xi_y2, xi_c, N_y1, N_y2, N_c, xi_cw, sigma_c, rhoc, xi_cs, N=10)
  J = xi_c*(delta==2) +cbind(xi_y1,0)*(delta==0) + cbind(0,xi_y2)*(delta==1)
  #xi_c <- true_xic
  #  #update w_c, Sigma_c parameters)
  xi_cs <- pgas_sc_cpp(xi_c, xi_cw, Sigma_c, xi_cs, N=10) %>% as.vector()
  #xi_cs <- true_xics
#print(paste0("sigma_c=", sigma_c[1],",", sigma_c[2], ", rhoc=",rhoc,", wc=",xi_cw[1],",", xi_cw[2]))
  
  #update delta (c++)
  delta <- update_delta(y,x,omega=v,xiy1=xi_y1, xiy2=xi_y2, xic=xi_c,mu,theta,phi,sigma_v,rho,lambda)
  #delta <- true_delta
  
  N_y1 <- as.numeric(delta == 0)
  N_y2 <- as.numeric(delta == 1)
  N_c <- as.numeric(delta == 2)
  #print(table(delta, true_delta))
  J = xi_c*(delta==2) +cbind(xi_y1,0)*(delta==0) + cbind(0,xi_y2)*(delta==1)
  
  if (fix == FALSE){
  # #update lambda (R)
  lambda <- update_lambda(c(sum(N_y1),sum(N_y2),sum(N_c),T-sum(c(N_y1,N_y2,N_c))),c(10,10,10,170))
  #lambda <- c(0.05,0.05,0.05,0.85)
  xi_y1eta <- update_eta(xi_y1,xi_y1w,xi_y1eta,xi_y1s,0.25)
  #xi_y1eta <- true_eta_y1
  xi_y1w <- update_w(xi_y1,xi_y1w,xi_y1eta,xi_y1s,0.25)
  xi_y2eta <- update_eta(xi_y2,xi_y2w,xi_y2eta,xi_y2s,0.25)
  #xi_y2eta <- true_eta_y2
  xi_y2w <- update_w(xi_y2,xi_y2w,xi_y2eta,xi_y2s,0.25)
  xi_cw <- update_w_c(xi_c, xi_cw, sigma_c, rhoc, xi_cs, tune_wsd = 0.25) %>% as.vector()
  #xi_cw <- true_w_c
  sigma_c <- update_sigma_c(xi_c, xi_cw, sigma_c, rhoc, xi_cs,tune_wsd = 0.25)
  #sigma_c <- true_sigma_c
  rhoc <- update_rho_c(xi_c, xi_cw, sigma_c, rhoc, xi_cs,tune_wsd = 0.02)
  #rhoc <- true_rhoc
  Sigma_c <- matrix(c(sigma_c[1]^2,rhoc*prod(sigma_c),rhoc*prod(sigma_c),sigma_c[2]^2),nrow=2)
  
  #xi_y2w <- true_w_y2
  #update mu
  mu <- update_mu(y,x,omega=v,J,theta,phi,sigma_v,rho) %>% as.vector
  
  #update theta
  theta <- update_theta(y,x,omega=v,J,mu,theta,phi,sigma_v,rho) %>% as.vector
  
  #update phi
  phi <- update_phi(y,x,omega=v,J,mu,theta,phi,sigma_v,rho) %>% as.vector
  
  # #update sigma_v and rho
  sigma_v <- update_sigma_v(y,x,omega=v,J,mu,theta,phi,sigma_v,rho,tune_sigmasq = 0.05) %>% as.vector
  rho <- update_rho(y,x,omega=v,J,mu,theta,phi,sigma_v,rho,tune_sigmasq = 0.04) %>% as.vector
  #tmp <- update_rho_sigma_v(); sigma_v = tmp[1]; rho = tmp[2];
  #print(paste0("theta=", theta[1],",", theta[2]))
  #print(paste0("phi=", phi[1],",", phi[2]))
  #print(paste0("sigma_v=", sigma_v[1],",", sigma_v[2]))
  #print(paste0("rho=", rho[1],",", rho[2], ",", rho[3],",", rho[4]))
  }
  if (i %% 10 == 0){
     #print(paste0("Simulation: ",s, " Completed Iteration: ",i))
  #   print(qplot(1:T,true_omega[-(T+1)], geom = "line") + 
  #           geom_line(aes(x = 1:T, y = v[-(T+1)]), colour = "purple") +
  #           geom_line(aes(x = 1:T, y = w_start[-(T+1)]), colour = "red"))
  }
  
  #store after burn in
  if (i > B) {
    
    #mean jump sizes (incorporates binary indicator delta draws)
    J_mean <- J_mean + J/R
    #mean stochastic volatility
    v_mean <- v_mean + v/R
    
    j = R*(chn - 1) + i - B
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
    keeps$xi_y1[j,] <- xi_y1
    keeps$xi_y2[j,] <- xi_y2
    keeps$xi_y1c[j,] <- xi_c[,1]
    keeps$xi_y2c[j,] <- xi_c[,2]
  }
}
}

library(reshape2)

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
