#keeps is a list that stores kept draws (everything after burnin B)
keeps <- list(
  chain = rep(0,n_chns*R),
  v = array(NA, dim = c(n_chns*R, T+1)),
  delta =  array(NA, dim = c(n_chns*R, T)),
  xi_y = array(NA, dim=c(n_chns*R,T)),
  J = array(NA, dim = c(n_chns*R, T)),
  sigma_v = rep(NA, n_chns*R),
  lambda = rep(NA, n_chns*R),
  
  #sigma_c = array(0, dim = c(R,2)),
  #rhoc = rep(NA, R),
  #xi_cw = array(NA, dim = c(R, 2)),
  
  xi_yeta = rep(NA, n_chns*R),
  xi_yw = rep(NA, n_chns*R),
  
  #xi_veta = rep(NA, R),
  #xi_vw = rep(NA, R),
  
  phi = rep(NA, n_chns*R),
  theta = rep(NA, n_chns*R),
  
  mu = rep(NA, n_chns*R),
  rho = rep(NA, n_chns*R),
  
  lnp_y_mid_ztheta = rep(0,n_chns*R),
  lnp_y_mid_theta = rep(0,n_chns*R)
  
  #xi_cy =array(NA, dim = c(R, T)),
  #xi_cv =array(NA, dim = c(R, T)),
  #delta =array(NA, dim = c(R, T))
)


#(in total, we're running R + B iterations of n_chns independent chains)
for (chn in 1:n_chns){
source("starting_values.R") #initialize values

for (i in 1:(R + B)){
  #update stochastic volatility using pgas
  v <- pgas_v_cpp(y, x, omega=v, J, mu, theta, phi, sigma_v, rho, N=10) %>% as.vector
  #sample xi^y, associated parameters s, eta, w
  xi_y<- pgas_xiy_cpp(y, x, omega=v, mu, theta, phi, sigma_v, rho, xi_y, N_y, xi_yw, xi_yeta, xi_ys, N=10) %>% as.vector()
  J = xi_y*N_y
  xi_ys <- pgas_s_cpp(xi_y, xi_yw, xi_yeta, xi_ys, N=10) %>%as.vector()
  delta <- update_delta(y,x,omega=v,xiy=xi_y,mu,theta,phi,sigma_v,rho,lambda)
  #update N_y,J
  N_y <- as.numeric(delta == 0)
  J = xi_y*N_y
  #updete THETA
  lambda <- update_lambda(c(sum(N_y),T-sum(N_y)),c(2,38))
  mu <- update_mu(y,x,omega=v,J,theta,phi,sigma_v,rho)
  theta <- update_theta(y,x,omega=v,J,mu,phi,sigma_v,rho)
  phi <- update_phi(y,x,omega=v,J,mu,theta,sigma_v,rho)
  sigma_v <- update_sigma_v(y,x,omega=v,J,mu,theta,phi,sigma_v,rho,tune_sigmasq = 0.1)
  rho <- update_rho(y,x,omega=v,J,mu,theta,phi,sigma_v,rho,tune_rhosd = 0.02)
  xi_yeta <- update_eta(xi_y,xi_yw,xi_yeta,xi_ys,0.25)
  xi_yw <- update_w(xi_y,xi_yw,xi_yeta,xi_ys,0.25)

  if (i %% 10 == 0){
     print(paste0("Running Chain: ",chn,", Completed Iteration: ",i))
  #   print(qplot(1:T,true_omega[-(T+1)], geom = "line") + 
  #           geom_line(aes(x = 1:T, y = v[-(T+1)]), colour = "purple") +
  #           geom_line(aes(x = 1:T, y = w_start[-(T+1)]), colour = "red"))
  }
  
  #store after burn in
  if (i > B) {
    
    #mean jump sizes (incorporates binary indicator delta draws)
    #J_mean <- J_mean + J/R
    #mean stochastic volatility
    #v_mean <- v_mean + v/R
    
    j = R*(chn - 1) + i - B
    keeps$chain[j] <- chn
    keeps$sigma_v[j] <- sigma_v
    keeps$v[j,] <- v
    keeps$J[j,] <- J
    keeps$lambda[j] <- lambda[1]
    keeps$mu[j] <- mu
    keeps$phi[j] <- phi
    keeps$theta[j] <- theta
    keeps$rho[j] <- rho
    keeps$xi_yeta[j] <- xi_yeta
    keeps$xi_yw[j] <- xi_yw
    keeps$delta[j,] <- delta
    keeps$xi_y[j,] <- xi_y
    #keeps$lnp_y_mid_ztheta[j] <- dnorm(y,
    #                                 x+mu+J+rho/sigma_v*c(v[-1]-theta-phi*(v[-(T)]-theta),0),
    #                                 sqrt(v*c(rep(1-rho^2,T-1),1)),log=TRUE) %>% sum
    #keeps$lnp_y_mid_theta[j] <- pg_lnp_y_mid_theta(y,x,mu,theta,phi,sigma_v,rho,xi_yw,xi_yeta,lambda[1],N=10000)
  }
}
}

library(reshape2)

# ggplot()+
#   geom_line(aes(x = Var2, y = value, group = Var1), alpha = I(.2), data = melt(keeps$v[seq(1,R, by = 50),-(T-1)]))+
#   geom_line(aes(1:T,v_mean[-(T+1)]), colour = "blue")+
#   geom_line(aes(x = 1:T, y = true_omega[-(T+1)]), colour = "red") 
# ggsave(paste0("Plots/volPlots/volatility_estimates_",s,".pdf"))
# 
# 
# qplot(1:T,true_xiy + true_xic[,1], geom = "line") + 
#   geom_point(aes(x = 1:T, y = J_mean[,1]), colour = "purple")
# ggsave(paste0("Plots/J1Plots/J1_estimates_",s,".pdf"))
# 
# qplot(1:(T-1),true_xiw[-1000] + true_xic[-1000,2], geom = "line") + 
#   geom_point(aes(x = 1:(T-1), y = J_mean[-1000,2]), colour = "purple")
# ggsave(paste0("Plots/J2Plots/J2_estimates_",s,".pdf"))
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
