postmeans <- rep(0,17)
postrmses <- rep(0,17)

for (s in 1:100){
  filename <- paste0("keeps/keeps",s,".rds")
  keeps <- readRDS(filename)
  postmeans[1] <- postmeans[1] + mean(keeps$mu)/100
  postmeans[2] <- postmeans[2] + mean(keeps$theta)/100
  postmeans[3] <- postmeans[3] + mean(keeps$phi)/100
  postmeans[4] <- postmeans[4] + mean(keeps$sigma_v)/100
  postmeans[5] <- postmeans[5] + mean(keeps$rho)/100
  postmeans[6] <- postmeans[6] + mean(keeps$xi_yw)/100
  postmeans[7] <- postmeans[7] + mean(keeps$xi_yeta)/100
  postmeans[8] <- postmeans[8] + mean(keeps$xi_vw)/100
  postmeans[9] <- postmeans[9] + mean(keeps$xi_veta)/100
  postmeans[10] <- postmeans[10] + mean(keeps$xi_cw[,1])/100
  postmeans[11] <- postmeans[11] + mean(keeps$xi_cw[,2])/100
  postmeans[12] <- postmeans[12] + mean(keeps$sigma_c[,1])/100
  postmeans[13] <- postmeans[13] + mean(keeps$sigma_c[,2])/100
  postmeans[14] <- postmeans[14] + mean(keeps$rhoc)/100
  postmeans[15] <- postmeans[15] + mean(keeps$lambda[,1])/100
  postmeans[16] <- postmeans[16] + mean(keeps$lambda[,2])/100
  postmeans[17] <- postmeans[17] + mean(keeps$lambda[,3])/100
  postrmses[1] <- postrmses[1] + (mean(keeps$mu) - true_mu)^2/100
  postrmses[2] <- postrmses[2] + (mean(keeps$theta) - true_theta)^2/100
  postrmses[3] <- postrmses[3] + (mean(keeps$phi) - true_phi)^2/100
  postrmses[4] <- postrmses[4] + (mean(keeps$sigma_v) - true_sigma_v)^2/100
  postrmses[5] <- postrmses[5] + (mean(keeps$rho) - true_rho)^2/100
  postrmses[6] <- postrmses[6] + (mean(keeps$xi_yw) - true_w_y)^2/100
  postrmses[7] <- postrmses[7] + (mean(keeps$xi_yeta) - true_eta_y)^2/100
  postrmses[8] <- postrmses[8] + (mean(keeps$xi_vw) - true_w_v)^2/100
  postrmses[9] <- postrmses[9] + (mean(keeps$xi_veta) - true_eta_v)^2/100
  postrmses[10] <- postrmses[10] + (mean(keeps$xi_cw[,1]) - true_w_c[1])^2/100
  postrmses[11] <- postrmses[11] + (mean(keeps$xi_cw[,2]) - true_w_c[2])^2/100
  postrmses[12] <- postrmses[12] + (mean(keeps$sigma_c[,1]) - true_sigma_c[1])^2/100
  postrmses[13] <- postrmses[13] + (mean(keeps$sigma_c[,2]) - true_sigma_c[2])^2/100
  postrmses[14] <- postrmses[14] + (mean(keeps$rhoc) - true_rhoc)^2/100
  postrmses[15] <- postrmses[15] + (mean(keeps$lambda[,1]) - 0.05)^2/100
  postrmses[16] <- postrmses[16] + (mean(keeps$lambda[,2]) - 0.05)^2/100
  postrmses[17] <- postrmses[17] + (mean(keeps$lambda[,3]) - 0.05)^2/100
}

