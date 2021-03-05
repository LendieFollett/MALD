
#information criterion calculations FOR 2 D ONLY
##---> RDIC and DIC
#where RDIC = D(theta-hat) + 2tr(I(theta-hat)*VCOV(theta-hat))
#and   DIC = -4E(ln(p(y|z,theta))) + 2ln(p(y|z-hat, theta-hat))
library(ggplot2)
library(grid)
library(gridExtra)
library(MASS)
library(tmvtnorm)
library(ald)
library(ggplot2)
library(dplyr)
library(truncnorm)
library(mvtnorm)
library(tmvtnorm)
library(Rcpp)
library(MCMCpack)
library(quantmod)
library(RcppTN)


#download full MCMC results (obtained from run_mcmc(_2d).R with fix = FALSE)
keepsBTCSP <- readRDS("keepsBTCSP.rds")
#need to re-run MCMC with parameters fixed at MAP estimates (find code in empirical_data_study.R)
#this is needed for D(theta-hat), Information matrix
keepsBTCSP_MAP  <- readRDS("keepsBTCSP_MAP.rds") #still need to run

#download data 
getSymbols("BTC-USD",from = "2014-09-17",to = "2020-09-30")
BTC <- as.data.frame(`BTC-USD`)
BTC$Date <- rownames(BTC)
getSymbols("^GSPC",from = "2014-09-17",to = "2020-09-30")
SP500 <- as.data.frame(`GSPC`)
SP500$Date <- rownames(SP500)
S <- merge(BTC,SP500)
T <- nrow(S) - 1
y <- data.frame(
  BTC = log(S$`GSPC.Close`[-1]) - log(S$`GSPC.Close`[-(T+1)]),
  SP = log(S$`GSPC.Close`[-1]) - log(S$`GSPC.Close`[-(T+1)])
)

#-----collect MAP estimates---------
#2d model
THETA_map2d <- list()

THETA_map2d$lambda<-  keepsBTCSP$lambda %>%  apply(2, mean) %>% (function(x){x/sum(x)})
THETA_map2d$sigma_v <- keepsBTCSP$sigma_v %>% apply(2, mean)
THETA_map2d$sigma_c <- keepsBTCSP$sigma_c %>% apply(2, mean)
THETA_map2d$rhoc <- keepsBTCSP$rhoc %>%mean
THETA_map2d$phi <- keepsBTCSP$phi %>%apply(2, mean)
THETA_map2d$theta<- keepsBTCSP$theta %>%apply(2, mean)
THETA_map2d$mu <- keepsBTCSP$mu %>%apply(2, mean)
THETA_map2d$rho<- keepsBTCSP$rho %>%apply(2, mean)
THETA_map2d$xi_yeta <- cbind(keepsBTCSP$xi_y1eta,keepsBTCSP$xi_y2eta )%>%apply(2, mean)
THETA_map2d$xi_yw <- cbind(keepsBTCSP$xi_y1w,keepsBTCSP$xi_y2w )%>%apply(2, mean)








#RDIC CALCULATIONS  (1)------- D(theta-hat) = 2ln(p(y|theta-hat)) = 2*(Q-H)-----
#ln(p(y|theta-hat)) = Q - H


#calculate Q

#calculate H

#RDIC CALCULATIONS  (2)------- 2tr(I(theta-hat)*V(theta-hat))------

#----GET VCOV MATRIX

#THIS NEEDS TO BE IN THE SAME ORDER AS THE DERIVATIVES MATRIX!

VCOV_2d <- keepsBTCSP[c("lambda", "sigma_v", "sigma_c", "rhoc", "phi", "theta", "mu", "rho")] %>%
  as.data.frame() %>%cov()

#---FIRST DERIVATIVE FUNCTIONS----

#---SECOND DERIVATIVE FUNCTIONS----


#Need E(1st ds | y, theta-hat), E(1st ds SQUARED | y, theta-hat), E(2nd ds | y, theta-hat)-----
#here, create functions (of the functions created above) to construct appropriately
#orrdered vectors of 1st derivatives and 2nd derivatives. These will be used in the next section
#----E(S), E(S^2), E(2nd derivative matrix)------


#---DIC calculations-----
