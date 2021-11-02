rm(list = ls()) 
library(ald)
library(ggplot2)
library(grid)
library(gridExtra)
library(dplyr)
library(truncnorm)
library(mvtnorm)
library(tmvtnorm)
library(Rcpp)
library(MCMCpack)
library(quantmod)
library(RcppTN)
library(xtable)
library(tibble)#rownames_to_column()
library(tidyr)
library(reshape2)

domean<- function(x, total){
  R <- total
  if(length(dim(x)) == 2){ #if it's a data frame
    return(apply(x[1:R,], 2, median))
  }else if (length(dim(x)) > 2){
    return(apply(x[1:R,,], 2:3, function(x){(median(x))}))
  }else{
    return(median(x[1:R]))
  }
}

docreds <- function(x, total,q){
  R <- total
  if(length(dim(x)) == 2){ #if it's a data frame
    return(apply(x[1:R,], 2, quantile, q))
  #}#else if (length(dim(x)) > 2){
    #return(apply(x[1:R,,], 2:3, function(x){(median(x))}))
  }else{
    return(quantile(x[1:R],q))
  }
}

library(LaplacesDemon)

total <- 20000 #number of mcmc iterations saved after burn-in, thinning
doESS <- function(x, total){
  R <- total
  if(!is.null(dim(x))){ #if it's a data frame
    return(apply(x[1:R,], 2, ESS))
  }else{
    return(ESS(x[1:R]))
  }
}



############################################################
#----LONG TIME SERIES ONLY---------------------------------
############################################################
getSymbols("BTC-USD",from = "2014-09-15",to = "2020-09-30")
BTC <- as.data.frame(`BTC-USD`)
BTC$Date <- seq(as.Date("2014-09-17"),as.Date("2020-09-30"),by="days")
BTC$`BTC-USD.Close`[BTC$Date=="2020-04-17"] <- 7096.18
getSymbols("^GSPC",from = "2014-09-15",to = "2020-09-30")
SP500 <- as.data.frame(`GSPC`)
SP500$Date <- as.Date(rownames(SP500))
S <- merge(BTC,SP500)
T <- nrow(S) - 1
Date <- S$Date

keeps <- readRDS(paste0("keeps/keepsBTCSP_MALD_threshold.rds"))
names <- lapply(keeps[c(4,6:17)], domean, total = 10) %>%unlist %>%names

m = "MALD"
    #keeps <- readRDS(paste0("keeps_long/keepsBTCSP_",m , ".rds"))
    keeps_v1 <- apply(keeps$v[1:20000,,1], 2, function(x){(mean(x))}) #alternative sv
    keeps_v2 <- apply(keeps$v[1:20000,,2], 2, function(x){(mean(x))}) #sp sv
    
    keeps_j1 <- apply(keeps$J[1:20000,,1], 2, function(x){(mean(x))}) #alternative sv
    keeps_j2 <- apply(keeps$J[1:20000,,2], 2, function(x){(mean(x))}) #sp sv
    # keeps_delta[[j]] <- keeps$delta %>%melt %>%group_by(Var2) %>%
    #   summarise(prop0 = mean(value == 0),prop1 = mean(value == 1),prop2 = mean(value == 2),prop3 = mean(value == 3)) %>%
    #   ungroup()
    temp <- rbind(lapply(keeps[c(4,6:17)], docreds,q=.025, total = 20000) %>%unlist,
                               lapply(keeps[c(4,6:17)], docreds,q=.975, total = 20000) %>%unlist)%>% t() %>%
       as.data.frame %>%
       rownames_to_column(var = "parameter")%>%
       mutate_if(is.numeric, round, 2)  %>%mutate(parameter = gsub(".2.5%", "", parameter))
    keeps_summary <- lapply(keeps[c(4,6:17)], domean, total = 20000) %>%unlist %>%melt()%>%rownames_to_column(var = "parameter")%>%
      merge(temp,  by = "parameter" ) %>%
      mutate(parameter = paste0("$\\", parameter, "$"),
             summary = paste0(round(value,3), ", (", V1, ", ", V2, ")"),
             model = m) %>%
      dplyr::select(parameter, summary)
    
    #rho_probs[[j]] <- keeps$rho[,c(3,4)] %>%apply(2, function(x){mean(x > 0)})
    
    #model[j] <-m
  
  

#TABLE XXX POSTERIOR MEANS OF PARAMETERS------
keeps_summary%>%
  xtable() %>%
  print(sanitize.text.function = function(x){x}, type = "latex", include.rownames=FALSE)

#FIGURE XXX STOCHASTIC VOLATILITY--------
keeps_v1_long <- data.frame(Date=Date,
                       BTC=S$`BTC-USD.Close`,
                       Volatility=keeps_v1) %>% gather("Variable","Value",-Date) %>%
      filter(Date <= "2017-01-01")
  
  keeps_v1_long %>%
  ggplot() +
  geom_line(aes(x = Date, y = Value)) +
    facet_grid(Variable~., scales = "free_y") +
    theme_bw() +
    scale_colour_grey() +
    labs(x = "Date",y="")+ theme(text = element_text(size = 20))
  ggsave("volatility.pdf", height = 10, width = 14) 
  
  
  #FIGURE XXX JUMP SIZES--------
  keeps_j1_long <- data.frame(Date=Date[-1],
                              BTC=log(S$`BTC-USD.Close`[-1]) - 
                                log(S$`BTC-USD.Close`[-(T+1)]),
                              Jumps=keeps_j1/100) %>% gather("Variable","Value",-Date) 
  
  
  keeps_j1_long %>%
    ggplot() +
    geom_line(aes(x = Date, y = Value)) +
    facet_grid(Variable~., scales = "free_y") +
    theme_bw() +
    scale_colour_grey() +
    labs(x = "Date",y="") + theme(text = element_text(size = 20))
  ggsave("jump_sizes.pdf", height = 10, width = 14) 

  
  
  # qqplots-----------
  get_qq <- function(keeps, model){## Function to make the QQ Plots
    delta <- rep(0,T)
    V <- array(0, dim = c(T+1, 2))
    J <- y <- x <- array(0, dim = c(T, 2))
    #V[1,] <- apply(keeps$v[,1,], 2, mean)
    V <- apply(keeps$v, 2:3, mean)
    #overwrite J
    J <- apply(keeps$J, 2:3, mean)
    sim <- 0
    for (t in 1:T){
      print(t)
      set.seed(t + 4922035+ sim)
      delta <- sample(c(0:3),prob=apply(keeps$lambda, 2, mean), 1)
      
      set.seed(15866245 + t + sim)
      if(model == "MVN"){
        B <- 1
      }else{
        B <- rexp(1)
      }
      Sigma <- matrix(c(mean(keeps$sigma_c[,1])^2,
                        mean(keeps$rhoc)*mean(keeps$sigma_c[,1])*mean(keeps$sigma_c[,2]),
                        mean(keeps$rhoc)*mean(keeps$sigma_c[,1])*mean(keeps$sigma_c[,2]),
                        mean(keeps$sigma_c[,2]^2)),
                      nrow=2)
      # xi_c <- apply(keeps$xi_cw, 2, mean)*B+
      #   sqrt(B)*rtmvnorm(n = 1, mean = c(0,0), sigma = Sigma)
      # 
      # if(model == "MVN"){
      #   B <- 1
      # }else{
      #   B <- rexp(1)
      # }
      # xi_y1 <- mean(keeps$xi_y1w)*B + sqrt(B)*rnorm(1,0,mean(keeps$xi_y1eta)) #SHOULD THIS BE SQRT(ETA)?
      # if(model == "MVN"){
      #   B <- 1
      # }else{
      #   B <- rexp(1)
      # }
      # xi_y2 <- mean(keeps$xi_y2w)*B + sqrt(B)*rnorm(1,0,mean(keeps$xi_y2eta)) #SHOULD THIS BE SQRT(ETA)?
      # 
      # 
      # J = xi_c*(delta==2) +cbind(xi_y1,0)*(delta==0) + cbind(0,xi_y2)*(delta==1)
      # 
      Sigma <- matrix(c(V[t,1],
                        mean(keeps$rho[,1])*sqrt(prod(V[t,])),
                        mean(keeps$rho[,3])*mean(keeps$sigma_v[,1])*V[t,1],0,
                        mean(keeps$rho[,1])*sqrt(prod(V[t,])),V[t,2],0,mean(keeps$rho[,4])*mean(keeps$sigma_v[,2])*V[t,2],
                        mean(keeps$rho[,3])*mean(keeps$sigma_v[,1])*V[t,1],0,mean(keeps$sigma_v[,1])^2*V[t,1],mean(keeps$rho[,2])*prod(apply(keeps$sigma_v, 2, mean))*sqrt(prod(V[t,])),
                        0,mean(keeps$rho[,4])*mean(keeps$sigma_v[,2])*V[t,2],mean(keeps$rho[,2])*prod(apply(keeps$sigma_v, 2, mean))*sqrt(prod(V[t,])),mean(keeps$sigma_v[,2])^2*V[t,2]),nrow=4)
      

      
      set.seed(463468+t)
      temp <- rtmvnorm(n = 1,
                       mean = c(apply(keeps$mu, 2 ,mean) + J[t,],
                                apply(keeps$theta,2, mean) + apply(keeps$phi, 2, mean)*(V[t,] - apply(keeps$theta, 2, mean))),
                       sigma = Sigma, lower=c(-Inf,-Inf, 0, 0))
      
      #V[t+1,] <- temp[c(3:4)]
      y[t,] <- temp[c(1,2)]
      if( t+1 <= T){ x[t+1] <- 0 }
    }
    
    QQdat = data.frame(SP=100*(log(S[-1,c( "GSPC.Close") ])-
                                 log(S[-nrow(S),c( "GSPC.Close") ])),
                       BTC=100*(log(S[-1,c( "BTC-USD.Close") ])-
                                  log(S[-nrow(S),c( "BTC-USD.Close") ])),
                       y)
    
    p1 <- ggplot() +
      geom_point(aes(x=quantile(QQdat$BTC,seq(0.01,0.99,0.01)),y=quantile(QQdat$X1,seq(0.01,0.99,0.01)))) +
      geom_abline(slope=1,intercept=0) +
      #xlim(c(-15,15)) + ylim(c(-15,15)) +
      xlab("Actual Quantiles") + ylab("Simulated Quantiles") + theme_bw() + ggtitle(paste0("SV-", model))
    
    p2 <- ggplot() +
      geom_point(aes(x=quantile(QQdat$SP,seq(0.01,0.99,0.01)),y=quantile(QQdat$X2,seq(0.01,0.99,0.01)))) +
      geom_abline(slope=1,intercept=0) +
      #xlim(c(-15,15)) + ylim(c(-15,15)) +
      xlab("Actual Quantiles") + ylab("Simulated Quantiles") + theme_bw() + ggtitle("S&P")
    p3 <- grid.arrange(p1,p2, nrow = 1)
    return(p1 + theme(plot.title = element_text(hjust = 0.5,size = 20)))
  }
  
  model <- "MALD"
  plot1 <- get_qq(keeps,model)

  
  ggsave(paste0("QQ_leverage.pdf"),plot1, width = 12, height = 12)
  

