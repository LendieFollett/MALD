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

get_qq_short <- function(keeps,data,i,m){## Function to make the QQ Plots
  #browser()
  delta <- rep(0,T)
  V <- array(0, dim = c(T+1, 2))
  J <- y <- x <- array(0, dim = c(T, 2))
  V[1,] <- apply(keeps$v[,1,], 2, mean)
  
  sim <- 0
  for (t in 1:T){
    print(t)
    set.seed(t + 4922035+ sim)
    delta <- sample(c(0:3),prob=apply(keeps$lambda, 2, mean), 1)
    
    set.seed(15866245 + t + sim)
    if (m == "MVN"){
      B <- 1
    } else {
      B <- rexp(1)
    }
    xi_y1 <- mean(keeps$xi_y1w)*B + sqrt(B)*rnorm(1,0,mean(keeps$xi_y1eta)) #SHOULD THIS BE SQRT(ETA)? No
    if (m == "MVN"){
      B <- 1
    } else {
      B <- rexp(1)
    }
    xi_y2 <- mean(keeps$xi_y2w)*B + sqrt(B)*rnorm(1,0,mean(keeps$xi_y2eta)) #SHOULD THIS BE SQRT(ETA)? No
    if (m == "IND"){
      xi_c = cbind(xi_y1,xi_y2)
    } else {
      if (m == "MVN"){
        B <- 1
      } else {
        B <- rexp(1)
      }
      Sigma <- matrix(c(mean(keeps$sigma_c[,1])^2,
                        mean(keeps$rhoc)*mean(keeps$sigma_c[,1])*mean(keeps$sigma_c[,2]),
                        mean(keeps$rhoc)*mean(keeps$sigma_c[,1])*mean(keeps$sigma_c[,2]),
                        mean(keeps$sigma_c[,2]^2)),
                      nrow=2)
      xi_c <- apply(keeps$xi_cw, 2, mean)*B+
        sqrt(B)*rtmvnorm(n = 1, mean = c(0,0), sigma = Sigma)
    }

    J = xi_c*(delta==2) +cbind(xi_y1,0)*(delta==0) + cbind(0,xi_y2)*(delta==1)
    
    Sigma <- matrix(c(V[t,1],
                      mean(keeps$rho[,1])*sqrt(prod(V[t,])),
                      mean(keeps$rho[,3])*mean(keeps$sigma_v[,1])*V[t,1],0,
                      mean(keeps$rho[,1])*sqrt(prod(V[t,])),V[t,2],0,mean(keeps$rho[,4])*mean(keeps$sigma_v[,2])*V[t,2],
                      mean(keeps$rho[,3])*mean(keeps$sigma_v[,1])*V[t,1],0,mean(keeps$sigma_v[,1])^2*V[t,1],mean(keeps$rho[,2])*prod(apply(keeps$sigma_v, 2, mean))*sqrt(prod(V[t,])),
                      0,mean(keeps$rho[,4])*mean(keeps$sigma_v[,2])*V[t,2],mean(keeps$rho[,2])*prod(apply(keeps$sigma_v, 2, mean))*sqrt(prod(V[t,])),mean(keeps$sigma_v[,2])^2*V[t,2]),nrow=4)
    
    set.seed(463468+t)
    temp <- rtmvnorm(n = 1,
                     mean = c(apply(keeps$mu, 2 ,mean) + J,
                              apply(keeps$theta,2, mean) + apply(keeps$phi, 2, mean)*(V[t,] - apply(keeps$theta, 2, mean))),
                     sigma = Sigma, lower=c(-Inf,-Inf, 0, 0))
    
    V[t+1,] <- temp[c(3:4)]
    y[t,] <- temp[c(1,2)]
    if( t+1 <= T){ x[t+1] <- 0 }
  }
  
  QQdat = cbind(data,y)
  names(QQdat) = c("V1","V2","V3","V4")
  
  p1 <- ggplot() +
    geom_point(aes(x=quantile(QQdat$V1,seq(0.01,0.99,0.01)),y=quantile(QQdat$V3,seq(0.01,0.99,0.01)))) +
    geom_abline(slope=1,intercept=0) +
    #xlim(c(-15,15)) + ylim(c(-15,15)) +
    xlab("Actual Quantiles") + ylab("Simulated Quantiles") + theme_bw() + ggtitle(i)
  return(p1)
}

############################################################
#----SHORT TIME SERIES ONLY---------------------------------
#ALL CRYPTO, MEME STOCKS
############################################################
#load data
BTC <- read.csv("BTC_USD_2020-07-01_2021-06-30-CoinDesk.csv")
BTC <- BTC %>% dplyr::select(!Currency)
names(BTC) <- c("Date","BTC-USD.Close","BTC-USD.Open","BTC-USD.High","BTC-USD.Low")
BTC$Date <- as.Date(BTC$Date)
getSymbols("GME",from = "2020-10-01",to = "2021-06-30")
GME <- as.data.frame(GME)
GME$Date <- as.Date(rownames(GME))
getSymbols("AMC",from = "2020-10-01",to = "2021-06-30")
AMC <- as.data.frame(AMC)
AMC$Date <- as.Date(rownames(AMC))
# getSymbols("DOGE-USD",from = "2020-10-01",to = "2021-05-31")
# DOGE <- as.data.frame(`DOGE-USD`)
# DOGE$Date <- as.Date(rownames(DOGE))
DOGE <- read.csv("DOGE_USD_2020-07-01_2021-06-30-CoinDesk.csv")
DOGE <- DOGE %>% dplyr::select(!Currency)
DOGE$Date <- as.Date(DOGE$Date)
names(DOGE) <- c("Date","DOGE-USD.Close","DOGE-USD.Open","DOGE-USD.High","DOGE-USD.Low")
getSymbols("^GSPC",from = "2020-10-01",to = "2021-06-30")
SP500 <- as.data.frame(`GSPC`)
SP500$Date <- as.Date(rownames(SP500))

S <- BTC %>% merge(GME) %>% merge(AMC) %>% merge(DOGE) %>% merge(SP500)
T <- nrow(S) - 1
Date <- S$Date

keeps <- readRDS(paste0("keeps_","SVMALD" ,"_","BTC", ".rds"))
names <- lapply(keeps[c(4,6:17)], domean, total = 10) %>%unlist %>%names

keeps_creds <- list()

prob_rho_1 <- NULL
prob_rho_2 <- NULL

for (i in c("BTC", "DOGE", "AMC", "GME")){
  if (i == "BTC"){
    QQData <- 100*(log(S[-1,c("BTC-USD.Close","GSPC.Close")]) - 
      log(S[-nrow(S),c("BTC-USD.Close","GSPC.Close")]))
    RawData <- rep(S[,c("BTC-USD.Close")],each=4)
  } else if (i == "DOGE"){
    QQData <- 100*(log(S[-1,c("DOGE-USD.Close","GSPC.Close")]) - 
                     log(S[-nrow(S),c("DOGE-USD.Close","GSPC.Close")]))
    RawData <- rep(S[,c("DOGE-USD.Close")],each=4)
  } else if (i == "AMC"){
    QQData <- 100*(log(S[-1,c("AMC.Close","GSPC.Close")]) - 
                     log(S[-nrow(S),c("AMC.Close","GSPC.Close")]))
    RawData <- rep(S[,c("AMC.Close")],each=4)
  } else {
    QQData <- 100*(log(S[-1,c("GME.Close","GSPC.Close")]) - 
                     log(S[-nrow(S),c("GME.Close","GSPC.Close")]))
    RawData <- rep(S[,c("GME.Close")],each=4)
  }
  RawData <- as.data.frame(RawData)
  names(RawData) <- "Price"
  RawData$Date <- rep(Date,each=4)
  j = 0
  keeps_summary <- array(dim = c(4, length(names)))
  keeps_v1 <- array(dim = c(4,dim(keeps$v)[2]))
  keeps_v2 <- array(dim = c(4,dim(keeps$v)[2]))
  model <- rep(NA, 4)
  data <- rep(i, 4)
  colnames(keeps_summary) <- names
  for (m in c("SVIND", "SVLD", "SVMVN", "SVMALD")){
    j = j + 1
    keeps <- readRDS(paste0("keeps_",m ,"_",i, ".rds"))
    keeps_summary[j,] <- lapply(keeps[c(4,6:17)], domean, total = 20000) %>%unlist
    keeps_creds[[j]] <- rbind(lapply(keeps[c(4,6:17)], docreds, q=.025, total = 20000) %>%unlist,
                             lapply(keeps[c(4,6:17)], docreds, q=.975, total = 20000) %>%unlist)%>% t() %>%
      as.data.frame %>%
      rownames_to_column(var = "parameter")%>%
      mutate_if(is.numeric, round, 3) 
    keeps_v1[j,] <- apply(keeps$v[1:20000,,1], 2, function(x){(mean(x))}) #alternative sv
    keeps_v2[j,] <- apply(keeps$v[1:20000,,2], 2, function(x){(mean(x))}) #sp sv
    model[j] <- m
    keeps_summary[j,] <- paste0(round(as.numeric(keeps_summary[j,]),2)," (",keeps_creds[[j]][,2],",",keeps_creds[[j]][,3],")")
    plot <- get_qq_short(keeps,QQData,i,m)
    ggsave(paste0("QQ_Plots_Short_Time/QQ_", m,"_",i, ".pdf"),plot, width = 12, height = 6)
    #prob_rho_1 <- c(prob_rho_1,length(which(keeps$rho[,3] > 0))/20000)
    #prob_rho_2 <- c(prob_rho_2,length(which(keeps$rho[,4] > 0))/20000)
    # if (m == "SVMALD"){
    #   lapply(keeps[c(4,6:17)], doESS, total = 20000) %>% str()
    #   pdf(file=paste0("Convergence Check Plots/xi_y1w_",i,".pdf"))
    #   plot(keeps$xi_y1w,type="l")
    #   dev.off()
    #   pdf(file=paste0("Convergence Check Plots/xi_y2w_",i,".pdf"))
    #   plot(keeps$xi_y2w,type="l")
    #   dev.off()
    #   pdf(file=paste0("Convergence Check Plots/xi_y1eta_",i,".pdf"))
    #   plot(keeps$xi_y1eta,type="l")
    #   dev.off()
    #   pdf(file=paste0("Convergence Check Plots/xi_y2eta_",i,".pdf"))
    #   plot(keeps$xi_y2eta,type="l")
    #   dev.off()
    #   pdf(file=paste0("Convergence Check Plots/xi_cw1_",i,".pdf"))
    #   plot(keeps$xi_cw[,1],type="l")
    #   dev.off()
    #   pdf(file=paste0("Convergence Check Plots/xi_cw2_",i,".pdf"))
    #   plot(keeps$xi_cw[,2],type="l")
    #   dev.off()
    #   pdf(file=paste0("Convergence Check Plots/sigma_c1_",i,".pdf"))
    #   plot(keeps$sigma_c[,1],type="l")
    #   dev.off()
    #   pdf(file=paste0("Convergence Check Plots/sigma_c2_",i,".pdf"))
    #   plot(keeps$sigma_c[,2],type="l")
    #   dev.off()
    # }
  }
  RawData$model = rep(model,T+1)
  keeps_summary <- keeps_summary %>% as.data.frame()%>%
    #mutate_all(round, digits = 2) %>%
    #mutate(lambda = paste0("(",lambda1,", ", lambda2,", ", lambda3,", ", lambda4,")"))%>%
    #dplyr::select(-c(lambda1, lambda2, lambda3, lambda4)) %>%
    mutate(model = model) %>%t()
  
  #TABLE XXX POSTERIOR MEANS OF PARAMETERS------
  keeps_summary%>%
    xtable() %>%
    print()
  # 
  #FIGURE XXX STOCHASTIC VOLATILITY--------
  #ALTERNATIVE CURRENCY
  keeps_v1_long <- keeps_v1 %>% as.data.frame()%>%
    mutate(model = model,
           series = data) %>%
    melt(id.vars = c("model", "series")) %>%
    mutate(Date = rep(Date, each = 4)) #CHECK THIS STRUCTURE
  #S&P
  keeps_v2_long <- keeps_v2 %>%as.data.frame()%>%
    mutate(model = model,
           series = data) %>%
    melt(id.vars = c("model", "series")) %>%
    mutate(Date = rep(Date, each = 4), #CHECK THIS STRUCTURE
           series = "S&P")

  #not sure the best way to display this, may need to modify
  g1 <- ggplot() +
    geom_line(data = RawData,aes(x=Date,y=Price,linetype=model)) +
    xlab("") + ylab("Price") + theme_bw()
  
  g2 <- keeps_v1_long %>%
    ggplot() +
    geom_line(aes(x = Date, y = value, linetype = model)) +
    #facet_grid(series~model., scales = "free_y") +
    ylab("Volatility") +
    theme_bw() +
    scale_colour_grey()
  
  G1 <- ggplotGrob(g1)
  G2 <- ggplotGrob(g2)
  plotname <- paste0(i,"_Short_Vol.pdf")
  ggsave(plotname,grid.draw(rbind(G1,G2)))
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

keeps <- readRDS(paste0("keeps_long/keepsBTCSP" ,"_","MALD", ".rds"))
names <- lapply(keeps[c(4,6:17)], domean, total = 10) %>%unlist %>%names

keeps_summary0 <- list()
keeps_v1 <- array(dim = c(4,dim(keeps$v)[2]))
keeps_v2 <- array(dim = c(4,dim(keeps$v)[2]))

keeps_j1 <- array(dim = c(4,dim(keeps$J)[2]))
keeps_j2 <- array(dim = c(4,dim(keeps$J)[2]))
model <- rep(NA, 4)
keeps_delta <- list()
keeps_creds <- list()
rho_probs <- list()

j = 0
for (m in c("MALD", "IND", "MVN", "LD")){# "SVMVN", "SVLD",
    j = j + 1
    keeps <- readRDS(paste0("keeps_long/keepsBTCSP_",m , ".rds"))
    keeps_v1[j,] <- apply(keeps$v[1:20000,,1], 2, function(x){(mean(x))}) #alternative sv
    keeps_v2[j,] <- apply(keeps$v[1:20000,,2], 2, function(x){(mean(x))}) #sp sv
    
    keeps_j1[j,] <- apply(keeps$J[1:20000,,1], 2, function(x){(mean(x))}) #alternative sv
    keeps_j2[j,] <- apply(keeps$J[1:20000,,2], 2, function(x){(mean(x))}) #sp sv
    keeps_delta[[j]] <- keeps$delta %>%melt %>%group_by(Var2) %>%
      summarise(prop0 = mean(value == 0),prop1 = mean(value == 1),prop2 = mean(value == 2),prop3 = mean(value == 3)) %>%
      ungroup()
    temp <- rbind(lapply(keeps[c(4,6:17)], docreds,q=.025, total = 20000) %>%unlist,
                              lapply(keeps[c(4,6:17)], docreds,q=.975, total = 20000) %>%unlist)%>% t() %>%
      as.data.frame %>%
      rownames_to_column(var = "parameter")%>%
      mutate_if(is.numeric, round, 2)  %>%mutate(parameter = gsub(".2.5%", "", parameter))
    keeps_summary0[[j]] <- lapply(keeps[c(4,6:17)], domean, total = 20000) %>%unlist %>%melt()%>%rownames_to_column(var = "parameter")%>%
      merge(temp,  by = "parameter" ) %>%
      mutate(parameter = paste0("$\\", parameter, "$"),
             summary = paste0(round(value,3), ", (", V1, ", ", V2, ")"),
             model = m) %>%
      dplyr::select(model, parameter, summary)
    
    rho_probs[[j]] <- keeps$rho[,c(3,4)] %>%apply(2, function(x){mean(x > 0)})
    
    model[j] <-m
}


keeps_summary <- do.call(rbind,keeps_summary0) %>% spread(model, summary) %>%
  dplyr::select(parameter, IND, LD, MVN, MALD)
  
do.call(rbind, rho_probs) %>%t() %>%as.data.frame() %>%
  dplyr::rename(MALD = V1, IN = V2 , MVN = V3, LD = V4)
  
  

#TABLE XXX POSTERIOR MEANS OF PARAMETERS------
keeps_summary%>%
  xtable() %>%
  print(sanitize.text.function = function(x){x}, type = "latex", include.rownames=FALSE)

#FIGURE XXX STOCHASTIC VOLATILITY--------
keeps_v1_long <- keeps_v1 %>%as.data.frame()%>%
  mutate(model = model) %>%
  melt(id.vars = c("model")) %>%
  mutate(Date = rep(Date, each = 4), #change when add more models
         series = "BTC")

keeps_v2_long <- keeps_v2 %>%as.data.frame()%>%
  mutate(model = model) %>%
  melt(id.vars = c("model")) %>%
  mutate(Date = rep(Date, each = 4), #change when add more models
         series = "S&P")
  
  
  rbind(keeps_v1_long,keeps_v2_long) %>%
    mutate(Model = factor(model, levels = c("MALD", "IND", "MVN", "LD"), labels = c("SV-MALD", "SV-IND", "SV-MVN", "SV-LD"))) %>%
  ggplot() +
  geom_line(aes(x = Date, y = value, linetype = Model)) +
    facet_grid(series~., scales = "free_y") +
    theme_bw() +
    scale_colour_grey() +
    labs(x = "Date", y = "Volatility")
  ggsave("volatility.pdf", height = 8, width = 10) 
  
  tempdat = keeps_v1_long %>% subset(model == "MALD") %>%
    merge(data.frame(BTC=100*(log(S[-1,c( "BTC-USD.Close") ])-
                                log(S[-nrow(S),c( "BTC-USD.Close") ])),Date=Date[-1]) , by = "Date")
  
p1 <-   ggplot() +
  geom_rect(aes(xmin=as.Date("2016-01-01"),
                xmax=as.Date("2016-12-31"),
                ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue") +
  geom_rect(aes(xmin=as.Date("2017-07-01"),
                xmax=as.Date("2018-03-31"),
                ymin=-Inf,ymax=Inf),alpha=I(0.2),fill="blue") +
  geom_rect(aes(xmin=as.Date("2020-02-01"),
                xmax=as.Date("2021-01-31"),
                ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue")+
    geom_line(aes(x = Date, y = BTC), alpha= I(.5), data = tempdat) +
    geom_line(aes(x = Date, y = value^.5), data = tempdat) +
    theme_bw() +labs(x = "")

tempdat= keeps_v2_long %>% subset(model == "MALD") %>%
  merge(data.frame(SP=100*(log(S[-1,c( "GSPC.Close") ])-
                              log(S[-nrow(S),c( "GSPC.Close") ])),Date=Date[-1]) , by = "Date") 
p2 <- ggplot() +
  geom_rect(aes(xmin=as.Date("2016-01-01"),
                xmax=as.Date("2016-12-31"),
                ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue") +
  geom_rect(aes(xmin=as.Date("2017-07-01"),
                xmax=as.Date("2018-03-31"),
                ymin=-Inf,ymax=Inf),alpha=I(0.2),fill="blue") +
  geom_rect(aes(xmin=as.Date("2020-02-01"),
                xmax=as.Date("2021-01-31"),
                ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue")+
  geom_line(aes(x = Date, y = SP), alpha= I(.5), data = tempdat) +
  geom_line(aes(x = Date, y = value^.5), data = tempdat) +
  theme_bw()
  

p3 <-grid.arrange(p1, p2, nrow = 2)  
p3
ggsave("leverage.pdf",p3, height = 8, width = 10) 



#FIGURE XXX JUMP TYPES OVER TIME-------------------
  
  keeps_delta[[1]] %>%
    mutate(Date = Date[-1]) %>%
    dplyr::select(-c(Var2))%>%
    melt(id.vars = c("Date")) %>%
    mutate(Jump_Type = factor(variable, levels = c("prop0", "prop1", "prop2", "prop3"), 
                              labels = c("BTC only", "S&P only", "Both", "None"))) %>%
    ggplot() +
    geom_line(aes(x = Date, y = value, linetype = Jump_Type)) +
    facet_wrap(~Jump_Type) +
    theme_bw() +
    labs(y = "Posterior Probability")+
    theme(legend.position = "none")
  ggsave("jump_type_SVMALD.pdf", height = 10, width = 10)  
  keeps_delta[[2]] %>%
    mutate(Date = Date[-1]) %>%
    dplyr::select(-c(Var2))%>%
    melt(id.vars = c("Date")) %>%
    mutate(Jump_Type = factor(variable, levels = c("prop0", "prop1", "prop2", "prop3"), 
                              labels = c("BTC only", "S&P only", "Both", "None"))) %>%
    ggplot() +
    geom_line(aes(x = Date, y = value, linetype = Jump_Type)) +
    facet_wrap(~Jump_Type) +
    theme_bw() +
    labs(y = "Posterior Probability")+
    theme(legend.position = "none")
  ggsave("jump_type_SVIND.pdf", height = 10, width = 10)  
  
  
  
  
  #FIGURE XXX STOCHASTIC VOLATILITY--------
  keeps_j1_long <- keeps_j1 %>%as.data.frame()%>%
    mutate(model = model) %>%
    melt(id.vars = c("model")) %>%
    mutate(Date = rep(Date[-1], each = 4), #change when add more models
           series = "BTC")
  
  keeps_j2_long <- keeps_j2 %>%as.data.frame()%>%
    mutate(model = model) %>%
    melt(id.vars = c("model")) %>%
    mutate(Date = rep(Date[-1], each = 4), #change when add more models
           series = "S&P")
  
  
  rbind(keeps_j1_long,keeps_j2_long) %>%
    mutate(Model = factor(model, levels = c("MALD", "IND", "MVN", "LD"), labels = c("SV-MALD", "SV-IND", "SV-MVN", "SV-LD"))) %>%
    ggplot() +
    geom_line(aes(x = Date, y = value, linetype = Model)) +
    facet_grid(series~Model, scales = "free_y") +
    theme_bw() +
    scale_colour_grey() +
    labs(x = "Date", y = "Jump size")
  ggsave("jump_sizes.pdf", height = 8, width = 14) 
  
  
  
  # POSTERIOR DISTRIBUTIN PLOTS
  
  keeps <- readRDS(paste0("keeps_long/keepsBTCSP" ,"_","IND", ".rds"))
  keeps$rho %>%
    as.data.frame() %>%
    melt() %>%
    #mutate(variable = factor(variable, levels = c("V1", "V2"),labels = c("BTC", "S&P")))%>%
    ggplot() +
    geom_histogram(aes(x = value)) +
    facet_grid(~variable) +
    ggtitle("rho")
  
  apply(keeps$rho, 2, function(x){mean(x >0)})
  
  
  # qqplots-----------
  get_qq <- function(keeps){## Function to make the QQ Plots
    delta <- rep(0,T)
    V <- array(0, dim = c(T+1, 2))
    J <- y <- x <- array(0, dim = c(T, 2))
    V[1,] <- apply(keeps$v[,1,], 2, mean)
    
    sim <- 0
    for (t in 1:T){
      print(t)
      set.seed(t + 4922035+ sim)
      delta <- sample(c(0:3),prob=apply(keeps$lambda, 2, mean), 1)
      
      set.seed(15866245 + t + sim)
      B <- rexp(1)
      Sigma <- matrix(c(mean(keeps$sigma_c[,1])^2,
                        mean(keeps$rhoc)*mean(keeps$sigma_c[,1])*mean(keeps$sigma_c[,2]),
                        mean(keeps$rhoc)*mean(keeps$sigma_c[,1])*mean(keeps$sigma_c[,2]),
                        mean(keeps$sigma_c[,2]^2)),
                      nrow=2)
      xi_c <- apply(keeps$xi_cw, 2, mean)*B+
        sqrt(B)*rtmvnorm(n = 1, mean = c(0,0), sigma = Sigma)
      
      B <- rexp(1)
      xi_y1 <- mean(keeps$xi_y1w)*B + sqrt(B)*rnorm(1,0,mean(keeps$xi_y1eta)) #SHOULD THIS BE SQRT(ETA)?
      B <- rexp(1)
      xi_y2 <- mean(keeps$xi_y2w)*B + sqrt(B)*rnorm(1,0,mean(keeps$xi_y2eta)) #SHOULD THIS BE SQRT(ETA)?
      
      
      J = xi_c*(delta==2) +cbind(xi_y1,0)*(delta==0) + cbind(0,xi_y2)*(delta==1)
      
      Sigma <- matrix(c(V[t,1],
                        mean(keeps$rho[,1])*sqrt(prod(V[t,])),
                        mean(keeps$rho[,3])*mean(keeps$sigma_v[,1])*V[t,1],0,
                        mean(keeps$rho[,1])*sqrt(prod(V[t,])),V[t,2],0,mean(keeps$rho[,4])*mean(keeps$sigma_v[,2])*V[t,2],
                        mean(keeps$rho[,3])*mean(keeps$sigma_v[,1])*V[t,1],0,mean(keeps$sigma_v[,1])^2*V[t,1],mean(keeps$rho[,2])*prod(apply(keeps$sigma_v, 2, mean))*sqrt(prod(V[t,])),
                        0,mean(keeps$rho[,4])*mean(keeps$sigma_v[,2])*V[t,2],mean(keeps$rho[,2])*prod(apply(keeps$sigma_v, 2, mean))*sqrt(prod(V[t,])),mean(keeps$sigma_v[,2])^2*V[t,2]),nrow=4)
      
      set.seed(463468+t)
      temp <- rtmvnorm(n = 1,
                       mean = c(apply(keeps$mu, 2 ,mean) + J,
                                apply(keeps$theta,2, mean) + apply(keeps$phi, 2, mean)*(V[t,] - apply(keeps$theta, 2, mean))),
                       sigma = Sigma, lower=c(-Inf,-Inf, 0, 0))
      
      V[t+1,] <- temp[c(3:4)]
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
      xlab("Actual Quantiles") + ylab("Simulated Quantiles") + theme_bw() + ggtitle("BTC")
    
    p2 <- ggplot() +
      geom_point(aes(x=quantile(QQdat$SP,seq(0.01,0.99,0.01)),y=quantile(QQdat$X2,seq(0.01,0.99,0.01)))) +
      geom_abline(slope=1,intercept=0) +
      #xlim(c(-15,15)) + ylim(c(-15,15)) +
      xlab("Actual Quantiles") + ylab("Simulated Quantiles") + theme_bw() + ggtitle("S&P")
    p3 <- grid.arrange(p1,p2, nrow = 1)
    return(p3)
  }
  
  model <- "MALD"
  keeps <- readRDS(paste0("keeps_long/keepsBTCSP_",model , ".rds"))
  plot <- get_qq(keeps)
  ggsave(paste0("QQ_", model, ".pdf"),plot, width = 12, height = 6)
  
  model <- "IND"
  keeps <- readRDS(paste0("keeps_long/keepsBTCSP_",model , ".rds"))
  plot <- get_qq(keeps)
  ggsave(paste0("QQ_", model, ".pdf"),plot)
  
  model <- "MVN"
  keeps <- readRDS(paste0("keeps_long/keepsBTCSP_",model , ".rds"))
  plot <- get_qq(keeps)
  ggsave(paste0("QQ_", model, ".pdf"),plot, width = 12, height = 6)
  
  model <- "LD"
  keeps <- readRDS(paste0("keeps_long/keepsBTCSP_",model , ".rds"))
  plot <- get_qq(keeps)
  ggsave(paste0("QQ_", model, ".pdf"),plot)
  

