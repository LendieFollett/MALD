library(ggplot2)
library(grid)
library(gridExtra)
library(MASS)
library(tmvtnorm)
library(dplyr)
library(quantmod)

####### Get Data ########
getSymbols("ETH-USD",from = "2016-01-01",to = "2021-01-31")
ETH <- as.data.frame(`ETH-USD`)
ETH$Date <- as.Date(rownames(ETH))
ETH$`ETH-USD.Close`[ETH$Date=="2020-10-13"] <- 381.19
ETH$`ETH-USD.Close`[ETH$Date=="2020-10-12"] <- 387.73
ETH$`ETH-USD.Close`[ETH$Date=="2020-10-09"] <- 356.59
ETH$`ETH-USD.Close`[ETH$Date=="2020-04-17"] <- 171.64

getSymbols("XRP-USD",from = "2016-01-01",to = "2021-01-31")
XRP <- as.data.frame(`XRP-USD`)
XRP$Date <- as.Date(rownames(XRP))
XRP$`XRP-USD.Close`[XRP$Date=="2020-10-13"] <- 0.2563
XRP$`XRP-USD.Close`[XRP$Date=="2020-10-12"] <- 0.2564
XRP$`XRP-USD.Close`[XRP$Date=="2020-10-09"] <- 0.2543
XRP$`XRP-USD.Close`[XRP$Date=="2020-04-17"] <- 0.1902

getSymbols("BTC-USD",from = "2016-01-01",to = "2021-01-31")
BTC <- as.data.frame(`BTC-USD`)
BTC$Date <- as.Date(rownames(BTC))
BTC$`BTC-USD.Close`[BTC$Date=="2020-10-13"] <- 11425.90
BTC$`BTC-USD.Close`[BTC$Date=="2020-10-12"] <- 11555.36
BTC$`BTC-USD.Close`[BTC$Date=="2020-10-09"] <- 11064.46
BTC$`BTC-USD.Close`[BTC$Date=="2020-04-17"] <- 7096.18

getSymbols("GILD",from = "2016-01-01",to = "2021-01-31")
GILD <- as.data.frame(GILD)
GILD$Date <- as.Date(rownames(GILD))

getSymbols("VRTX",from = "2016-01-01",to = "2021-01-31")
VRTX <- as.data.frame(VRTX)
VRTX$Date <- as.Date(rownames(VRTX))

R <- 10000
n_chns <- 4

S <- right_join(ETH,XRP) %>%
  right_join(BTC) %>%
  right_join(GILD) %>%
  right_join(VRTX)
T <- nrow(S) - 1
S <- S[,c("Date",
                    "BTC-USD.Close",
                    "ETH-USD.Close",
                    "XRP-USD.Close",
                    "GILD.Close",
                    "VRTX.Close")]


names <- c("BTC")

for (i in 1:length(names)){
  #keeps <- readRDS(paste0("keeps/keeps",names[i],".rds"))
  dt <- data.frame(iter = rep(1:R,n_chns),
                   chain = factor(keeps$chain),
                   mu = keeps$mu,
                   theta = keeps$theta,
                   kappa = 1-keeps$phi,
                   sigma_v = keeps$sigma_v,
                   rho = keeps$rho,
                   mu_y = keeps$xi_yw,
                   nu_y = keeps$xi_yeta,
                   lambda = keeps$lambda)
  ### mu plot
  g1 <- ggplot(data=dt) +
    geom_histogram(aes(x=252*mu)) +
    geom_vline(aes(xintercept=252*mean(mu)),color="red",lwd=1.5) +
    geom_vline(aes(xintercept=252*quantile(mu,0.025)),linetype="dashed",color="red",lwd=1.5) +
    geom_vline(aes(xintercept=252*quantile(mu,0.975)),linetype="dashed",color="red",lwd=1.5) +
    xlab("") + ylab("") + theme_bw()
  
  g2 <- ggplot(data=dt) +
    geom_line(aes(x=iter,y=252*mu,color=chain)) +
    geom_hline(aes(yintercept=252*mean(mu)),color="red",lwd=1.5) +
    geom_hline(aes(yintercept=252*quantile(mu,0.025)),linetype="dashed",color="red",lwd=1.5) +
    geom_hline(aes(yintercept=252*quantile(mu,0.975)),linetype="dashed",color="red",lwd=1.5) +
    xlab("") + ylab("") + theme_bw() + theme(legend.position = "none")
  g <- grid.arrange(g1,g2,nrow=2)
  ggsave(paste0("Plots",names[i],"/mu.pdf"),g)
  
  ### theta plot
  g1 <- ggplot(data=dt) +
    geom_histogram(aes(x=sqrt(252*theta))) +
    geom_vline(aes(xintercept=mean(sqrt(252*theta))),color="red",lwd=1.5) +
    geom_vline(aes(xintercept=quantile(sqrt(252*theta),0.025)),linetype="dashed",color="red",lwd=1.5) +
    geom_vline(aes(xintercept=quantile(sqrt(252*theta),0.975)),linetype="dashed",color="red",lwd=1.5) +
    xlab("") + ylab("") + theme_bw()
  
  g2 <- ggplot(data=dt) +
    geom_line(aes(x=iter,y=sqrt(252*theta),color=chain)) +
    geom_hline(aes(yintercept=mean(sqrt(252*theta))),color="red",lwd=1.5) +
    geom_hline(aes(yintercept=quantile(sqrt(252*theta),0.025)),linetype="dashed",color="red",lwd=1.5) +
    geom_hline(aes(yintercept=quantile(sqrt(252*theta),0.975)),linetype="dashed",color="red",lwd=1.5) +
    xlab("") + ylab("") + theme_bw() + theme(legend.position = "none")
  g <- grid.arrange(g1,g2,nrow=2)
  ggsave(paste0("Plots",names[i],"/theta.pdf"),g)
  
  ### kappa plot
  g1 <- ggplot(data=dt) +
    geom_histogram(aes(x=kappa)) +
    geom_vline(aes(xintercept=mean(kappa)),color="red",lwd=1.5) +
    geom_vline(aes(xintercept=quantile(kappa,0.025)),linetype="dashed",color="red",lwd=1.5) +
    geom_vline(aes(xintercept=quantile(kappa,0.975)),linetype="dashed",color="red",lwd=1.5) +
    xlab("") + ylab("") + theme_bw()
  
  g2 <- ggplot(data=dt) +
    geom_line(aes(x=iter,y=kappa,color=chain)) +
    geom_hline(aes(yintercept=mean(kappa)),color="red",lwd=1.5) +
    geom_hline(aes(yintercept=quantile(kappa,0.025)),linetype="dashed",color="red",lwd=1.5) +
    geom_hline(aes(yintercept=quantile(kappa,0.975)),linetype="dashed",color="red",lwd=1.5) +
    xlab("") + ylab("") + theme_bw() + theme(legend.position = "none")
  g <- grid.arrange(g1,g2,nrow=2)
  ggsave(paste0("Plots",names[i],"/kappa.pdf"),g)
  
  ### sigma_v plot
  g1 <- ggplot(data=dt) +
    geom_histogram(aes(x=sqrt(252*sigma_v))) +
    geom_vline(aes(xintercept=mean(sqrt(252*sigma_v))),color="red",lwd=1.5) +
    geom_vline(aes(xintercept=quantile(sqrt(252*sigma_v),0.025)),linetype="dashed",color="red",lwd=1.5) +
    geom_vline(aes(xintercept=quantile(sqrt(252*sigma_v),0.975)),linetype="dashed",color="red",lwd=1.5) +
    xlab("") + ylab("") + theme_bw()
  
  g2 <- ggplot(data=dt) +
    geom_line(aes(x=iter,y=sqrt(252*sigma_v),color=chain)) +
    geom_hline(aes(yintercept=mean(sqrt(252*sigma_v))),color="red",lwd=1.5) +
    geom_hline(aes(yintercept=quantile(sqrt(252*sigma_v),0.025)),linetype="dashed",color="red",lwd=1.5) +
    geom_hline(aes(yintercept=quantile(sqrt(252*sigma_v),0.975)),linetype="dashed",color="red",lwd=1.5) +
    xlab("") + ylab("") + theme_bw() + theme(legend.position = "none")
  g <- grid.arrange(g1,g2,nrow=2)
  ggsave(paste0("Plots",names[i],"/sigma_v.pdf"),g)
  
  ### rho plot
  g1 <- ggplot(data=dt) +
    geom_histogram(aes(x=rho)) +
    geom_vline(aes(xintercept=mean(rho)),color="red",lwd=1.5) +
    geom_vline(aes(xintercept=quantile(rho,0.025)),linetype="dashed",color="red",lwd=1.5) +
    geom_vline(aes(xintercept=quantile(rho,0.975)),linetype="dashed",color="red",lwd=1.5) +
    xlab("") + ylab("") + theme_bw()
  
  g2 <- ggplot(data=dt) +
    geom_line(aes(x=iter,y=rho,color=chain)) +
    geom_hline(aes(yintercept=mean(rho)),color="red",lwd=1.5) +
    geom_hline(aes(yintercept=quantile(rho,0.025)),linetype="dashed",color="red",lwd=1.5) +
    geom_hline(aes(yintercept=quantile(rho,0.975)),linetype="dashed",color="red",lwd=1.5) +
    xlab("") + ylab("") + theme_bw() + theme(legend.position = "none")
  g <- grid.arrange(g1,g2,nrow=2)
  ggsave(paste0("Plots",names[i],"/rho.pdf"),g)
  
  ### mu_y plot
  g1 <- ggplot(data=dt) +
    geom_histogram(aes(x=mu_y)) +
    geom_vline(aes(xintercept=mean(mu_y)),color="red",lwd=1.5) +
    geom_vline(aes(xintercept=quantile(mu_y,0.025)),linetype="dashed",color="red",lwd=1.5) +
    geom_vline(aes(xintercept=quantile(mu_y,0.975)),linetype="dashed",color="red",lwd=1.5) +
    xlab("") + ylab("") + theme_bw()
  
  g2 <- ggplot(data=dt) +
    geom_line(aes(x=iter,y=mu_y,color=chain)) +
    geom_hline(aes(yintercept=mean(mu_y)),color="red",lwd=1.5) +
    geom_hline(aes(yintercept=quantile(mu_y,0.025)),linetype="dashed",color="red",lwd=1.5) +
    geom_hline(aes(yintercept=quantile(mu_y,0.975)),linetype="dashed",color="red",lwd=1.5) +
    xlab("") + ylab("") + theme_bw() + theme(legend.position = "none")
  g <- grid.arrange(g1,g2,nrow=2)
  ggsave(paste0("Plots",names[i],"/mu_y.pdf"),g)
  
  ### nu_y plot
  g1 <- ggplot(data=dt) +
    geom_histogram(aes(x=nu_y)) +
    geom_vline(aes(xintercept=mean(nu_y)),color="red",lwd=1.5) +
    geom_vline(aes(xintercept=quantile(nu_y,0.025)),linetype="dashed",color="red",lwd=1.5) +
    geom_vline(aes(xintercept=quantile(nu_y,0.975)),linetype="dashed",color="red",lwd=1.5) +
    xlab("") + ylab("") + theme_bw()
  
  g2 <- ggplot(data=dt) +
    geom_line(aes(x=iter,y=nu_y,color=chain)) +
    geom_hline(aes(yintercept=mean(nu_y)),color="red",lwd=1.5) +
    geom_hline(aes(yintercept=quantile(nu_y,0.025)),linetype="dashed",color="red",lwd=1.5) +
    geom_hline(aes(yintercept=quantile(nu_y,0.975)),linetype="dashed",color="red",lwd=1.5) +
    xlab("") + ylab("") + theme_bw() + theme(legend.position = "none")
  g <- grid.arrange(g1,g2,nrow=2)
  ggsave(paste0("Plots",names[i],"/nu_y.pdf"),g)
  
  ### lambda plot
  g1 <- ggplot(data=dt) +
    geom_histogram(aes(x=lambda*252)) +
    geom_vline(aes(xintercept=mean(lambda*252)),color="red",lwd=1.5) +
    geom_vline(aes(xintercept=quantile(lambda*252,0.025)),linetype="dashed",color="red",lwd=1.5) +
    geom_vline(aes(xintercept=quantile(lambda*252,0.975)),linetype="dashed",color="red",lwd=1.5) +
    xlab("") + ylab("") + theme_bw()
  
  g2 <- ggplot(data=dt) +
    geom_line(aes(x=iter,y=lambda*252,color=chain)) +
    geom_hline(aes(yintercept=mean(lambda*252)),color="red",lwd=1.5) +
    geom_hline(aes(yintercept=quantile(lambda*252,0.025)),linetype="dashed",color="red",lwd=1.5) +
    geom_hline(aes(yintercept=quantile(lambda*252,0.975)),linetype="dashed",color="red",lwd=1.5) +
    xlab("") + ylab("") + theme_bw() + theme(legend.position = "none")
  g <- grid.arrange(g1,g2,nrow=2)
  ggsave(paste0("Plots",names[i],"/lambda.pdf"),g)
  
  Vol = data.frame(date = S$Date,
                   S = S[,i+1],
                   vol = sqrt(apply(keeps$v,2,mean)*252))
  
  if (i==1) {
    g1 <- ggplot() + 
    geom_rect(aes(xmin=as.Date("2017-07-01"),
                  xmax=as.Date("2018-03-31"),
                  ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue") +
    geom_rect(aes(xmin=as.Date("2020-02-01"),
                  xmax=as.Date("2021-01-31"),
                  ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue") +
      geom_line(data = Vol,aes(x=date,y=S)) + 
      xlab("") + ylab("Price") + theme_bw()
    
    g2 <- ggplot() + 
      geom_rect(aes(xmin=as.Date("2017-07-01"),
                    xmax=as.Date("2018-03-31"),
                    ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue") +
      geom_rect(aes(xmin=as.Date("2020-02-01"),
                    xmax=as.Date("2021-01-31"),
                    ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue") +
      geom_line(data = Vol,aes(x=date,y=vol)) + 
      xlab("Date") + ylab("Estimated Volatility") + theme_bw()
  } else if (i==2){
    g1 <- ggplot() + 
    geom_rect(aes(xmin=as.Date("2016-01-01"),
                  xmax=as.Date("2016-08-31"),
                  ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue") +
      geom_rect(aes(xmin=as.Date("2017-03-01"),
                    xmax=as.Date("2018-01-31"),
                    ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue") +
      geom_rect(aes(xmin=as.Date("2020-02-01"),
                    xmax=as.Date("2021-01-31"),
                    ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue")+
      geom_line(data = Vol,aes(x=date,y=S)) + 
      xlab("") + ylab("Price") + theme_bw()
    
    g2 <- ggplot() + 
      geom_rect(aes(xmin=as.Date("2016-01-01"),
                    xmax=as.Date("2016-08-31"),
                    ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue") +
      geom_rect(aes(xmin=as.Date("2017-03-01"),
                    xmax=as.Date("2018-01-31"),
                    ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue") +
      geom_rect(aes(xmin=as.Date("2020-02-01"),
                    xmax=as.Date("2021-01-31"),
                    ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue")+
      geom_line(data = Vol,aes(x=date,y=vol)) + 
      xlab("Date") + ylab("Estimated Volatility") + theme_bw()
  } else if (i==3){
    g1 <- ggplot() + 
    geom_rect(aes(xmin=as.Date("2017-03-01"),
                  xmax=as.Date("2018-04-30"),
                  ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue") +
      geom_rect(aes(xmin=as.Date("2020-09-01"),
                    xmax=as.Date("2021-01-31"),
                    ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue")+
      geom_line(data = Vol,aes(x=date,y=S)) + 
      xlab("") + ylab("Price") + theme_bw()
    
    g2 <- ggplot() + 
      geom_rect(aes(xmin=as.Date("2017-03-01"),
                    xmax=as.Date("2018-04-30"),
                    ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue") +
      geom_rect(aes(xmin=as.Date("2020-09-01"),
                    xmax=as.Date("2021-01-31"),
                    ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue")+
      geom_line(data = Vol,aes(x=date,y=vol)) + 
      xlab("Date") + ylab("Estimated Volatility") + theme_bw()
  } else if (i==4) {
    g1 <- ggplot() + 
    geom_rect(aes(xmin=as.Date("2020-02-01"),
                  xmax=as.Date("2021-01-31"),
                  ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue")+
      geom_line(data = Vol,aes(x=date,y=S)) + 
      xlab("") + ylab("Price") + theme_bw()
    
    g2 <- ggplot() + 
      geom_rect(aes(xmin=as.Date("2020-02-01"),
                    xmax=as.Date("2021-01-31"),
                    ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue")+
      geom_line(data = Vol,aes(x=date,y=vol)) + 
      xlab("Date") + ylab("Estimated Volatility") + theme_bw()
  } else if (i==5){
    g1 <- ggplot() + 
    geom_rect(aes(xmin=as.Date("2016-01-01"),
                  xmax=as.Date("2016-12-31"),
                  ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue") +
      geom_rect(aes(xmin=as.Date("2020-02-01"),
                    xmax=as.Date("2021-01-31"),
                    ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue")+
      geom_line(data = Vol,aes(x=date,y=S)) + 
      xlab("") + ylab("Price") + theme_bw()
    
    g2 <- ggplot() + 
      geom_rect(aes(xmin=as.Date("2016-01-01"),
                    xmax=as.Date("2016-12-31"),
                    ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue") +
      geom_rect(aes(xmin=as.Date("2020-02-01"),
                    xmax=as.Date("2021-01-31"),
                    ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue")+
      geom_line(data = Vol,aes(x=date,y=vol)) + 
      xlab("Date") + ylab("Estimated Volatility") + theme_bw()
  }
  
  G1 <- ggplotGrob(g1)
  G2 <- ggplotGrob(g2)
  ggsave(paste0("Plots",names[i],"/Volatility.pdf"),grid.draw(rbind(G1,G2)))

  Jumps = data.frame(date = S$Date[-1],
                     S = log(S[-1,i+1]) - log(S[-(T+1),i+1]),
                     J = apply(keeps$J,2,mean))
  
  if (i==1) {
    g1 <- ggplot() + 
      geom_rect(aes(xmin=as.Date("2017-07-01"),
                  xmax=as.Date("2018-03-31"),
                  ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue") +
      geom_rect(aes(xmin=as.Date("2020-02-01"),
                    xmax=as.Date("2021-01-31"),
                    ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue") +
      ylim(-0.5,0.25) + geom_line(data = Jumps,aes(x=date,y=S)) + theme_bw()  +
      xlab("") + ylab("log-Returns")
    
    g2 <- ggplot() + 
      geom_rect(aes(xmin=as.Date("2017-07-01"),
                    xmax=as.Date("2018-03-31"),
                    ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue") +
      geom_rect(aes(xmin=as.Date("2020-02-01"),
                    xmax=as.Date("2021-01-31"),
                    ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue") +
      ylim(-0.5,0.25) + geom_line(data = Jumps,aes(x=date,y=J/100)) + theme_bw()  +
      xlab("Date") + ylab("Estimated Jumps")
  } else if (i==2){
    g1 <- ggplot() + 
    geom_rect(aes(xmin=as.Date("2016-01-01"),
                  xmax=as.Date("2016-08-31"),
                  ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue") +
      geom_rect(aes(xmin=as.Date("2017-03-01"),
                    xmax=as.Date("2018-01-31"),
                    ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue") +
      geom_rect(aes(xmin=as.Date("2020-02-01"),
                    xmax=as.Date("2021-01-31"),
                    ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue") +
      ylim(-0.6,0.6)+ geom_line(data = Jumps,aes(x=date,y=S)) + theme_bw()  +
      xlab("") + ylab("log-Returns")
    
    g2 <- ggplot() + 
      geom_rect(aes(xmin=as.Date("2016-01-01"),
                    xmax=as.Date("2016-08-31"),
                    ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue") +
      geom_rect(aes(xmin=as.Date("2017-03-01"),
                    xmax=as.Date("2018-01-31"),
                    ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue") +
      geom_rect(aes(xmin=as.Date("2020-02-01"),
                    xmax=as.Date("2021-01-31"),
                    ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue") +
      ylim(-0.6,0.6)+ geom_line(data = Jumps,aes(x=date,y=J/100)) + theme_bw()  +
      xlab("Date") + ylab("Estimated Jumps")
  } else if (i==3){
    g1 <- ggplot() + 
    geom_rect(aes(xmin=as.Date("2017-03-01"),
                  xmax=as.Date("2018-04-30"),
                  ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue") +
      geom_rect(aes(xmin=as.Date("2017-03-01"),
                    xmax=as.Date("2018-01-31"),
                    ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue") +
      geom_rect(aes(xmin=as.Date("2020-09-01"),
                    xmax=as.Date("2021-01-31"),
                    ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue") +
      ylim(-0.6,0.8)+ geom_line(data = Jumps,aes(x=date,y=S)) + theme_bw()  +
      xlab("") + ylab("log-Returns")
    
    g2 <- ggplot() + 
      geom_rect(aes(xmin=as.Date("2017-03-01"),
                    xmax=as.Date("2018-04-30"),
                    ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue") +
      geom_rect(aes(xmin=as.Date("2017-03-01"),
                    xmax=as.Date("2018-01-31"),
                    ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue") +
      geom_rect(aes(xmin=as.Date("2020-09-01"),
                    xmax=as.Date("2021-01-31"),
                    ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue") +
      ylim(-0.6,0.8)+ geom_line(data = Jumps,aes(x=date,y=J/100)) + theme_bw()  +
      xlab("Date") + ylab("Estimated Jumps")
  } else if (i==4) {
    g1 <- ggplot() + 
    geom_rect(aes(xmin=as.Date("2020-02-01"),
                  xmax=as.Date("2021-01-31"),
                  ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue") +
      ylim(-0.1,0.1)+ geom_line(data = Jumps,aes(x=date,y=S)) + theme_bw()  +
      xlab("") + ylab("log-Returns")
    
    g2 <- ggplot() + 
      geom_rect(aes(xmin=as.Date("2020-02-01"),
                    xmax=as.Date("2021-01-31"),
                    ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue") +
      ylim(-0.1,0.1)+ geom_line(data = Jumps,aes(x=date,y=J/100)) + theme_bw()  +
      xlab("Date") + ylab("Estimated Jumps")
  } else if (i==5){
    g1 <- ggplot() + 
    geom_rect(aes(xmin=as.Date("2016-01-01"),
                  xmax=as.Date("2016-12-31"),
                  ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue") +
      geom_rect(aes(xmin=as.Date("2020-02-01"),
                    xmax=as.Date("2021-01-31"),
                    ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue") +
      ylim(-0.3,0.2)+ geom_line(data = Jumps,aes(x=date,y=S)) + theme_bw()  +
      xlab("") + ylab("log-Returns")
    
    g2 <- ggplot() + 
      geom_rect(aes(xmin=as.Date("2016-01-01"),
                    xmax=as.Date("2016-12-31"),
                    ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue") +
      geom_rect(aes(xmin=as.Date("2020-02-01"),
                    xmax=as.Date("2021-01-31"),
                    ymin=-Inf,ymax=Inf),alpha=0.2,fill="blue") +
      ylim(-0.3,0.2)+ geom_line(data = Jumps,aes(x=date,y=J/100)) + theme_bw()  +
      xlab("Date") + ylab("Estimated Jumps")
  }
  
  G1 <- ggplotGrob(g1)
  G2 <- ggplotGrob(g2)
  ggsave(paste0("Plots",names[i],"/Jumps.pdf"),grid.draw(rbind(G1,G2)))

  # delta <- rep(0,T)
  # xiy1 <- xiy1s <- rep(0,T)
  # V <- rep(0,T+1)
  # J <- y <- x <- rep(0,T)
  # V[1] <- mean(keeps$v[,1])
  # 
  # sim <- 0
  # for (t in 1:T){
  #   print(t)
  #   set.seed(t + 4922035+ sim)
  #   whichone <- sample(c(0:1),prob=c(mean(keeps$lambda),1-mean(keeps$lambda)), 1)
  #   delta[t] <- whichone
  #   
  #   set.seed(15866245 + t + sim)
  #   xiy1s[t] <- rexp(1)
  #   xiy1[t] <- (mean(keeps$xi_yw)*xiy1s[t]+ 
  #                 sqrt(xiy1s[t])*mean(keeps$xi_yeta)*rnorm(1,0,1))
  #   
  #   J[t] = (delta[t] == 0)*xiy1[t]
  #   
  #   Sigma <- matrix(c(V[t],mean(keeps$rho)*mean(keeps$sigma_v)*V[t],mean(keeps$rho)*mean(keeps$sigma_v)*V[t],mean(keeps$sigma_v)^2*V[t]),
  #                   nrow=2)
  #   idx <- 0
  #   set.seed(463468+i + idx + sim)
  #   temp <- rtmvnorm(n = 1, 
  #                    mean = c(x[t] + mean(keeps$mu) + J[t],
  #                             mean(keeps$theta) + mean(keeps$phi)*(V[t] - mean(keeps$theta))),
  #                    sigma = Sigma, lower=c(-Inf,0))
  #   
  #   V[t+1] <- temp[2]
  #   y[t] <- temp[1]
  #   if( t+1 <= T){ x[t+1] <- 0 }  
  # }
  # 
  # QQdat = data.frame(S1 = 100*(log(S[-1,i+1]) - log(S[-(T+1),i+1])),
  #                    y1 = y)
  # 
  # g3 <- ggplot() +
  #   geom_point(aes(x=quantile(QQdat$S1,seq(0.01,0.99,0.01)),y=quantile(QQdat$y1,seq(0.01,0.99,0.01)))) +
  #   geom_abline(slope=1,intercept=0) +
  #   #xlim(c(-15,15)) + ylim(c(-15,15)) +
  #   xlab("Actual Quantiles") + ylab("Simulated Quantiles") + theme_bw()
  # ggsave(paste0("Plots",names[i],"/QQ.pdf"),g3)
}


