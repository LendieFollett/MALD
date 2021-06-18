rm(list = ls()) 
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
library(xtable)
library(tibble)#rownames_to_column()
library(tidyr)

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

docreds<- function(x, total,q){
  R <- total
  if(length(dim(x)) == 2){ #if it's a data frame
    return(apply(x[1:R,], 2, quantile, q))
  #}#else if (length(dim(x)) > 2){
    #return(apply(x[1:R,,], 2:3, function(x){(median(x))}))
  }else{
    return(quantile(x[1:R],q))
  }
}

############################################################
#----SHORT TIME SERIES ONLY---------------------------------
#ALL CRYPTO, MEME STOCKS
############################################################
#load data
getSymbols("BTC-USD",from = "2020-12-01",to = "2021-05-31")
BTC <- as.data.frame(`BTC-USD`)
BTC$Date <- seq(as.Date("2020-12-01"),as.Date("2021-05-31"),by="days")
BTC$`BTC-USD.Close`[BTC$Date=="2020-04-17"] <- 7096.18
getSymbols("GME",from = "2020-12-01",to = "2021-05-31")
GME <- as.data.frame(GME)
GME$Date <- as.Date(rownames(GME))
getSymbols("AMC",from = "2020-12-01",to = "2021-05-31")
AMC <- as.data.frame(AMC)
AMC$Date <- as.Date(rownames(AMC))
getSymbols("DOGE-USD",from = "2020-12-01",to = "2021-05-31")
DOGE <- as.data.frame(`DOGE-USD`)
DOGE$Date <- as.Date(rownames(DOGE))
getSymbols("^GSPC",from = "2020-12-01",to = "2021-05-31")
SP500 <- as.data.frame(`GSPC`)
SP500$Date <- as.Date(rownames(SP500))

S <- BTC %>%merge(GME) %>% merge(AMC) %>% merge(DOGE) %>% merge(SP500)
T <- nrow(S) - 1
Date <- S$Date[-1]



keeps <- readRDS(paste0("keeps_short/keeps_","SVMALD" ,"_","BTC", ".rds"))
names <- lapply(keeps[c(4,6:17)], domean, total = 10) %>%unlist %>%names

keeps_summary <- array(dim = c(16, length(names)))
keeps_v1 <- array(dim = c(16,dim(keeps$v)[2]))
keeps_v2 <- array(dim = c(16,dim(keeps$v)[2]))
model <- rep(NA, 16)
data <- rep(NA, 16)
colnames(keeps_summary) <- names
keeps_creds <- list()
j = 0
for (m in c("SVMALD", "SVMVN", "SVLD", "SVIND")){
  for (i in c("BTC", "DOGE", "AMC", "GMC")){
    j = j + 1
    keeps <- readRDS(paste0("keeps_short/keeps_",m ,"_",i, ".rds"))
    keeps_summary[j,] <- lapply(keeps[c(4,6:17)], domean, total = 20000) %>%unlist
    keeps_creds[[j]] <- rbind(lapply(keeps[c(4,6:17)], quantile,.025, total = 20000) %>%unlist,
                             lapply(keeps[c(4,6:17)], quantile,.975, total = 20000) %>%unlist)%>% t() %>%
      as.data.frame %>%
      rownames_to_column(var = "parameter")%>%
      mutate_if(is.numeric, round, 3) 
    keeps_v1[j,] <- apply(keeps$v[1:20000,,1], 2, function(x){(mean(x))}) #alternative sv
    keeps_v2[j,] <- apply(keeps$v[1:20000,,2], 2, function(x){(mean(x))}) #sp sv
    model[j] <-m
    data[j] <- i
  }
}


keeps_summary <- keeps_summary %>%as.data.frame()%>%
  mutate_all(round, digits = 2) %>%
  #mutate(lambda = paste0("(",lambda1,", ", lambda2,", ", lambda3,", ", lambda4,")"))%>%
  #dplyr::select(-c(lambda1, lambda2, lambda3, lambda4)) %>%
  mutate(model = model) %>%t()


#TABLE XXX POSTERIOR MEANS OF PARAMETERS------
keeps_summary%>%
  xtable() %>%
  print()

#FIGURE XXX STOCHASTIC VOLATILITY--------
#ALTERNATIVE CURRENCY
keeps_v1_long <- keeps_v1 %>%as.data.frame()%>%
  mutate(model = model,
         series = date) %>%
  melt(id.vars = c("model", "series")) %>%
  mutate(Date = rep(Date, each = 16)) #CHECK THIS STRUCTURE 
#S&P
keeps_v2_long <- keeps_v2 %>%as.data.frame()%>%
  mutate(model = model,
         series = data) %>%
  melt(id.vars = c("model", "series")) %>%
  mutate(Date = rep(Date, each = 16), #CHECK THIS STRUCTURE 
         series = "S&P")

#not sure the best way to display this, may need to modify
rbind(keeps_v1_long,keeps_v2_long) %>%
  ggplot() +
  geom_line(aes(x = Date, y = value, linetype = model)) +
  facet_grid(series~model., scales = "free_y") + 
  theme_bw() +
  scale_colour_grey()


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

keeps_summary0 <- array(dim = c(4, length(names)))
keeps_v1 <- array(dim = c(4,dim(keeps$v)[2]))
keeps_v2 <- array(dim = c(4,dim(keeps$v)[2]))
model <- rep(NA, 4)
keeps_delta <- list()
keeps_creds <- list()
colnames(keeps_summary0) <- names

j = 0
for (m in c("MALD", "IND", "MVN", "LD")){# "SVMVN", "SVLD",
    j = j + 1
    keeps <- readRDS(paste0("keeps_long/keepsBTCSP_",m , ".rds"))
    keeps_summary0[j,] <- lapply(keeps[c(4,6:17)], domean, total = 20000) %>%unlist
    keeps_v1[j,] <- apply(keeps$v[1:20000,,1], 2, function(x){(mean(x))}) #alternative sv
    keeps_v2[j,] <- apply(keeps$v[1:20000,,2], 2, function(x){(mean(x))}) #sp sv
    keeps_delta[[j]] <- keeps$delta %>%melt %>%group_by(Var2) %>%
      summarise(prop0 = mean(value == 0),prop1 = mean(value == 1),prop2 = mean(value == 2),prop3 = mean(value == 3)) %>%
      ungroup()
    keeps_creds[[j]] <- rbind(lapply(keeps[c(4,6:17)], docreds,q=.025, total = 20000) %>%unlist,
                              lapply(keeps[c(4,6:17)], docreds,q=.975, total = 20000) %>%unlist)%>% t() %>%
      as.data.frame %>%
      rownames_to_column(var = "parameter")%>%
      mutate_if(is.numeric, round, 2)  %>%mutate(parameter = gsub(".2.5%", "", parameter))
    model[j] <-m
}


keeps_ci <- do.call(rbind,keeps_creds) %>% 
  mutate(model = rep(model, each = length(names)), #tack on model
         ci = paste0("(", V1,", " ,V2, ")")) %>% #format to (l, u) for latex table
  dplyr::select(-c(V1, V2)) %>% spread(model, ci)%>%#make wide by model
  mutate(parameter = paste0("$\\", parameter, "$")) #math formatting start

keeps_summary <- keeps_summary0 %>%as.data.frame()%>%
  mutate_all(round, digits = 3) %>%
  #mutate(lambda = paste0("(",lambda1,", ", lambda2,", ", lambda3,", ", lambda4,")"))%>%
  #dplyr::select(-c(lambda1, lambda2, lambda3, lambda4)) %>%
  mutate(model = model) %>%t() %>% as.data.frame() %>%
  rownames_to_column(var = "parameter") %>%
  mutate(parameter = paste0("$\\", parameter, "$")) %>%
  rename(MALD =V1, IND = V2 , MVN = V3, LD = V4) %>% rbind(keeps_ci) %>%
  arrange(parameter)


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
  
  
p1 <-   keeps_v1_long %>% subset(model == "MALD") %>%
    merge(data.frame(BTC=100*(log(S[-1,c( "BTC-USD.Close") ])-
                       log(S[-nrow(S),c( "BTC-USD.Close") ])),Date=Date[-1]) , by = "Date") %>%
    ggplot() +
    geom_line(aes(x = Date, y = BTC), alpha= I(.5)) +
    geom_line(aes(x = Date, y = value^.5)) +
    theme_bw() +labs(x = "")

p2 <- keeps_v2_long %>% subset(model == "MALD") %>%
  merge(data.frame(SP=100*(log(S[-1,c( "GSPC.Close") ])-
                              log(S[-nrow(S),c( "GSPC.Close") ])),Date=Date[-1]) , by = "Date") %>%
  ggplot() +
  geom_line(aes(x = Date, y = SP), alpha= I(.5)) +
  geom_line(aes(x = Date, y = value^.5)) +
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
  
