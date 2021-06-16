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



# MEANS  ----------


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


keeps <- readRDS(paste0("keeps_",SVMALD ,"_",BTC, ".rds"))
names <- lapply(keeps[c(4,6:17)], domean, total = 10) %>%unlist %>%names

keeps_summary <- array(dim = c(16, length(names)))
keeps_v <- array(dim = c(16,dim(keeps$v)[2]))
model <- rep(NA, 16)
data <- rep(NA, 16)
colnames(keeps_summary) <- names

j = 0
for (m in c("SVMALD", "SVMVN", "SVLD", "SVIND")){
  for (i in c("BTC", "DOGE", "AMC", "GMC")){
    j = j + 1
    keeps <- readRDS(paste0("keeps_",m ,"_",i, ".rds"))
    keeps_summary[j,] <- lapply(keeps[c(4,6:17)], domean, total = 10) %>%unlist
    keeps_v[j,] <- apply(keeps$v[1:R,,], 2:3, function(x){(mean(x))})
    model[j] <-m
    data[j] <- i
  }
}

keeps_summary <- keeps_summary %>%as.data.frame()%>%
  mutate_all(round, digits = 2) %>%
  mutate(lambda = paste0("(",lambda1,", ", lambda2,", ", lambda3,", ", lambda4,")"))%>%
  dplyr::select(-c(lambda1, lambda2, lambda3, lambda4))
keeps_summary$model <- model; keeps_summary$series <- data

#TABLE XXX POSTERIOR MEANS OF PARAMETERS------
keeps_summary[,c(ncol(keeps_summary), ncol(keeps_summary)-1, 1:(ncol(keeps_summary)-2))] %>%
  arrange(series) %>%
  xtable() %>%
  print()

#FIGURE XXX STOCHASTIC VOLATILITY--------
keeps_v %>%
  mutate(model = model, series = data) %>%
  melt(id.vars = c(model, series)) %>%
  ggplot() +
  geom_line(aes(x = Var2, y = value, colour = model)) +
  facet_grid(series~., scales ="free_y")

