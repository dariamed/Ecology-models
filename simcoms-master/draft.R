
### FIT ALL MODELS


##SET UP ##########################################################################################################
library(arm)
library(jagsUI)
library(ggplot2)
library(gridExtra)
library(parallel)
library(Rcpp)
library(magrittr)
library(purrr)
library(readr)
library(corrplot)
library(gjam)
library(Hmsc)
library(rlist)
library(ggmcmc)
library(devtools)
#install_github("sbfnk/fitR")
library(fitR) ###load through github
library(reshape)
library(AUC)
library(ggthemes)
library(mcclust.ext)
Rcpp::sourceCpp('/Users/dariabystrova/Documents/GitHub/gjamed/src/cppFns.cpp')
source("/Users/dariabystrova/Documents/GitHub/gjamed/R/gjamHfunctions_mod.R")
#setwd("C:/Users/giaru/Desktop/Documents/GitHub/Ecology-models/simcoms-master")

setwd("/Users/dariabystrova/Documents/GitHub/Ecology-models/simcoms-master/ExampleFiles")
#setwd("~/Tesi/Code/Ecology-models-master_1/simcoms-master")

# lapply(list.files(path = "."),load,.GlobalEnv)
#setwd("~/Tesi/Code/Ecology-models-master/simcoms-master")
load("params.rds")
load("sim_names.rds")
load("comp_inter.rds")
load("fac_inter.rds")
#sim_data<-readRDS("sim_data.rds")

load_object <- function(file) {
  tmp <- new.env()
  load(file = file, envir = tmp)
  tmp[[ls(tmp)[1]]]
}

setwd("~/Documents/GitHub/Ecology-models/simcoms-master")
sim_data<-readRDS("sim_data.rds")

##SET UP ##########################################################################################################
##JSDM set up  ##########################################################################################################
prob_cooccur_es <- function(Y) {
  K <- ncol(Y)
  ans <- matrix(0, K, K)
  
  for (k in 1:K) {
    for (kk in 1:K) {
      N1 <- sum(Y[, k])
      N2 <- sum(Y[, kk])
      N <- nrow(Y)
      j <- max(0, N1 + N2 - N):min(N1, N2)
      p <- vector("numeric", length(j))
      for (i in seq_along(j)) {
        p[i] <- (
          choose(N, j[i]) * choose(N - j[i], N2 - j[i]) *
            choose(N - N2, N1 - j[i])
        ) / (
          choose(N, N2) * choose(N, N1)
        )
      }  
      ans[k, kk] <- (sum(Y[, k] + Y[, kk] == 2) - sum(p * j)) / N
    }
  }
  ans
}

me5 <- load_object("model-2019-04-09-19-02-16.rda")
jsdm_conv<-function(mod) {
  
  # cat(sprintf("Maximum Rhat value for Beta: %s\n", round(max(mod$R$B),4)))
  #  cat(sprintf("Maximum Rhat value for Rho: %s\n", round(max(mod$R$Rho),4)))
  #  cat(sprintf("Maximum Rhat value for EnvRho: %s\n", round(max(mod$R$EnvRho),4)))
  n_eff_beta<- as.data.frame(t(mod$n.eff$B))
  n_eff_beta$parameter<- "beta"
  n_eff_sigma<- as.data.frame(mod$n.eff$Rho)
  n_eff_sigma$parameter<- "correlation"
  neff<- rbind(n_eff_beta,n_eff_sigma)
  neff_mod<- melt(neff, id=c("parameter"), variable_name = "effective size")
  p<- ggplot(neff_mod, aes(x=value, color=parameter,fill=parameter)) +
    geom_histogram( alpha=0.4, position="identity",binwidth=100) + xlab("effective size") +
    scale_color_brewer(palette="Dark2")+scale_fill_brewer(palette="Dark2")+
    ggtitle("Effective size for the parameters beta and coefficients of correlation matrix")
  plot(p)
  Rhat_beta<- as.data.frame(t(mod$Rhat$B))
  Rhat_beta$parameter<- "beta"
  Rhat_sigma<- as.data.frame(mod$Rhat$Rho)
  Rhat_sigma$parameter<- "correlation"
  Rhat<- rbind(Rhat_beta,Rhat_sigma)
  Rhat_mod<- melt(Rhat, id=c("parameter"), variable_name = "Rhat")
  Rhat_mod_nonNA<- Rhat_mod[!is.na(Rhat_mod$value),]
  p2<- ggplot(Rhat_mod_nonNA, aes(x=value, color=parameter,fill=parameter)) +
    geom_histogram( alpha=0.4, position="identity",binwidth = 0.01) +
    scale_color_brewer(palette="Dark2")+scale_fill_brewer(palette="Dark2")+ xlab("Rhat") +
    ggtitle("Rhat for the parameters beta and the coefficients of correlation matrix")
  plot(p2)
  
}
jsdm_conv(me5)
#me5$mcmc.info[1:7]

metrics_jsdm<-function(model,fac=NULL,comp=NULL,only_env=T){
  
  Rho<-(!model$overlap0$Rho)*model$mean$Rho
  #{
  if(!is.null(comp) & !is.null(fac)){
    fac_comp <- fac-comp
    TP_comp <- FN_comp <-TP_fac <- FN_fac <- FP <- TN <- wrong <- 0
    for(i in 1:(nrow(Rho)-1)){
      for(j in (i+1):ncol(Rho)){
        if(fac_comp[i,j]>0 & Rho[i,j] >0)
          TP_fac=TP_fac+1
        if(fac_comp[i,j]>0 & Rho[i,j] == 0)
          FN_fac=FN_fac+1
        if(fac_comp[i,j]==0 & Rho[i,j] == 0)
          TN=TN+1
        if(fac_comp[i,j]==0 & Rho[i,j] != 0)
          FP=FP+1
        if(fac_comp[i,j]< 0 & Rho[i,j] < 0)
          TP_comp=TP_comp+1
        if(fac_comp[i,j]<0 & Rho[i,j] == 0)
          FN_comp=FN_comp+1
        if((fac_comp[i,j]<0 & Rho[i,j] > 0) | (fac_comp[i,j]>0 & Rho[i,j] < 0))
          wrong=wrong+1
      }
    }
    
  }else{
    if(!is.null(comp)){
      
      TP_comp <- FN_comp <- FP <- TN <- wrong <- 0
      TP_fac <- FN_fac<-NULL
      for(i in 1:(nrow(Rho)-1)){
        for(j in (i+1):ncol(Rho)){
          if(comp[i,j]>0 & Rho[i,j] <0)
            TP_comp=TP_comp+1
          if(comp[i,j]>0 & Rho[i,j] >0)
            wrong=wrong+1
          if(comp[i,j]>0 & Rho[i,j] == 0)
            FN_comp=FN_comp+1
          if(comp[i,j]==0 & Rho[i,j] == 0)
            TN=TN+1
          if(comp[i,j]==0 & Rho[i,j] != 0)
            FP=FP+1
          
        }
        #    }
        #    }
      }
    }
    if(!is.null(fac)){
      
      TP_fac <- FN_fac <- FP <- TN <- wrong <- 0
      TP_comp <- FN_comp<-NULL
      for(i in 1:(nrow(Rho)-1)){
        for(j in (i+1):ncol(Rho)){
          if(fac[i,j]>0 & Rho[i,j] >0)
            TP_fac=TP_fac+1
          if(fac[i,j]>0 & Rho[i,j] <0)
            wrong = wrong+1
          if(fac[i,j]>0 & Rho[i,j] == 0)
            FN_fac=FN_fac+1
          if(fac[i,j]==0 & Rho[i,j] == 0)
            TN=TN+1
          if(fac[i,j]==0 & Rho[i,j] != 0)
            FP=FP+1
        }
      }
    }
  }
  if(only_env){
    
    FP <- TN <- 0
    TP_fac <- FN_fac<-TP_comp <- FN_comp<-wrong<-NULL
    for(i in 1:(nrow(Rho)-1)){
      for(j in (i+1):nrow(Rho)){
        if(Rho[i,j] != 0)
          FP=FP+1
        if(Rho[i,j] == 0)
          TN=TN+1
      }
    }
    
  }
  success_env<-TN/(TN+FP)
  
  if(!is.null(TP_comp)){success_comp <- TP_comp/(TP_comp+FN_comp)}else{ success_comp = NULL}
  
  if(!is.null(TP_fac)) {success_fac <- TP_fac/(TP_fac+FN_fac)}else{ success_fac = NULL}
  
  list("FP"=FP,"TN"=TN,"TP_comp"=TP_comp,"FN_comp"=FN_comp,"TP_fac"=TP_fac,"FN_fac"=FN_fac,
       "wrong"=wrong,"success_env"=success_env,"success_comp"=success_comp,"success_fac"=success_fac)
}
predict_y_jsdm <- function(x, ns = 1000) {
  data   <- if(x$parallel) x$model[[1]]$data() else x$model$data()
  n      <- data$n
  nsp    <- data$J
  n.samples  <- x$mcmc.info$n.samples
  B <- x$sims.list$B
  X <- data$X
  Y <- data$Y
  inprod_mat <- cppFunction(
    "NumericMatrix inprod_mat(NumericMatrix X, NumericMatrix B) {
    
    int nX = X.nrow();
    int nB = B.nrow();
    int nK = X.ncol();
    
    NumericMatrix theta(nX, nB);
    
    for (int i = 0; i < nX; i++) {
    for (int j = 0; j < nB; j++) {
    theta(i, j) = 0;
    for (int k = 0; k < nK; k++) {
    theta(i, j) += X(i, k) * B(j, k);
    }
    }
    }
    
    return(theta);
    
}"
  )
  vapply(
    floor(seq(1, n.samples, length.out = ns)),
    function(x){ inprod_mat(data$X, B[x, ,])},
    matrix(NA_real_, n, nsp)
  )
}

############## JSDM5env########################################################################################################


data<-sim_data$EnvEvenSp5
data <- list(
  Y = subset(data, select = -env),
  X = cbind(1, scale(poly(data$env, 2))),
  covx = cov(cbind(1, scale(poly(data$env, 2)))),
  K = 3,
  J = ncol(data) - 1,
  n = nrow(data),
  I = diag(ncol(data) - 1),
  df = ncol(data)
)

Y_cor<-cor(data$Y)

to_prec<-function(m){
  n<-dim(m)[1]
  Tau_n<-matrix(nrow=n, ncol=n)
  for (j in 1:n) {
    for (k in 1:n){
      Tau_n[j, k] <-  -m[j, k]/sqrt((m[j,j]*m[k,k]))
    }
  }
  return(Tau_n)
}

#Tau_n<-matrix(nrow=dim(model$mean$Tau)[1], ncol=dim(model$mean$Tau)[1])
#Tau_n<-to_prec(me5$mean$Tau)
#Tau_k<-Tau_n*(!(model$q97.5$Tau>0 & model$q2.5$Tau<0))

par(mfrow=c(2,4),oma = c(3, 1, 2, 1))
cols = colorRampPalette(c("dark blue","white","red"))
col2 <- colorRampPalette(c("#4393C3", "#2166AC", "#053061",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#67001F", "#B2182B", "#D6604D", "#F4A582"))

gcols = colorRampPalette(c( "White", "White", "Black"))

corrplot(Y_cor, diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
title("Correlation cor(Y)")
corrplot(me5$mean$EnvRho, diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
title("EnvRho")
corrplot(me5$mean$EnvRho*(!me5$overlap0$EnvRho), diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
title("EnvRho signif")
corrplot(me5$mean$Rho, diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
title("Rho")
corrplot(me5$mean$Rho*(!me5$overlap0$Rho), diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
title("Rho signif")
# corrplot(Tau_n, diag = FALSE, order ="original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
# title("Tau")
# corrplot(Tau_n*(!me5$overlap0$Tau), diag = FALSE, order ="original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
# title("Tau signif")

j_metric_e5<-metrics_jsdm(me5)
cat(sprintf("Success rate: %s\n", metrics_jsdm(me5)$success_env))
mod_list_Rho<-list()
mod_list_Rho<-list(jsdm =as.matrix(me5$mean$Rho*(!me5$overlap0$Rho))) 
mod_list_Tau<-list()
mod_list_Tau<-list(jsdm =as.matrix(me5$mean$Rho*(!me5$overlap0$Rho))) 

####Prediction
pred_j<-pnorm(predict_y_jsdm(me5))
pred_j_mean <- apply(pred_j, 1:2, mean)

pred_j_05 <- apply(pred_j, 1:2, quantile,0.05)
pred_j_95 <- apply(pred_j, 1:2, quantile,0.95)
pred_jsdm<-list(pred_j_mean=pred_j_mean,pred_j_05=pred_j_05,pred_j_95=pred_j_95)

#prepare AUC
AUC_j_env<-vector()
for(i in 1:data$J) AUC_j_env<-c(AUC_j_env,auc(roc(pred_j_mean[,i],factor(data$Y[,i]))))

##JSDM5env ##########################################################################################################
###GJAM set up ##########################################################################################################


###GJAM set up ##########################################################################################################
