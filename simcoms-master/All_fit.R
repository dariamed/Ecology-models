  ####SETUP########################################################################
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
library(fitR) ###load through github
library(reshape)
library(AUC)
library(ggthemes)
library(mcclust.ext)
Rcpp::sourceCpp('/Users/dariabystrova/Documents/GitHub/gjamed/src/cppFns.cpp')
source("/Users/dariabystrova/Documents/GitHub/gjamed/R/gjamHfunctions_mod.R")


#setwd("C:/Users/giaru/Desktop/Documents/GitHub/Ecology-models/simcoms-master")

setwd("/Users/dariabystrova/Documents/GitHub/Ecology-models/simcoms-master/ExampleFiles")
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

####SETUP########################################################################



### JSDM Functions##############################################################################################################################################

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

jsdm_conv<-function(mod) {
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


### JSDM Functions##############################################################################################################################################




### GJAM Functions##############################################################################################################################################
makeSymm <- function(m) {
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  return(m)
}

expandSigma_rmd <- function(sigma, S){
  
  ss <- diag(S)
  ss[lower.tri(ss,diag=T)] <- sigma
  ss[upper.tri(ss)] <- t(ss)[upper.tri(ss)]
  ss
}


convert_to_m<-function(ar){
  d <-floor((sqrt(length(ar)*8+1)-1)/2)
  C <- matrix(0,d,d)
  i.lwr <- which(lower.tri(C, diag = TRUE), arr.ind=TRUE)
  C[i.lwr] <- ar
  C<-makeSymm(C)
  return(t(C))
}


metrics_gjam<-function(Rho, comp=NULL,fac=NULL,only_env=T){
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


fit_gjam<-function(data, it=2500,burn=500 , name="./gjam_models/temp.rda",interact=diag(ncol(data$Y))){
  #setup parameters
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
  xdata<-as.data.frame(data$X[,-1])
  colnames(xdata)<- c("env","env2")
  ydata<-as.data.frame(data$Y)
  ###formula
  formula<-as.formula( ~env+ env2)
  ml   <- list(ng = it, burnin = burn, typeNames = 'PA', PREDICTX = F)
  ####fit
  
  mod_gjam1  <- gjam(formula, xdata = xdata, ydata = ydata, modelList = ml)
  mod_gjam2  <- gjam(formula, xdata = xdata, ydata = ydata, modelList = ml)
  gjam_mods<- list(m1=mod_gjam1,m2=mod_gjam2)
  
  save(gjam_mods, file = name)
  gjam_bs<- mcmc.list(mcmc(gjam_mods$m1$chains$bgibbsUn,start=burn, end=it),mcmc(gjam_mods$m2$chains$bgibbsUn,start=burn, end=it))
  gjam_sigma<- mcmc.list(mcmc(gjam_mods$m1$chains$sgibbs,start=burn, end=it),mcmc(gjam_mods$m2$chains$sgibbs,start=burn, end=it))
  
  
  hist(effectiveSize(gjam_bs), main="ess(beta)",lwd=2,col=gray(.6))
  hist(gelman.diag(gjam_bs, multivariate=FALSE)$psrf,lwd=2,col=gray(.6), main="psrf(beta)")
  #gelman.plot(gjam_bs)
  #xyplot(gjam_bs)
  #traceplot(gjam_sigma)
  hist(effectiveSize(gjam_sigma), main="ess(beta)",lwd=2,col=gray(.6))
  hist(gelman.diag(gjam_sigma,multivariate=FALSE)$psrf,lwd=2,col=gray(.6), main="psrf(sigma)")
  #gelman.plot(gjam_sigma)
  #xyplot(gjam_sigma)
  
  gjam_mods_2bgibbs<-  abind(gjam_mods$m1$chains$bgibbsUn,gjam_mods$m2$chains$bgibbsUn, along = 3)
  gjam_mods_2sgibbs<-abind(gjam_mods$m1$chains$sgibbs,gjam_mods$m2$chains$sgibbs, along = 3)
  
  #postH<-apply(mod_gjam1$chains$sgibbs, 2, quantile,0.95)
  #postL<-apply(mod_gjam1$chains$sgibbs, 2, quantile,0.05)
  postH<- apply(gjam_mods_2sgibbs, 2, quantile,0.95)
  postL<-apply(gjam_mods_2sgibbs, 2, quantile,0.05)
  post_mean<-apply(gjam_mods_2sgibbs, 2, mean)
  
  
  pH<-convert_to_m(postH)
  pL<-convert_to_m(postL)
  post_mean_s<-convert_to_m(post_mean)
  #R_sign<-cov2cor(mod_gjam1$parameters$sigMu)*(!(pH>0 & pL<0))
  R_sign<-cov2cor(post_mean_s)*(!(pH>0 & pL<0))
  
  sgibbs<-abind(mod_gjam1$chains$sgibbs[-(1:burn),],mod_gjam2$chains$sgibbs[-(1:burn),],along=1)
  
  tau<-array(NA,dim=c(data$J,data$J,dim(sgibbs)[1]))
  for(j in 1:dim(sgibbs)[1]){
    ss <- expandSigma_rmd(sgibbs[j,], S = data$J)
    si <- solve(ss)
    tau[,,j] <- -cov2cor(si)
  }
  
  tau_mean<-apply(tau,c(1,2), mean)
  tau_HI<-apply(tau,c(1,2),quantile,0.95)
  tau_LO<-apply(tau,c(1,2),quantile,0.05)
  
  Tau_sign<-tau_mean*(!(tau_HI>0 & tau_LO<0))
  
  #newdata   <- list(xdata = xdata, nsim = 50 )
  
  x<- data$X
  #p1<-gjamPredict(mod_gjam1)
  mu<-array(NA,dim=c(data$n,data$J,it))
  for(k in 1:it){
    for(j in 1:data$J){
      
      mu[,j,k] <- pnorm(x%*%mod_gjam1$chains$bgibbs[k,(3*(j-1)+1):(3*j)])
      
    }
  }
  
  return(list(Rho_sign=R_sign,Tau_sign=Tau_sign,predict=mu))
  
}




###this function only loads the model and return the R and T significant
load_gjam<-function(data,it=2500,burn=500,name="./gjam_models/temp.rda",interact=diag(ncol(data$Y))){
  #setup parameters
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
  
  xdata<-as.data.frame(data$X[,-1])
  
  gj_mod<-load_object(name)
  S<-ncol(data$Y)
  gjam_bs<- mcmc.list(mcmc(gj_mod$m1$chains$bgibbsUn[-(1:burn),]),mcmc(gj_mod$m2$chains$bgibbsUn[-(1:burn),]))
  gjam_sigma<- mcmc.list(mcmc(gj_mod$m1$chains$sgibbs[-(1:burn),]),mcmc(gj_mod$m2$chains$sgibbs[-(1:burn),]))
  ###NEW plot for effective size
  n_eff_beta<- as.data.frame(effectiveSize(gjam_bs))
  colnames(n_eff_beta)<- c("value")
  n_eff_beta$parameter<- "beta"
  n_eff_sigma<- as.data.frame(effectiveSize(gjam_sigma))
  colnames(n_eff_sigma)<- c("value")
  n_eff_sigma$parameter<- "correlation"
  neff<- rbind(n_eff_beta,n_eff_sigma)
  p<- ggplot(neff, aes(x=value, color=parameter,fill=parameter)) +
    geom_histogram( alpha=0.4, position="identity") + xlab("effective size") +
    scale_color_brewer(palette="Dark2")+scale_fill_brewer(palette="Dark2")+
    ggtitle("Effective size for the parameters beta and coefficients of correlation matrix")
  plot(p)
  
  Rhat_beta<-as.data.frame(gelman.diag(gjam_bs, multivariate=FALSE)$psrf)
  Rhat_beta$parameter<- "beta"
  Rhat_sigma<- as.data.frame(gelman.diag(gjam_sigma, multivariate=FALSE)$psrf)
  Rhat_sigma$parameter<- "correlation"
  Rhat<- rbind(Rhat_beta,Rhat_sigma)
  p2<- ggplot(Rhat, aes(x= Rhat$`Point est.`, color=parameter,fill=parameter)) +
    geom_histogram( alpha=0.4, position="identity",binwidth =1) +
    scale_color_brewer(palette="Dark2")+scale_fill_brewer(palette="Dark2")+ xlab("Rhat") +
    ggtitle("Rhat for the parameters beta and the coefficients of correlation matrix")
  plot(p2)
  
  #gelman.plot(gjam_bs)
  #xyplot(gjam_bs)
  #traceplot(gjam_sigma)
  #gelman.plot(gjam_sigma)
  #xyplot(gjam_sigma)
  
  gjam_mods_2bgibbs<-  abind(gj_mod$m1$chains$bgibbsUn,gj_mod$m2$chains$bgibbsUn, along = 3)
  gjam_mods_2sgibbs<-abind(gj_mod$m1$chains$sgibbs,gj_mod$m2$chains$sgibbs, along = 3)
  
  #postH<-apply(mod_gjam1$chains$sgibbs, 2, quantile,0.95)
  #postL<-apply(mod_gjam1$chains$sgibbs, 2, quantile,0.05)
  postH<- apply(gjam_mods_2sgibbs, 2, quantile,0.95)
  postL<-apply(gjam_mods_2sgibbs, 2, quantile,0.05)
  post_mean<-apply(gjam_mods_2sgibbs, 2, mean)
  
  
  pH<-convert_to_m(postH)
  pL<-convert_to_m(postL)
  post_mean_s<-convert_to_m(post_mean)
  S_mean<-cov2cor(post_mean_s)
  #R_sign<-cov2cor(mod_gjam1$parameters$sigMu)*(!(pH>0 & pL<0))
  R_sign<-cov2cor(post_mean_s)*(!(pH>0 & pL<0))
  
  sgibbs<-abind(gj_mod$m1$chains$sgibbs[-(1:burn),],gj_mod$m2$chains$sgibbs[-(1:burn),],along=1)
  
  tau<-array(NA,dim=c(data$J,data$J,dim(sgibbs)[1]))
  for(j in 1:dim(sgibbs)[1]){
    ss <- expandSigma_rmd(sgibbs[j,], S = data$J)
    si <- solve(ss)
    tau[,,j] <- -cov2cor(si)
  }
  
  tau_mean<-apply(tau,c(1,2), mean)
  tau_HI<-apply(tau,c(1,2),quantile,0.95)
  tau_LO<-apply(tau,c(1,2),quantile,0.05)
  
  Tau_sign<-tau_mean*(!(tau_HI>0 & tau_LO<0))
  
  x<- data$X
  #p1<-gjamPredict(mod_gjam1)
  mu<-array(NA,dim=c(data$n,data$J,it))
  for(k in 1:it){
    for(j in 1:data$J){
      
      mu[,j,k] <- pnorm(x%*%gj_mod$m1$chains$bgibbs[k,(3*(j-1)+1):(3*j)])
      
    }
  }
  
  par(mfrow=c(2,3),oma = c(1, 1, 1, 1))
  corrplot(cor(data$Y), diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
  title("Correlation cor(Y)")
  corrplot(S_mean, diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
  title("R")
  corrplot(gj_mod$m1$parameters$ematrix, diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
  title("E matrix")
  corrplot(Tau_sign, diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
  title("Tau signif")
  corrplot(R_sign, diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
  title("R signif")
  corrplot(interact, diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
  title("True interactions")
  
  return(list(Rho_sign=R_sign,Tau_sign=Tau_sign,predict=mu))
}


### GJAM Functions##############################################################################################################################################


### HMSC Functions##############################################################################################################################################


fit_hmsc<-function(data,label="Fit",nsamples = 1000,nchains=2,name="./HMmodels/hmtemp.rda" ){
  if (label=="Fit"){
    Y_data = subset(data, select = -env)
    ns<- ncol(Y_data)
    np <- nrow(Y_data)
    X<-scale(poly(data$env[1:np], 2))
    colnames(X)<-c("env","env2")
    studyDesign = data.frame(sample = as.factor(1:np))
    rL = HmscRandomLevel(units = studyDesign$sample)
    m = Hmsc(Y=as.matrix(Y_data), XData=as.data.frame(X), XFormula=~env+env2, distr="probit",
             studyDesign = studyDesign, ranLevels = list(sample = rL))
    m = sampleMcmc(m, nsamples, thin=10, adaptNf=c(200,200), transient=500,nChains=nchains ,verbose=F)
    save(m, file = name)
    return(m)
  }
  if (label=="Load"){
    return(load_object(name))
  }
}

hm_conv<-function(mod){
  codaList = convertToCodaObject(mod)
  #convergence histograms
  ##NEW convergence histograms
  n_eff_beta<- as.data.frame(effectiveSize(codaList$Beta))
  colnames(n_eff_beta)<- c("value")
  n_eff_beta$parameter<- "beta"
  n_eff_sigma<- as.data.frame(effectiveSize(codaList$Omega[[1]]))
  colnames(n_eff_sigma)<- c("value")
  n_eff_sigma$parameter<- "correlation"
  neff<- rbind(n_eff_beta,n_eff_sigma)
  p<- ggplot(neff, aes(x=value, color=parameter,fill=parameter)) +
    geom_histogram( alpha=0.4, position="identity",binwidth =50) + xlab("effective size") +
    scale_color_brewer(palette="Dark2")+scale_fill_brewer(palette="Dark2")+
    ggtitle("Effective size for the parameters beta and coefficients of correlation matrix")
  plot(p)
  
  Rhat_beta<-as.data.frame(gelman.diag(codaList$Beta,multivariate=FALSE)$psrf)
  Rhat_beta$parameter<- "beta"
  Rhat_sigma<- as.data.frame(gelman.diag(codaList$Omega[[1]], multivariate=FALSE)$psrf)
  Rhat_sigma$parameter<- "correlation"
  Rhat<- rbind(Rhat_beta,Rhat_sigma)
  p2<- ggplot(Rhat, aes(x= Rhat$`Point est.`, color=parameter,fill=parameter)) +
    geom_histogram( alpha=0.4, position="identity",binwidth =0.01) +
    scale_color_brewer(palette="Dark2")+scale_fill_brewer(palette="Dark2")+ xlab("Rhat") +
    ggtitle("Rhat for the parameters beta and the coefficients of correlation matrix")
  plot(p2)
}


metrics_hmsc<-function(Rho, comp=NULL,fac=NULL,only_env=T){
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
      }
      #  }
      #  }
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


hm_inter<-function(mod, nchains=2,nsamples = 1000, interact=diag(ns)){
  getOmega = function(a,r=1)
    return(crossprod(a$Lambda[[r]]))
  ns<-mod$ns
  postOmega1 = array(unlist(lapply(mod$postList[[1]],getOmega)),c(ns,ns,mod$samples))
  postOmega2 = array(unlist(lapply(mod$postList[[2]],getOmega)),c(ns,ns,mod$samples))
  
  postOmega<-abind(postOmega1,postOmega2,along=3)
  postOmegaMean = apply(postOmega,c(1,2),mean)
  postOmegaUp=apply(postOmega,c(1,2),quantile,0.95)
  postOmegaLo=apply(postOmega,c(1,2),quantile,0.05)
  
  postR<-array(dim=c(ns,ns,nchains*nsamples))
  postT<-array(dim=c(ns,ns,nchains*nsamples))
  for(i in 1:dim(postOmega)[3]){
    postR[,,i]<-stats::cov2cor(postOmega[,,i])
    postT[,,i]<- -stats::cov2cor(solve(postOmega[,,i]+diag(ns)))
  }
  
  postRMean = apply(postR,c(1,2),mean)
  postRUp=apply(postR,c(1,2),quantile,0.95)
  postRLo=apply(postR,c(1,2),quantile,0.05)
  
  postTMean = apply(postT,c(1,2),mean)
  postTUp=apply(postT,c(1,2),quantile,0.95)
  postTLo=apply(postT,c(1,2),quantile,0.05)
  
  Toplot_R<-postRMean*(!(postRUp>0 & postRLo<0))
  Toplot_T<-postTMean*(!(postTUp>0 & postTLo<0))
  
  # Omegacor<- computeAssociations(m)
  # supportLevel<- 0.95
  # toPlot<- ((Omegacor[[1]]$support>supportLevel)+ (Omegacor[[1]]$support<(1-supportLevel))>0)*Omegacor[[1]]$mean
  # corrplot(toPlot, method="color", col=colorRampPalette(c("blue", "white", "red"))(200))
  
  par(mfrow=c(2,3),oma = c(1, 1, 1, 1))
  corrplot(cor(hm_mod$Y), diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
  title("Correlation cor(Y)")
  corrplot(postRMean, diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
  title("R")
  corrplot(Toplot_R, diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
  title("Plot only non zero value")
  corrplot(Toplot_T, diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
  title("Partial correlation matrix")
  corrplot(interact, diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
  title("True interactions")
  
  return(list(Rho_sign=Toplot_R,Tau_sign=Toplot_T)) 
}





### HMSC Functions##############################################################################################################################################

### DR-GJAM Functions##############################################################################################################################################
gjam_dim_red_5<-function(datab,ng=2500, burnin=500,name="./gjam_models/tmpdr.rda", regime="L"){
  ydata<-subset(datab, select = -env)
  xdata<-scale(poly(datab$env, 2))
  ns<- ncol(ydata)
  N<-ncol(ydata)-1
  n<-nrow(datab)
  colnames(xdata)<- c("env","env2")
  formula<- ~env + env2
  r<-3
  
  rl <-list(N=N, r=r)
  ml  <- list(ng = ng, burnin = burnin, typeNames = 'PA',reductList=rl,PREDICTX = F )
  if(regime=="L"){
    gjam_dr_mods<-load_object(name)
    mod_gjam_red_1<-gjam_dr_mods$m1
    mod_gjam_red_2<-gjam_dr_mods$m2
  }else{
    mod_gjam_red_1 <- gjam(formula, xdata, ydata, modelList = ml)
    mod_gjam_red_2 <- gjam(formula, xdata, ydata, modelList = ml)
    gjam_dr_mods<- list(m1=mod_gjam_red_1,m2=mod_gjam_red_2)
    save(gjam_dr_mods, file = name)
  }
  gjam_bs<- mcmc.list(mcmc(mod_gjam_red_1$chains$bgibbs[-(1:burnin),]),mcmc(mod_gjam_red_2$chains$bgibbsUn[-(1:burnin),]))
  gjam_sigma<- mcmc.list(mcmc(mod_gjam_red_1$chains$sgibbs[-(1:burnin),]),mcmc(mod_gjam_red_2$chains$sgibbs[-(1:burnin),]))
  n_eff_beta<- as.data.frame(effectiveSize(gjam_bs))
  colnames(n_eff_beta)<- c("value")
  n_eff_beta$parameter<- "beta"
  n_eff_sigma<- as.data.frame(effectiveSize(gjam_sigma))
  colnames(n_eff_sigma)<- c("value")
  n_eff_sigma$parameter<- "correlation"
  neff<- rbind(n_eff_beta,n_eff_sigma)
  p<- ggplot(neff, aes(x=value, color=parameter,fill=parameter)) +
    geom_histogram( alpha=0.4, position="identity",binwidth = 100) + xlab("effective size") +
    scale_color_brewer(palette="Dark2")+scale_fill_brewer(palette="Dark2")+
    ggtitle("Effective size for the parameters beta and coefficients of correlation matrix")
  
  plot(p)
  
  Rhat_beta<-as.data.frame(gelman.diag(gjam_bs, multivariate=FALSE)$psrf)
  Rhat_beta$parameter<- "beta"
  Rhat_sigma<- as.data.frame(gelman.diag(gjam_sigma, multivariate=FALSE)$psrf)
  Rhat_sigma$parameter<- "correlation"
  Rhat<- rbind(Rhat_beta,Rhat_sigma)
  p2<- ggplot(Rhat, aes(x= Rhat$`Point est.`, color=parameter,fill=parameter)) +
    geom_histogram( alpha=0.4, position="identity",binwidth =0.1) +
    scale_color_brewer(palette="Dark2")+scale_fill_brewer(palette="Dark2")+ xlab("Rhat") +
    ggtitle("Rhat for the parameters beta and the coefficients of correlation matrix")
  plot(p2)
  
  
  #gelman.plot(gjam_bs)
  #xyplot(gjam_bs)
  #traceplot(gjam_sigma)
  #gelman.plot(gjam_sigma)
  #xyplot(gjam_sigma)
  
  sgibbs<-abind(mod_gjam_red_1$chains$sgibbs[-(1:burnin),],mod_gjam_red_2$chains$sgibbs[-(1:burnin),],along=1)
  sigErrGibbs<-abind(mod_gjam_red_1$chains$sigErrGibbs[-(1:burnin)],mod_gjam_red_2$chains$sigErrGibbs[-(1:burnin)],along=1)
  kgibbs<-abind(mod_gjam_red_1$chains$kgibbs[-(1:burnin),],mod_gjam_red_2$chains$kgibbs[-(1:burnin),],along=1)
  sigma<-invsigma<-array(NA,dim=c(ns,ns,2*(ng-burnin)))
  
  
  #sgibbs<-mod_gjam_red$chains$sgibbs
  #sigErrGibbs<-mod_gjam_red$chains$sigErrGibbs
  #kgibbs<-mod_gjam_red$chains$kgibbs
  N<-mod_gjam_red_1$modelList$reductList$N
  r<-mod_gjam_red_1$modelList$reductList$r
  N_dim<-2*(ng-burnin)
  for(j in 1:N_dim){
    Z  <- matrix(sgibbs[j,],N,r)
    sigma[,,j] <- .expandSigma(sigErrGibbs[j], ns, Z = Z, kgibbs[j,], REDUCT = T) #sigma
    invsigma[,,j] <- invWbyRcpp(sigErrGibbs[j], Z[kgibbs[j,],]) #inverse sigma
  } 
  
  sigma_mean<-apply(sigma,c(1,2),mean) 
  sigma_q05<-apply(sigma,c(1,2),quantile,0.05) 
  sigma_q95<-apply(sigma,c(1,2),quantile,0.95) 
  Sigma_sign<--cov2cor(sigma_mean*(!(sigma_q95>0 & sigma_q05<0)))
  
  invsigma_mean<-apply(invsigma,c(1,2),mean) 
  invsigma_q05<-apply(invsigma,c(1,2),quantile,0.05) 
  invsigma_q95<-apply(invsigma,c(1,2),quantile,0.95) 
  INVSigma_sign<--cov2cor(sigma_mean*(!(invsigma_q95>0 & invsigma_q05<0)))
  x <- cbind(1, scale(poly(datab$env, 2)))
  #x<- xdata
  J<- ncol(datab) - 1
  
  bgibbspr<-abind(mod_gjam_red_1$chains$bgibbs[-(1:burnin),],mod_gjam_red_2$chains$bgibbs[-(1:burnin),],along=1)
  N_dim<-2*(ng-burnin)
  mu<-array(NA,dim=c(n,J,N_dim))
  for(k in 1:N_dim){
    for(j in 1:J){
      
      mu[,j,k] <- pnorm(x%*%bgibbspr[k,(3*(j-1)+1):(3*j)])
      
    }
  }
  
  return(list(Rho_sign_d=Sigma_sign,Tau_sign_d=INVSigma_sign,k=kgibbs, predict=mu, bchain=bgibbspr))
}






### DR-GJAM Functions##############################################################################################################################################





###  Common functions##############################################################################################################################################
### Plot species responses


create_plot_response<- function(data, nsp=5,pred_gjam, pred_jsdm,pred_hmsc,gjam_dr){
  np<-nrow(data)
  nspecies<-nsp
  niche_optima  = seq(2, 98, length.out=nspecies)
  niche_breadth = 20
  
  #dataframe for fundamental niches
  table_fundamental<-data.frame()
  
  for(i in 1:nspecies){
    tmp<-data.frame(xx=0:100,niche=dnorm(0:100,mean=niche_optima[i],sd=niche_breadth)/dnorm(niche_optima[i],mean=niche_optima[i],sd=niche_breadth),species=i)
    table_fundamental<-rbind(table_fundamental,tmp)
  }
  
  ######## JSDM
  xx<-data$env
  
  table_jsdm<-data.frame()
  pred_j_mean<-pred_jsdm$pred_j_mean
  pred_j_05<-pred_jsdm$pred_j_05
  pred_j_95<-pred_jsdm$pred_j_95
  
  
  for(i in 1:nspecies) {
    tmp<-data.frame(xx=xx,mean=pred_j_mean[,i],q_95=pred_j_95[,i],q_05=pred_j_05[,i],type=rep("jsdm",np),species=rep(i,np))
    table_jsdm<-rbind(table_jsdm,tmp)
    
  }
  ######## GJAM
  table_gjam<-data.frame()
  
  pred_gj_mean<-pred_gjam$pred_gj_mean
  pred_gj_95<-pred_gjam$pred_gj_95
  pred_gj_05<-pred_gjam$pred_gj_05
  
  for(i in 1:nspecies) {
    tmp<-data.frame(xx=xx,mean=pred_gj_mean[,i],q_95=pred_gj_95[,i],q_05=pred_gj_05[,i],type=rep("gjam",np),species=rep(i,np))
    table_gjam<-rbind(table_gjam,tmp)
  }
  
  
  ######## GJAM_DR
  table_gjam_dr<-data.frame()
  
  pred_gj_dr_mean<-pred_gjam_dr$pred_gj_dr_mean
  pred_gj_dr_95<-pred_gjam_dr$pred_gj_dr_95
  pred_gj_dr_05<-pred_gjam_dr$pred_gj_dr_05
  
  for(i in 1:nspecies) {
    tmp<-data.frame(xx=xx,mean=pred_gj_dr_mean[,i],q_95=pred_gj_dr_95[,i],q_05=pred_gj_dr_05[,i],type=rep("gjam_dr",np),species=rep(i,np))
    table_gjam_dr<-rbind(table_gjam_dr,tmp)
  }
  
  
  #HMSC
  table_hmsc<-data.frame()
  pred_hm_mean<-pred_hmsc$pred_hm_mean
  pred_hm_95<-pred_hmsc$pred_hm_95
  pred_hm_05<-pred_hmsc$pred_hm_05
  
  for(i in 1:nspecies) {
    tmp<-data.frame(xx=xx,mean=pred_hm_mean[,i],q_95=pred_hm_95[,i],q_05=pred_hm_05[,i],type=rep("hmsc",np),species=rep(i,np))
    table_hmsc<-rbind(table_hmsc,tmp)
  }
  
  #table for predictions
  table<-rbind(table_jsdm,table_gjam,table_gjam_dr,table_hmsc)
  
  #table for observations
  Y_data = subset(data, select = -env)
  table_obs<-data.frame()
  
  for(i in 1:nspecies){
    tmp<-data.frame(xx=xx,obs=Y_data[,i],species=rep(i,np))
    table_obs<-rbind(table_obs,tmp)
  }
  
  
  ######## FIRST TYPE OF PLOT
  # 1 plot for each species, so one figure with 5 plots, in total 3 figures, one for each model
  p_1<-list()
  #jsdm 
  for(i in 1:nsp)
    local({
      i<-i
      tmp<-table_jsdm[which(table_jsdm$species==i),]
      tmp_obs<-table_obs[which(table_obs$species==i),]
      tmp_fund<-table_fundamental[which(table_fundamental$species==i),]
      g<<-ggplot()+
        geom_ribbon(aes(x=tmp$xx,ymin=tmp$q_05,ymax=tmp$q_95),alpha=0.5)+
        geom_line(aes(x=tmp$xx,y=tmp$mean,color = "Predicted probability"),lwd=1.5)+
        geom_point(aes(x=tmp_obs$xx,y=tmp_obs$obs),col="#000066",size=0.5) +xlab("Environmental gradient")+ylab("Probability of presence")+
        geom_line(data=tmp_fund,aes(x=tmp_fund$xx,y=tmp_fund$niche,color = "Fundamental niche"),lwd=1)+
        labs(title=paste0("JSDM, Species ",i))+
        scale_color_manual(name = c("Legend"), values = c("Predicted probability" = "#FF6666","Fundamental niche"="#9999FF"))
      p_1<<-list.append(p_1,g)
      assign(paste0("p",i), g, pos =1)
    })
  
  if(nsp==5){
    grid.arrange(p1,p2,p3,p4,p5,nrow=ceiling(nspecies/2))
  } else{ p_1}
  
  
  #gjam
  p_2<-list()
  for(i in 1:nsp)
    local({
      i<-i
      tmp<-table_gjam[which(table_gjam$species==i),]
      tmp_obs<-table_obs[which(table_obs$species==i),]
      tmp_fund<-table_fundamental[which(table_fundamental$species==i),]
      g<<-ggplot()+
        geom_ribbon(aes(x=tmp$xx,ymin=tmp$q_05,ymax=tmp$q_95),alpha=0.5)+
        geom_line(aes(x=tmp$xx,y=tmp$mean,color = "Predicted probability"),lwd=1.5)+
        geom_point(aes(x=tmp_obs$xx,y=tmp_obs$obs),col="#000066",size=0.5) +xlab("Environmental gradient")+ylab("Probability of presence")+
        geom_line(data=tmp_fund,aes(x=tmp_fund$xx,y=tmp_fund$niche,color = "Fundamental niche"),lwd=1)+
        labs(title=paste0("GJAM, Species ",i))+
        scale_color_manual(name = c("Legend"), values = c("Predicted probability" = "#FF6666","Fundamental niche"="#9999FF"))
      
      p_2<<-list.append(p_2,g)
      assign(paste0("p",i), g, pos =1)
    })
  
  if (nsp==5){
    grid.arrange(p1,p2,p3,p4,p5,nrow=ceiling(nspecies/2))
  }else{ 
    p_2
  }
  
  #gjam_dr
  p_2_1<-list()
  for(i in 1:nsp)
    local({
      i<-i
      tmp<-table_gjam_dr[which(table_gjam_dr$species==i),]
      tmp_obs<-table_obs[which(table_obs$species==i),]
      tmp_fund<-table_fundamental[which(table_fundamental$species==i),]
      g<<-ggplot()+
        geom_ribbon(aes(x=tmp$xx,ymin=tmp$q_05,ymax=tmp$q_95),alpha=0.5)+
        geom_line(aes(x=tmp$xx,y=tmp$mean,color = "Predicted probability"),lwd=1.5)+
        geom_point(aes(x=tmp_obs$xx,y=tmp_obs$obs),col="#000066",size=0.5) +xlab("Environmental gradient")+ylab("Probability of presence")+
        geom_line(data=tmp_fund,aes(x=tmp_fund$xx,y=tmp_fund$niche,color = "Fundamental niche"),lwd=1)+
        labs(title=paste0("GJAM DR, Species ",i))+
        scale_color_manual(name = c("Legend"), values = c("Predicted probability" = "#FF6666","Fundamental niche"="#9999FF"))
      
      p_2_1<<-list.append(p_2_1,g)
      assign(paste0("p",i), g, pos =1)
    })
  
  if (nsp==5){
    grid.arrange(p1,p2,p3,p4,p5,nrow=ceiling(nspecies/2))
  }else{ 
    p_2_1
  }
  
  
  
  #hmsc
  p_3<-list()
  for(i in 1:nsp)
    local({
      i<-i
      tmp<-table[which(table_hmsc$species==i),]
      tmp_obs<-table_obs[which(table_obs$species==i),]
      tmp_fund<-table_fundamental[which(table_fundamental$species==i),]
      g<<-ggplot()+
        geom_ribbon(aes(x=tmp$xx,ymin=tmp$q_05,ymax=tmp$q_95),alpha=0.5)+
        geom_line(aes(x=tmp$xx,y=tmp$mean,color = "Predicted probability"),lwd=1.5)+
        geom_point(aes(x=tmp_obs$xx,y=tmp_obs$obs),col="#000066",size=0.5) +xlab("Environmental gradient")+ylab("Probability of presence")+
        geom_line(data=tmp_fund,aes(x=tmp_fund$xx,y=tmp_fund$niche,color = "Fundamental niche"),lwd=1)+
        labs(title=paste0("HMSC, Species ",i))+
        scale_color_manual(name = c("Legend"), values = c("Predicted probability" = "#FF6666","Fundamental niche"="#9999FF"))
      p_3<<-list.append(p_3,g)
      assign(paste0("p",i), g, pos =1)
    })
  
  #x11()
  if (nsp==5){
    grid.arrange(p1,p2,p3,p4,p5,nrow=ceiling(nspecies/2))
  }else{ 
    p_3
  }
  
  
  
  
  ########SECOND TYPE OF PLOT
  p_4<-list()
  for(i in 1:nsp)
    local({
      i<-i
      tmp<-table[which(table$species==i),]
      tmp_obs<-table_obs[which(table_obs$species==i),]
      tmp_fund<-table_fundamental[which(table_fundamental$species==i),]
      g<<-ggplot()+
        geom_ribbon(aes(x=tmp$xx,ymin=tmp$q_05,ymax=tmp$q_95,col=as.factor(tmp$type)),alpha=0.5)+
        geom_line(aes(x=tmp$xx,y=tmp$mean,col=as.factor(tmp$type)),lwd=1.5)+
        geom_point(aes(x=tmp_obs$xx,y=tmp_obs$obs),col="#000066",size=0.5) +xlab("Environmental gradient")+ylab("Probability of presence")+
        geom_line(data=tmp_fund,aes(x=tmp_fund$xx,y=tmp_fund$niche,color = "Fundamental niche"),lwd=1)+
        labs(title=paste0("All models, Species ",i))+
        scale_color_manual(name = c("Legend"), values = c("jsdm" = "#FF6666","gjam" = "#66CC99","gjam_dr" = "#FFFF10","hmsc" = "#FFB266","Fundamental niche"="#9999FF"))
      p_4<<-list.append(p_4,g)
      assign(paste0("p",i), g, pos =1)
    })
  
  if (nsp==5){
    grid.arrange(p1,p2,p3,p4,p5,nrow=ceiling(nspecies/2))
  }else{ 
    p_4
  }
  #x11()
  #grid.arrange(p1,p2,p3,p4,p5,nrow=ceiling(nspecies/2))
  
}



prepare_table<- function(datab,name,class){
  data_new<-as.data.frame(datab) 
  colnames(data_new) <- c("inter","env","env2")
  data_new$sp<-as.character(seq.int(nrow(data_new)))
  data_out<- melt(data_new, id=c("sp"), variable_name = "beta")
  colnames(data_out)[colnames(data_out)=="value"] <-name
  data_out$model<-class
  return(data_out)
}

########GJAM####################################################
prepare_beta_gjam<-function(datab,name,class){
  beta_out<-melt(datab)
  beta_out$rn<-rownames(beta_out)
  beta_out$sp<- as.character(as.numeric(substr(beta_out$rn,3,4)))
  beta_out$beta<- sapply(beta_out$rn,function(x) strsplit(x,"_")[[1]][2])
  colnames(beta_out)[colnames(beta_out)=="value"] <-name
  beta_out$model<- class
  beta_out$beta[which(beta_out$beta=="intercept")] <-"inter"
  return(beta_out[,c("beta","sp","model",name)])
}



####hmsc

prepare_beta_hmsc<-function(datab,name,class,vcol,vrow){
  colnames(datab)<- vcol
  rownames(datab)<-vrow
  beta_out<-melt(datab,value.name = "mean",varnames=c('beta', 'sp'))
  beta_out$sp<- as.character(as.numeric(substr(beta_out$sp,4,5)))
  colnames(beta_out)[colnames(beta_out)=="value"] <-name
  beta_out$model<- class
  beta_out$beta<-as.character(beta_out$beta)
  beta_out$beta[which(beta_out$beta=="(Intercept)")] <-"inter"
  colnames(beta_out)[colnames(beta_out)=="value"] <-name
  return(beta_out[,c("beta","sp","model",name)])
}



create_beta_plot<- function(gjmod, jsdmmod, hmmod, gjam_dr_bchain){
  nsp<-hmmod$ns
  beta_mean<-prepare_table(jsdmmod$mean$B,name="mean",class="jsdm")
  beta_q25<-prepare_table(datab=jsdmmod$q2.5$B,name="low",class="jsdm")
  beta_q97<-prepare_table(jsdmmod$q97.5$B,name="high",class="jsdm")
  beta_jsdm<-merge(beta_mean,beta_q25, by=c("beta","sp","model"))
  beta_jsdm<-merge(beta_jsdm,beta_q97,by=c("beta","sp","model"))
  
  ########GJAM###################################################
  gjmod_1<-gjmod$m1
  beta_mean_raw<- apply(gjmod_1$chains$bgibbs,2,mean)
  beta_q25_raw<- apply(gjmod_1$chains$bgibbs,2,quantile,0.025)
  beta_q97_raw<- apply(gjmod_1$chains$bgibbs,2,quantile,0.975)
  
  beta_mean<- prepare_beta_gjam(beta_mean_raw,"mean", "gjam")
  beta_q25<- prepare_beta_gjam(beta_q25_raw,"low", "gjam")
  beta_q97<- prepare_beta_gjam(beta_q97_raw,"high", "gjam")
  
  beta_gjam<-merge(beta_mean,beta_q25, by=c("beta","sp","model"))
  beta_gjam<-merge(beta_gjam,beta_q97,by=c("beta","sp","model"))
  
  
  ########GJAM_DR###################################################
  gjdr_bchain<-gjam_dr_bchain
  gjdr_bchain<-S$bchain
  beta_mean_raw<- apply(gjdr_bchain,2,mean)
  beta_q25_raw<- apply(gjdr_bchain,2,quantile,0.025)
  beta_q97_raw<- apply(gjdr_bchain,2,quantile,0.975)
  
  beta_mean<- prepare_beta_gjam(beta_mean_raw,"mean", "gjam_dr")
  beta_q25<- prepare_beta_gjam(beta_q25_raw,"low", "gjam_dr")
  beta_q97<- prepare_beta_gjam(beta_q97_raw,"high", "gjam_dr")
  
  beta_gjam_dr<-merge(beta_mean,beta_q25, by=c("beta","sp","model"))
  beta_gjam_dr<-merge(beta_gjam_dr,beta_q97,by=c("beta","sp","model"))
  
  
  
  ####hmsc
  
  beta_hm<-array(NA,dim=c(hmmod$nc,hmmod$ns,2*hmmod$samples)) #2 is nchains
  for(k in 1:(2*hmmod$samples)){
    if(k<=hmmod$samples){beta_hm[,,k] <- hmmod$postList[[1]][[k]]$Beta
    }else {beta_hm[,,k] <- hmmod$postList[[2]][[k-hmmod$samples]]$Beta}
    
  }
  beta_mean_raw<- apply(beta_hm,c(1,2),mean)
  beta_q25_raw<- apply(beta_hm,c(1,2),quantile,0.025)
  beta_q97_raw<- apply(beta_hm,c(1,2),quantile,0.975)
  
  
  beta_mean<- prepare_beta_hmsc(datab=beta_mean_raw,name="mean", class="hmsc",vcol=hmmod$spNames,vrow=hmmod$covNames)
  beta_q25<- prepare_beta_hmsc(beta_q25_raw,"low", "hmsc",vcol=hmmod$spNames,vrow=hmmod$covNames)
  beta_q97<- prepare_beta_hmsc(beta_q97_raw,"high", "hmsc",vcol=hmmod$spNames,vrow=hmmod$covNames)
  
  beta_hmsc<-merge(beta_mean,beta_q25, by=c("beta","sp","model"))
  beta_hmsc<-merge(beta_hmsc,beta_q97,by=c("beta","sp","model"))
  
  
  
  data_beta<-rbind(beta_gjam,beta_gjam_dr,beta_jsdm,beta_hmsc)
  data_beta$name<- "Specie"
  data_beta$name<- paste(data_beta$name,data_beta$sp)
  data_beta[, 'model'] <- as.factor(data_beta[, 'model'])
  
  
  if(nsp==5){
    p <- ggplot(data_beta, aes(mean, beta, colour = model))
    q<- p + geom_point( size = 2) +xlab("Environmental gradient")+ylab("Coefficients")+scale_y_discrete(
      labels=c("inter" =expression(beta[0]),"env" = expression(beta[1]),"env2" = expression(beta[2])))+
      geom_errorbarh(aes(xmax = high, xmin = low, height = .3))+ facet_grid(name ~.,scales = "free")+
      geom_vline(aes(xintercept = 0),size = 0.2, colour = "red") + theme(strip.text.x = element_text(size = 10, colour = "black")) +
      ggtitle("Comparison of the estimated covariates coefficiens")
    q
  }else{
    for(i in 1:nsp){
      tmp<-data_beta[which(data_beta$sp==i),]
      p <- ggplot(tmp, aes(mean, beta, colour = model))
      q<- p + geom_point( size = 2) +xlab("Environmental gradient")+ylab("Coefficients")+scale_y_discrete(
        labels=c("inter" =expression(beta[0]),"env" = expression(beta[1]),"env2" = expression(beta[2])))+
        geom_errorbarh(aes(xmax = high, xmin = low, height = .3))+
        geom_vline(aes(xintercept = 0),size = 0.2, colour = "red") + theme(strip.text.x = element_text(size = 10, colour = "black")) +
        ggtitle(paste("Comparison of the estimated covariates coefficiens for specie ",i))
      print(q)
    }
  }
}





ALL3<-function(jsdm_mod,gjam_mod,hmsc_mod,interact=diag(5)){
  par(mfrow=c(2,2),oma = c(1, 1, 1, 1))
  corrplot(jsdm_mod, diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
  title("R  JSDM ")
  corrplot(gjam_mod, diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
  title("R GJAM")
  corrplot(hmsc_mod, diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
  title("R HMSC")
  corrplot(interact, diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
  title("True interactions")
}



ALL4<-function(jsdm_mod,gjam_mod,hmsc_mod,interact=diag(5)){
  par(mfrow=c(2,2),oma = c(1, 1, 1, 1))
  #corrplot(jsdm_mod, diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
  #title("Partial correlation  JSDM ")
  corrplot(gjam_mod, diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
  title("Partial correlation GJAM")
  corrplot(hmsc_mod, diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
  title("Partial correlation HMSC")
  corrplot(interact, diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
  title("True interactions")
}







gjam_dim_red<-function(datab,ng=2500, burnin=500,r=1, name="./gjam_models/gjamDR5env.rda", regime="F"){
  ydata<-subset(datab, select = -env)
  xdata<-scale(poly(datab$env, 2))
  ns<- ncol(ydata)
  n<-nrow(datab)
  N<-ncol(ydata)-1
  colnames(xdata)<- c("env","env2")
  formula<- ~env + env2
  if (r==1){
    rseq<-seq(1,N,by=2)
    print(rseq)
    DIC_array<-array(NA,dim=length(rseq)-2)
    RSMPE_array<-array(NA,dim=length(rseq)-2)
    for (j in (2:(length(rseq)-1))){
      rl <-list(N=N, r=rseq[j])
      ml  <- list(ng = ng, burnin = burnin, typeNames = 'PA',reductList=rl, PREDICTX = F )
      mod_gjam_c <- gjam(formula, xdata, ydata, modelList = ml)
      DIC_array[j-1]<-mod_gjam_c$fit$DIC
      RSMPE_array[j-1]<-mod_gjam_c$fit$rmspeAll
    }
    r<- rseq[which.min(RSMPE_array)+1]
    print(DIC_array)
  }
  
  rl <-list(N=N, r=r)
  ml  <- list(ng = ng, burnin = burnin, typeNames = 'PA',reductList=rl)
  if(regime=="L"){
    gjam_dr_mods<-load_object(name)
    mod_gjam_red_1<-gjam_dr_mods$m1
    mod_gjam_red_2<-gjam_dr_mods$m2
  }else{
    mod_gjam_red_1 <- gjam(formula, xdata, ydata, modelList = ml)
    mod_gjam_red_2 <- gjam(formula, xdata, ydata, modelList = ml)
    gjam_dr_mods<- list(m1=mod_gjam_red_1,m2=mod_gjam_red_2)
    save(gjam_dr_mods, file = name)
  }
  gjam_bs<- mcmc.list(mcmc(mod_gjam_red_1$chains$bgibbsUn[-(1:burnin),]),mcmc(mod_gjam_red_2$chains$bgibbsUn[-(1:burnin),]))
  gjam_sigma<- mcmc.list(mcmc(mod_gjam_red_1$chains$sgibbs[-(1:burnin),]),mcmc(mod_gjam_red_2$chains$sgibbs[-(1:burnin),]))
  n_eff_beta<- as.data.frame(effectiveSize(gjam_bs))
  colnames(n_eff_beta)<- c("value")
  n_eff_beta$parameter<- "beta"
  n_eff_sigma<- as.data.frame(effectiveSize(gjam_sigma))
  colnames(n_eff_sigma)<- c("value")
  n_eff_sigma$parameter<- "correlation"
  neff<- rbind(n_eff_beta,n_eff_sigma)
  p<- ggplot(neff, aes(x=value, color=parameter,fill=parameter)) +
    geom_histogram( alpha=0.4, position="identity",binwidth = 50) + xlab("effective size") +
    scale_color_brewer(palette="Dark2")+scale_fill_brewer(palette="Dark2")+
    ggtitle("Effective size for the parameters beta and coefficients of correlation matrix")
  plot(p)
  
  Rhat_beta<-as.data.frame(gelman.diag(gjam_bs, multivariate=FALSE)$psrf)
  Rhat_beta$parameter<- "beta"
  Rhat_sigma<- as.data.frame(gelman.diag(gjam_sigma, multivariate=FALSE)$psrf)
  Rhat_sigma$parameter<- "correlation"
  Rhat<- rbind(Rhat_beta,Rhat_sigma)
  p2<- ggplot(Rhat, aes(x= Rhat$`Point est.`, color=parameter,fill=parameter)) +
    geom_histogram( alpha=0.4, position="identity",binwidth =0.1) +
    scale_color_brewer(palette="Dark2")+scale_fill_brewer(palette="Dark2")+ xlab("Rhat") +
    ggtitle("Rhat for the parameters beta and the coefficients of correlation matrix")
  plot(p2)
  
  #gelman.plot(gjam_bs)
  #xyplot(gjam_bs)
  #traceplot(gjam_sigma)
  # hist(effectiveSize(gjam_sigma), main="ess(sigma)",lwd=2,col=gray(.6),breaks = 50)
  # hist(gelman.diag(gjam_sigma,multivariate=FALSE)$psrf,lwd=2,col=gray(.6), main="psrf(sigma)",breaks = 50)
  #gelman.plot(gjam_sigma)
  #xyplot(gjam_sigma)
  sgibbs<-abind(mod_gjam_red_1$chains$sgibbs[-(1:burnin),],mod_gjam_red_2$chains$sgibbs[-(1:burnin),],along=1)
  sigErrGibbs<-abind(mod_gjam_red_1$chains$sigErrGibbs[-(1:burnin)],mod_gjam_red_2$chains$sigErrGibbs[-(1:burnin)],along=1)
  kgibbs<-abind(mod_gjam_red_1$chains$kgibbs[-(1:burnin),],mod_gjam_red_2$chains$kgibbs[-(1:burnin),],along=1)
  sigma<-invsigma<-array(NA,dim=c(ns,ns,2*(ng-burnin)))
  #.expandSigmaChains(snames, sgibbs, otherpar, simIndex=simIndex, sigErrGibbs, kgibbs, REDUCT)
  sigma<-invsigma<-array(NA,dim=c(ns,ns,2*(ng-burnin)))
  
  # sgibbs<-mod_gjam_red$chains$sgibbs
  # sigErrGibbs<-mod_gjam_red$chains$sigErrGibbs
  # kgibbs<-mod_gjam_red$chains$kgibbs
  N<-mod_gjam_red_1$modelList$reductList$N
  r<-mod_gjam_red_1$modelList$reductList$r
  N_dim<-2*(ng-burnin)
  for(j in 1:N_dim){
    Z  <- matrix(sgibbs[j,],N,r)
    sigma[,,j] <- .expandSigma(sigErrGibbs[j], ns, Z = Z, kgibbs[j,], REDUCT = T) #sigma
    invsigma[,,j] <- invWbyRcpp(sigErrGibbs[j], Z[kgibbs[j,],]) #inverse sigma
  } 
  sigma_mean<-apply(sigma,c(1,2),mean) 
  sigma_q05<-apply(sigma,c(1,2),quantile,0.05) 
  sigma_q95<-apply(sigma,c(1,2),quantile,0.95) 
  Sigma_sign<--cov2cor(sigma_mean*(!(sigma_q95>0 & sigma_q05<0)))
  
  invsigma_mean<-apply(invsigma,c(1,2),mean) 
  invsigma_q05<-apply(invsigma,c(1,2),quantile,0.05) 
  invsigma_q95<-apply(invsigma,c(1,2),quantile,0.95) 
  INVSigma_sign<--cov2cor(sigma_mean*(!(invsigma_q95>0 & invsigma_q05<0)))
  x <- cbind(1, scale(poly(datab$env, 2)))
  J<- ncol(datab) - 1
  bgibbspr<-abind(mod_gjam_red_1$chains$bgibbs[-(1:burnin),],mod_gjam_red_2$chains$bgibbs[-(1:burnin),],along=1)
  N_dim<-2*(ng-burnin)
  mu<-array(NA,dim=c(n,J,N_dim))
  for(k in 1:N_dim){
    for(j in 1:J){
      
      mu[,j,k] <- pnorm(x%*%bgibbspr[k,(3*(j-1)+1):(3*j)])
      
    }
  }
  
  return(list(Rho_sign_d=Sigma_sign,Tau_sign_d=INVSigma_sign,k=kgibbs, predict=mu, bchain=bgibbspr))
}


###  Common functions##############################################################################################################################################



###  Run analysis 5##############################################################################################################################################





analysis_chunk_5<- function(datab,gjam, hmsc, jsdm, gjam_dr){
  ###jsdm######################################################################
  jsdm_mod <- load_object(jsdm)
  summary(jsdm_mod)
  jsdm_mod$mcmc.info[1:7]
  j_metric_FacCompSparse20<-metrics_jsdm(cmp_fds20,comp = comp_inter[[21]],fac=fac_inter[[21]],only_env = F)
  # ######################################################Prepare data
  data_mod <- list(
    Y = subset(datab, select = -env),
    X = cbind(1, scale(poly(datab$env, 2))),
    covx = cov(cbind(1, scale(poly(datab$env, 2)))),
    K = 3,
    J = ncol(datab) - 1,
    n = nrow(datab),
    I = diag(ncol(datab) - 1),
    df = ncol(datab)
  )
  mod_list_Rho<-list()
  mod_list_Rho<-list(jsdm =jsdm_mod$mean$Rho*(!jsdm_mod$overlap0$Rho)) 
  mod_list_Tau<-list()
  mod_list_Tau<-list(jsdm =jsdm_mod$mean$Rho*(!jsdm_mod$overlap0$Rho)) 
  
  pred_j<-array(NA,dim=c(data_mod$n,data_mod$J,jsdm_mod$mcmc.info$n.samples))
  for(i in 1:jsdm_mod$mcmc.info$n.samples){
    pred_j[,,i]<-pnorm(data_mod$X%*%t(jsdm_mod$sims.list$B[i,,]))
  }
  pred_j_mean <- apply(pred_j, 1:2, mean)
  pred_j_05 <- apply(pred_j, 1:2, quantile,0.05)
  pred_j_95 <- apply(pred_j, 1:2, quantile,0.95)
  pred_jsdm<-list(pred_j_mean=pred_j_mean,pred_j_05=pred_j_05,pred_j_95=pred_j_95)
  for(i in 1:data_mod$J) AUC_j_comp_fac_sparse<-c(AUC_j_comp_fac_sparse,auc(roc(pred_j_mean[,i],factor(data_mod$Y[,i]))))
  #####gjam#################################################################
  gjam_mod<-fit_gjam(datab,10000,1000,gjam,interact= (-1)*comp_inter[[21]]+fac_inter[[21]])
  gjam_mod<-load_gjam(datab,10000,1000,name=gjam, interact= (-1)*comp_inter[[21]]+fac_inter[[21]])
  g_metric_FacCompSparse20<-metrics_gjam(gjam_mod$Rho_sign,comp=comp_inter[[21]], fac=fac_inter[[21]],only_env = F)
  g_metric_FacCompSparse20_p<-metrics_gjam(gjam_mod$Tau_sign,comp=comp_inter[[21]], fac=fac_inter[[21]],only_env = F)
  mod_list_Rho<-list.append(mod_list_Rho, gjam=gjam_mod$Rho_sign) 
  mod_list_Tau<-list.append(mod_list_Tau, gjam=gjam_mod$Tau_sign) 
  pred_gj_mean <-apply(gjam_mod$predict, 1:2,mean)
  pred_gj_05 <- apply(gjam_mod$predict, 1:2,quantile,0.05)
  pred_gj_95 <- apply(gjam_mod$predict, 1:2,quantile,0.95)
  pred_gjam<- list(pred_gj_mean=pred_gj_mean,pred_gj_05=pred_gj_05,pred_gj_95=pred_gj_95)
  for(i in 1:(ncol(data)-1)) AUC_g_comp_fac_sparse<-c(AUC_g_comp_fac_sparse,auc(roc(pred_gj_mean[,i],factor(data[,i]))))
  #####gjam_dr############################################################
  gjam_dr_mod<- gjam_dim_red_5(datab,ng=10000, burnin=1000,name=gjam_dr, regime="L")
  g_dr_metric_FacCompSparse5<-metrics_gjam(gjam_dr_mod$Rho_sign_d)
  g_dr_metric_FacCompSparse5_p<-metrics_gjam(gjam_dr_mod$Tau_sign_d)
  AUC_gdr_comp_fac_sparse<-vector()
  ##Prediction
  pred_gj_dr_mean <-apply(gjam_dr_mod$predict, 1:2,mean)
  pred_gj_dr_05 <- apply(gjam_dr_mod$predict, 1:2,quantile,0.05)
  pred_gj_dr_95 <- apply(gjam_dr_mod$predict, 1:2,quantile,0.95)
  for(i in 1:(ncol(datab)-1)) AUC_gdr_comp_fac_sparse<-c(AUC_gdr_comp_fac_sparse,auc(roc(pred_gj_dr_mean[,i],factor(datab[,i]))))
  pred_gjam_dr<- list(pred_gj_dr_mean=pred_gj_dr_mean,pred_gj_dr_05=pred_gj_dr_05,pred_gj_dr_95=pred_gj_dr_95)
  #####hmsc###############################################################
  hm_mod<-fit_hmsc(datab,"Fit",nsamples=10000, nchains=2,name=hmsc )
  hm_mod<-fit_hmsc(datab,"Load",nsamples=10000, nchains=2,name=hmsc )
  hm_conv(hm_mod)
  hm_mod_R<-hm_inter(hm_mod,nsamples=10000, nchains=2,interact =  (-1)*comp_inter[[21]] +fac_inter[[21]])
  h_metric_FacCompSparse20<-metrics_hmsc(hm_mod_R$Rho_sign,comp=comp_inter[[21]],fac=fac_inter[[21]],only_env = F)
  h_metric_FacCompSparse20_p<-metrics_hmsc(hm_mod_R$Tau_sign,comp=comp_inter[[21]],fac=fac_inter[[21]],only_env = F)
  mod_list_Rho<-list.append(mod_list_Rho,hmsc=hm_mod_R$Rho_sign)
  mod_list_Tau<-list.append(mod_list_Tau,hmsc=hm_mod_R$Tau_sign)
  ####Prediction
  x = cbind(1, scale(poly(datab$env, 2)))
  pred_hm<-array(NA,dim=c(hm_mod$ny,hm_mod$ns,2*hm_mod$samples)) #2 is nchains
  for(k in 1:(2*hm_mod$samples)){
    if(k<=hm_mod$samples){pred_hm[,,k] <- pnorm(x%*%hm_mod$postList[[1]][[k]]$Beta)
    }else {pred_hm[,,k] <- pnorm(x%*%hm_mod$postList[[2]][[k-hm_mod$samples]]$Beta)}
    
  }
  pred_hmsc<- list(pred_hm_mean=apply(pred_hm, c(1,2), mean),
                   pred_hm_05=apply(pred_hm, c(1,2), quantile,0.05),
                   pred_hm_95=apply(pred_hm, c(1,2), quantile,0.95))
  
  for(i in 1:(ncol(datab)-1)) AUC_h_comp_fac_sparse<-c(AUC_h_comp_fac_sparse,auc(roc(pred_hmsc$pred_hm_mean[,i],factor(datab[,i]))))
}


analysis_chunk<- function(datab,gjam, hmsc, jsdm, gjam_dr){
  ###jsdm######################################################################
  jsdm_mod <- load_object(jsdm)
  summary(jsdm_mod)
  jsdm_mod$mcmc.info[1:7]
  j_metric<-metrics_jsdm(jsdm_mod,comp = comp_inter[[21]],fac=fac_inter[[21]],only_env = F)
  # ######################################################Prepare data
  data_mod <- list(
    Y = subset(datab, select = -env),
    X = cbind(1, scale(poly(datab$env, 2))),
    covx = cov(cbind(1, scale(poly(datab$env, 2)))),
    K = 3,
    J = ncol(datab) - 1,
    n = nrow(datab),
    I = diag(ncol(datab) - 1),
    df = ncol(datab)
  )
  mod_list_Rho<-list()
  mod_list_Rho<-list(jsdm =jsdm_mod$mean$Rho*(!jsdm_mod$overlap0$Rho)) 
  mod_list_Tau<-list()
  mod_list_Tau<-list(jsdm =jsdm_mod$mean$Rho*(!jsdm_mod$overlap0$Rho)) 
  
  pred_j<-array(NA,dim=c(data_mod$n,data_mod$J,jsdm_mod$mcmc.info$n.samples))
  for(i in 1:jsdm_mod$mcmc.info$n.samples){
    pred_j[,,i]<-pnorm(data_mod$X%*%t(jsdm_mod$sims.list$B[i,,]))
  }
  pred_j_mean <- apply(pred_j, 1:2, mean)
  pred_j_05 <- apply(pred_j, 1:2, quantile,0.05)
  pred_j_95 <- apply(pred_j, 1:2, quantile,0.95)
  pred_jsdm<-list(pred_j_mean=pred_j_mean,pred_j_05=pred_j_05,pred_j_95=pred_j_95)
  for(i in 1:data_mod$J) AUC_j<-auc(roc(pred_j_mean[,i],factor(data_mod$Y[,i])))
  #####gjam#################################################################
  gjam_mod<-fit_gjam(datab,10000,1000,gjam,interact= (-1)*comp_inter[[21]]+fac_inter[[21]])
  gjam_mod<-load_gjam(datab,10000,1000,name=gjam, interact= (-1)*comp_inter[[21]]+fac_inter[[21]])
  g_metric<-metrics_gjam(gjam_mod$Rho_sign,comp=comp_inter[[21]], fac=fac_inter[[21]],only_env = F)
  g_metric_p<-metrics_gjam(gjam_mod$Tau_sign,comp=comp_inter[[21]], fac=fac_inter[[21]],only_env = F)
  mod_list_Rho<-list.append(mod_list_Rho, gjam=gjam_mod$Rho_sign) 
  mod_list_Tau<-list.append(mod_list_Tau, gjam=gjam_mod$Tau_sign) 
  pred_gj_mean <-apply(gjam_mod$predict, 1:2,mean)
  pred_gj_05 <- apply(gjam_mod$predict, 1:2,quantile,0.05)
  pred_gj_95 <- apply(gjam_mod$predict, 1:2,quantile,0.95)
  pred_gjam<- list(pred_gj_mean=pred_gj_mean,pred_gj_05=pred_gj_05,pred_gj_95=pred_gj_95)
  for(i in 1:(ncol(data)-1)) AUC_g<-auc(roc(pred_gj_mean[,i],factor(data[,i])))
  #####gjam_dr############################################################
  gjam_dr_mod<- gjam_dim_red(datab,ng=10000, burnin=1000,r=3,name=gjam_dr, regime="F")
  g_dr_metric<-metrics_gjam(gjam_dr_mod$Rho_sign_d)
  g_dr_metric_p<-metrics_gjam(gjam_dr_mod$Tau_sign_d)
  
  ##Prediction
  pred_gj_dr_mean <-apply(gjam_dr_mod$predict, 1:2,mean)
  pred_gj_dr_05 <- apply(gjam_dr_mod$predict, 1:2,quantile,0.05)
  pred_gj_dr_95 <- apply(gjam_dr_mod$predict, 1:2,quantile,0.95)
  for(i in 1:(ncol(datab)-1)) AUC_gdr<-auc(roc(pred_gj_dr_mean[,i],factor(datab[,i])))
  pred_gjam_dr<- list(pred_gj_dr_mean=pred_gj_dr_mean,pred_gj_dr_05=pred_gj_dr_05,pred_gj_dr_95=pred_gj_dr_95)
  #####hmsc###############################################################
  hm_mod<-fit_hmsc(datab,"Fit",nsamples=10000, nchains=2,name=hmsc )
  hm_mod<-fit_hmsc(datab,"Load",nsamples=10000, nchains=2,name=hmsc )
  hm_conv(hm_mod)
  hm_mod_R<-hm_inter(hm_mod,nsamples=10000, nchains=2,interact =  (-1)*comp_inter[[21]] +fac_inter[[21]])
  h_metric<-metrics_hmsc(hm_mod_R$Rho_sign,comp=comp_inter[[21]],fac=fac_inter[[21]],only_env = F)
  h_metric_p<-metrics_hmsc(hm_mod_R$Tau_sign,comp=comp_inter[[21]],fac=fac_inter[[21]],only_env = F)
  mod_list_Rho<-list.append(mod_list_Rho,hmsc=hm_mod_R$Rho_sign)
  mod_list_Tau<-list.append(mod_list_Tau,hmsc=hm_mod_R$Tau_sign)
  ####Prediction
  x = cbind(1, scale(poly(datab$env, 2)))
  pred_hm<-array(NA,dim=c(hm_mod$ny,hm_mod$ns,2*hm_mod$samples)) #2 is nchains
  for(k in 1:(2*hm_mod$samples)){
    if(k<=hm_mod$samples){pred_hm[,,k] <- pnorm(x%*%hm_mod$postList[[1]][[k]]$Beta)
    }else {pred_hm[,,k] <- pnorm(x%*%hm_mod$postList[[2]][[k-hm_mod$samples]]$Beta)}
    
  }
  pred_hmsc<- list(pred_hm_mean=apply(pred_hm, c(1,2), mean),
                   pred_hm_05=apply(pred_hm, c(1,2), quantile,0.05),
                   pred_hm_95=apply(pred_hm, c(1,2), quantile,0.95))
  
  for(i in 1:(ncol(datab)-1)) AUC_h<-auc(roc(pred_hmsc$pred_hm_mean[,i],factor(datab[,i])))
  return(list(AUC_l=list(AUC_h,AUC_gdr,AUC_g,AUC_j), metric=list(h_metric,g_dr_metric,g_metric,j_metric), metric_p=list(h_metric_p,g_dr_metric_p,g_metric_p),
              Rho=mod_list_Rho,Tau=mod_list_Tau, prediction=list(pred_hmsc,pred_gjam_dr,pred_gjam,pred_jsdm))
}


################################################################################################################################
################################################################################################################################
################################################################################################################################

## Environmental filtering + Competition +Facilitation 20 sparse species

### JSDM 

cmp_fds20 <- load_object("model-2019-04-23-09-42-18.rda")
jsdm_conv(cmp_fds20)
cmp_fds20$mcmc.info[1:7]
# 
j_metric_FacCompSparse20<-metrics_jsdm(cmp_fds20,comp = comp_inter[[21]],fac=fac_inter[[21]],only_env = F)
# ######################################################Prepare data
data<-sim_data$FacCompSparseSp20
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

mod_list_Rho<-list()
mod_list_Rho<-list(jsdm =cmp_fds20$mean$Rho*(!cmp_fds20$overlap0$Rho)) 

mod_list_Tau<-list()
mod_list_Tau<-list(jsdm =cmp_fds20$mean$Rho*(!cmp_fds20$overlap0$Rho)) 

pred_j<-array(NA,dim=c(data$n,data$J,cmp_fds20$mcmc.info$n.samples))
for(i in 1:cmp_fds20$mcmc.info$n.samples){
  
  pred_j[,,i]<-pnorm(data$X%*%t(cmp_fds20$sims.list$B[i,,]))
  
}

pred_j_mean <- apply(pred_j, 1:2, mean)
pred_j_05 <- apply(pred_j, 1:2, quantile,0.05)
pred_j_95 <- apply(pred_j, 1:2, quantile,0.95)

pred_jsdm<-list(pred_j_mean=pred_j_mean,pred_j_05=pred_j_05,pred_j_95=pred_j_95)

for(i in 1:data$J) AUC_j_comp_fac_sparse<-c(AUC_j_comp_fac_sparse,auc(roc(pred_j_mean[,i],factor(data$Y[,i]))))





### Gjam

data<-sim_data$FacCompSparseSp20
#gjam_mod<-fit_gjam(data,10000,1000,"./gjam_models/gjamcmp_facs20.rda",interact= (-1)*comp_inter[[21]]+fac_inter[[21]])
gjam_mod<-load_gjam(data,10000,1000,name="./gjam_models/gjamcmp_facs20.rda", interact= (-1)*comp_inter[[21]]+fac_inter[[21]])
g_metric_FacCompSparse20<-metrics_gjam(gjam_mod$Rho_sign,comp=comp_inter[[21]], fac=fac_inter[[21]],only_env = F)
g_metric_FacCompSparse20_p<-metrics_gjam(gjam_mod$Tau_sign,comp=comp_inter[[21]], fac=fac_inter[[21]],only_env = F)

cat("\n")
cat(sprintf("Success rate for non-interacting: %s\n", round(g_metric_FacCompSparse20$success_env,3)))
cat(sprintf("Success rate for competition: %s\n", round(g_metric_FacCompSparse20$success_comp,3)))
cat(sprintf("Success rate for facilitation: %s\n", round(g_metric_FacCompSparse20$success_fac,3)))
mod_list_Rho<-list.append(mod_list_Rho, gjam=gjam_mod$Rho_sign) 
mod_list_Tau<-list.append(mod_list_Tau, gjam=gjam_mod$Tau_sign) 


pred_gj_mean <-apply(gjam_mod$predict, 1:2,mean)
pred_gj_05 <- apply(gjam_mod$predict, 1:2,quantile,0.05)
pred_gj_95 <- apply(gjam_mod$predict, 1:2,quantile,0.95)

pred_gjam<- list(pred_gj_mean=pred_gj_mean,pred_gj_05=pred_gj_05,pred_gj_95=pred_gj_95)

for(i in 1:(ncol(data)-1)) AUC_g_comp_fac_sparse<-c(AUC_g_comp_fac_sparse,auc(roc(pred_gj_mean[,i],factor(data[,i]))))




### HMSC


data<-sim_data$FacCompSparseSp20
hm_mod<-fit_hmsc(data,"Load",nsamples=2000, nchains=2,name="./HMmodels/hmcmp_facs20.rda" )
hm_conv(hm_mod)
hm_mod_R<-hm_inter(hm_mod, nsamples=2000, nchains=2,interact =  (-1)*comp_inter[[21]] +fac_inter[[21]])
h_metric_FacCompSparse20<-metrics_hmsc(hm_mod_R$Rho_sign,comp=comp_inter[[21]],fac=fac_inter[[21]],only_env = F)
h_metric_FacCompSparse20_p<-metrics_hmsc(hm_mod_R$Tau_sign,comp=comp_inter[[21]],fac=fac_inter[[21]],only_env = F)

cat("\n")
cat(sprintf("Success rate for non-interacting: %s\n", round(h_metric_FacCompSparse20$success_env,3)))
cat(sprintf("Success rate for competition: %s\n", round(h_metric_FacCompSparse20$success_comp,3)))
cat(sprintf("Success rate for facilitation: %s\n", round(h_metric_FacCompSparse20$success_fac,3)))
#mod_list<-list.append(mod_list,hmsc=Rho_sign)
mod_list_Rho<-list.append(mod_list_Rho,hmsc=hm_mod_R$Rho_sign)
mod_list_Tau<-list.append(mod_list_Tau,hmsc=hm_mod_R$Tau_sign)
####Prediction
x = cbind(1, scale(poly(data$env, 2)))
pred_hm<-array(NA,dim=c(hm_mod$ny,hm_mod$ns,2*hm_mod$samples)) #2 is nchains
for(k in 1:(2*hm_mod$samples)){
  if(k<=hm_mod$samples){pred_hm[,,k] <- pnorm(x%*%hm_mod$postList[[1]][[k]]$Beta)
  }else {pred_hm[,,k] <- pnorm(x%*%hm_mod$postList[[2]][[k-hm_mod$samples]]$Beta)}
  
}
pred_hmsc<- list(pred_hm_mean=apply(pred_hm, c(1,2), mean),
                 pred_hm_05=apply(pred_hm, c(1,2), quantile,0.05),
                 pred_hm_95=apply(pred_hm, c(1,2), quantile,0.95))

for(i in 1:(ncol(data)-1)) AUC_h_comp_fac_sparse<-c(AUC_h_comp_fac_sparse,auc(roc(pred_hmsc$pred_hm_mean[,i],factor(data[,i]))))




### Dimension reduction

################################################################################################################################

#r-opt=15
S<- gjam_dim_red(sim_data$FacCompSparseSp20,ng=10000, burnin=1000,r=3,name="./gjam_models/gjamDR20compfacs.rda", regime="L")
g_dr_metric_FacCompSparse20<-metrics_gjam(S$Rho_sign_d)
g_dr_metric_FacCompSparse20_p<-metrics_gjam(S$Tau_sign_d)


##Prediction
pred_gj_dr_mean <-apply(S$predict, 1:2,mean)
pred_gj_dr_05 <- apply(S$predict, 1:2,quantile,0.05)
pred_gj_dr_95 <- apply(S$predict, 1:2,quantile,0.95)
data<-sim_data$FacCompSparseSp20
for(i in 1:(ncol(data)-1)) AUC_gdr_comp_fac_sparse<-c(AUC_gdr_comp_fac_sparse,auc(roc(pred_gj_dr_mean[,i],factor(data[,i]))))
pred_gjam_dr<- list(pred_gj_dr_mean=pred_gj_dr_mean,pred_gj_dr_05=pred_gj_dr_05,pred_gj_dr_95=pred_gj_dr_95)


par(mfrow=c(2,2),oma = c(1, 1, 1, 1))
corrplot(S$Rho_sign_d, diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
title("Correlation from dimension reduction")
corrplot(S$Tau_sign_d, diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
title("Partial Correlation from dimension reduction")
corrplot(comp.psm(S$k), diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=gcols(200),cl.lim=c(0, 1), type = "lower")
title("Partial Correlation from dimension reduction")
corrplot((-1)*comp_inter[[21]]+fac_inter[[21]], diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
title("True interactions")


################################################################################################################################


### Summary plots
create_beta_plot(load_object("./gjam_models/gjamcmp_facs20.rda"),cmp_fds20,load_object("./HMmodels/hmcmp_facs20.rda"), S$bchain)
create_plot_response(sim_data$FacCompSparseSp20, nsp=20,pred_gjam, pred_jsdm,pred_hmsc,pred_gjam_dr )
ALL3(mod_list_Rho$jsdm,mod_list_Rho$gjam,mod_list_Rho$hmsc, (-1)*comp_inter[[21]]+ fac_inter[[21]])
ALL4(mod_list_Tau$jsdm,mod_list_Tau$gjam,mod_list_Tau$hmsc, (-1)*comp_inter[[21]]+ fac_inter[[21]])
################################################################################################################################






