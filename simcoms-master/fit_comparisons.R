##SET UP ##########################################################################################################
rm(list=ls())
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
library(LearnBayes)
Rcpp::sourceCpp('/Users/dariabystrova/Documents/GitHub/gjamed/src/cppFns.cpp')
source("/Users/dariabystrova/Documents/GitHub/gjamed/R/gjamHfunctions_mod.R")
#setwd("C:/Users/giaru/Desktop/Documents/GitHub/Ecology-models/simcoms-master")
####################################Set up###############################################
##Set the directory
setwd("/Users/dariabystrova/Documents/GitHub/Ecology-models/simcoms-master/ExampleFiles")
#setwd("~/Tesi/Code/Ecology-models-master_1/simcoms-master")

# lapply(list.files(path = "."),load,.GlobalEnv)
#setwd("~/Tesi/Code/Ecology-models-master/simcoms-master")
load("params.rds")
load("sim_names.rds")
load("comp_inter.rds")
load("fac_inter.rds")

setwd("~/Documents/GitHub/Ecology-models/simcoms-master")
sim_data<-readRDS("sim_data.rds")


#####################################Functions###########################################

load_object <- function(file) {
  tmp <- new.env()
  load(file = file, envir = tmp)
  tmp[[ls(tmp)[1]]]
}


############Functions JSDM###################################################################################################

jsdm_conv<-function(mod) {
  n_eff_beta<- as.data.frame(t(mod$n.eff$B))
  n_eff_beta$parameter<- "beta"
  n_eff_sigma<- as.data.frame(mod$n.eff$Rho)
  n_eff_sigma$parameter<- "correlation"
  neff<- rbind(n_eff_beta,n_eff_sigma)
  neff_mod<- melt(neff, id=c("parameter"), variable_name = "effective size")
  Rhat_beta<- as.data.frame(t(mod$Rhat$B))
  Rhat_beta$parameter<- "beta"
  Rhat_sigma<- as.data.frame(mod$Rhat$Rho)
  Rhat_sigma$parameter<- "correlation"
  Rhat<- rbind(Rhat_beta,Rhat_sigma)
  Rhat_mod<- melt(Rhat, id=c("parameter"), variable_name = "Rhat")
  Rhat_mod_nonNA<- Rhat_mod[!is.na(Rhat_mod$value),]
  return(list(neff=neff_mod,Rhat=Rhat_mod_nonNA))
}


metrics_jsdm<-function(model,fac=NULL,comp=NULL,only_env=T){
  
  if(is.null(dim(fac))&is.null(dim(comp))){
    only_env=T;
    fac=NULL
    comp=NULL
  }else{
    if(is.null(dim(fac))){
      only_env=F
      fac=NULL
    }else{
      
      only_env=F
      comp=NULL
    }
    
    
  }
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

jsdm_tau<-function(jsdm_mod){
  nsamples<-jsdm_mod$mcmc.info$n.samples
  ns<- ncol(jsdm_mod$model$cluster1$data()$Y)
  postT<-array(dim=c(ns,ns,nsamples))
  for(i in 1:nsamples){
    postT[,,i]<- -stats::cov2cor(solve(jsdm_mod$sims.list$Rho[i,,]+0.1*diag(ns)))
  }
  postTMean = apply(postT,c(1,2),mean)
  postTUp=apply(postT,c(1,2),quantile,0.95)
  postTLo=apply(postT,c(1,2),quantile,0.05)
  
  Toplot_T<-postTMean*(!(postTUp>0 & postTLo<0))
  return(list(Tau=postTMean,Tau_sign=Toplot_T)) 
}

jsdm_rho<-function(jsdm_mod){
  return(list(Rho=jsdm_mod$mean$Rho,Rho_sign= (!jsdm_mod$overlap0$Rho)*jsdm_mod$mean$Rho))
}

dataprep<- function(datap){
  data <- list(
    Y = subset(datap, select = -env),
    X = cbind(1, scale(poly(datap$env, 2))),
    covx = cov(cbind(1, scale(poly(datap$env, 2)))),
    K = 3,
    J = ncol(datap) - 1,
    n = nrow(datap),
    I = diag(ncol(datap) - 1),
    df = ncol(datap)
  )
  return(data)
}
predict_cm<- function(datab, mod){
  data<- dataprep(datab)
  pred_j<-array(NA,dim=c(data$n,data$J,mod$mcmc.info$n.samples))
  for(i in 1:mod$mcmc.info$n.samples){
    pred_j[,,i]<-pnorm(data$X%*%t(mod$sims.list$B[i,,]))
  }
  pred_j_mean <- apply(pred_j, 1:2, mean)
  pred_j_05 <- apply(pred_j, 1:2, quantile,0.05)
  pred_j_95 <- apply(pred_j, 1:2, quantile,0.95)
  pred_jsdm<-list(pred_j_mean=pred_j_mean,pred_j_05=pred_j_05,pred_j_95=pred_j_95)
  return(pred_jsdm)
}


############Functions GJAM###################################################################################################

gj_conv<-function(name){
  gj_mod<-load_object(name)
  burn<-gj_mod$m1$modelList$burnin
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
  Rhat_beta<-as.data.frame(gelman.diag(gjam_bs, multivariate=FALSE)$psrf)
  Rhat_beta$parameter<- "beta"
  Rhat_sigma<- as.data.frame(gelman.diag(gjam_sigma, multivariate=FALSE)$psrf)
  Rhat_sigma$parameter<- "correlation"
  Rhat<- rbind(Rhat_beta,Rhat_sigma)
  
  return(list(neff=neff,Rhat=Rhat ))
}


convert_to_m<-function(ar){
  d <-floor((sqrt(length(ar)*8+1)-1)/2)
  C <- matrix(0,d,d)
  i.lwr <- which(lower.tri(C, diag = TRUE), arr.ind=TRUE)
  C[i.lwr] <- ar
  C<-makeSymm(C)
  return(t(C))
}

###this function only loads the model and return the R and T significant
load_gjam<-function(datab,name="./gjam_models/temp.rda"){
  #setup parameters
  data <- dataprep(datab)
  xdata<-as.data.frame(data$X[,-1])
  gj_mod<-load_object(name)
  burn <-  gj_mod$m1$modelList$burnin
  it<- gj_mod$m1$modelList$ng
  S<-ncol(data$Y)
  gjam_bs<- mcmc.list(mcmc(gj_mod$m1$chains$bgibbsUn[-(1:burn),]),mcmc(gj_mod$m2$chains$bgibbsUn[-(1:burn),]))
  gjam_sigma<- mcmc.list(mcmc(gj_mod$m1$chains$sgibbs[-(1:burn),]),mcmc(gj_mod$m2$chains$sgibbs[-(1:burn),]))

  gjam_mods_2bgibbs<-  abind(gj_mod$m1$chains$bgibbsUn,gj_mod$m2$chains$bgibbsUn, along = 3)
  gjam_mods_2sgibbs<-abind(gj_mod$m1$chains$sgibbs,gj_mod$m2$chains$sgibbs, along = 3)
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
  mu<-array(NA,dim=c(data$n,data$J,it))
  for(k in 1:it){
    for(j in 1:data$J){
      mu[,j,k] <- pnorm(x%*%gj_mod$m1$chains$bgibbs[k,(3*(j-1)+1):(3*j)])
    }
  }
  
  
  pred_gj_mean <-apply(mu, 1:2,mean)
  pred_gj_05 <- apply(mu, 1:2,quantile,0.05)
  pred_gj_95 <- apply(mu, 1:2,quantile,0.95)
  
  pred_gjam<- list(pred_gj_mean=pred_gj_mean,pred_gj_05=pred_gj_05,pred_gj_95=pred_gj_95)

  return(list(Rho_sign=R_sign,Tau_sign=Tau_sign,predict=pred_gjam))
}



############Functions HMSC###################################################################################################
hm_conv<-function(name){
  mod<-load_object(name)
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
  Rhat_beta<-as.data.frame(gelman.diag(codaList$Beta,multivariate=FALSE)$psrf)
  Rhat_beta$parameter<- "beta"
  Rhat_sigma<- as.data.frame(gelman.diag(codaList$Omega[[1]], multivariate=FALSE)$psrf)
  Rhat_sigma$parameter<- "correlation"
  Rhat<- rbind(Rhat_beta,Rhat_sigma)
  return(list(neff=neff,Rhat=Rhat))
}

hm_inter<-function(name, nchains=2,nsamples = 1000, interact=diag(ns)){
  mod<- load_object(name)
  getOmega = function(a,r=1)
    return(crossprod(a$Lambda[[r]]))
  ns<-mod$ns
  nsamples<-mod$samples
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

hm_pred<- function(data, name){
  hm_mod<- load_object(name)
  x = cbind(1, scale(poly(data$env, 2)))
  pred_hm<-array(NA,dim=c(hm_mod$ny,hm_mod$ns,2*hm_mod$samples)) #2 is nchains
  for(k in 1:(2*hm_mod$samples)){
    if(k<=hm_mod$samples){pred_hm[,,k] <- pnorm(x%*%hm_mod$postList[[1]][[k]]$Beta)
    }else {pred_hm[,,k] <- pnorm(x%*%hm_mod$postList[[2]][[k-hm_mod$samples]]$Beta)}
  }
  pred_hm_mean <- apply(pred_hm, c(1,2), mean)
  pred_hm_05 <- apply(pred_hm, c(1,2), quantile,0.05)
  pred_hm_95 <- apply(pred_hm, c(1,2), quantile,0.95)
  pred_hmsc<- list(pred_hm_mean=pred_hm_mean,pred_hm_05=pred_hm_05,pred_hm_95=pred_hm_95)
 return(pred_hmsc)
}

################################################Models for this test#############

jsdm_files_list<-list(EnvEvenSp5="model-2019-04-09-19-02-16.rda",EnvEvenSp10="model-2019-04-10-08-26-20.rda",EnvEvenSp20 ="model-2019-04-11-19-06-02.rda",
                      FacDenseSp5="model-jsdm-2019-05-06-18-41-20.rda",FacDenseSp10="model-2019-04-12-06-59-42.rda",FacDenseSp20="model-2019-04-13-17-22-12.rda",
                      FacSparseSp5 ="model-2019-04-13-17-51-16.rda",FacSparseSp10="model-2019-04-14-05-17-58.rda",FacSparseSp20 ="model-2019-04-15-15-34-20.rda",
                      CompDenseSp5 ="model-2019-04-15-16-03-39.rda",CompDenseSp10="model-2019-04-16-03-39-17.rda",CompDenseSp20="model-2019-04-17-14-00-45.rda",
                      CompSparseSp5="model-2019-04-17-14-30-01.rda",CompSparseSp10 ="model-2019-04-18-01-55-45.rda",CompSparseSp20="model-2019-04-19-12-40-30.rda"
                      , FacCompDenseSp5="model-2019-04-19-13-08-47.rda",FacCompDenseSp10="model-2019-04-20-00-37-25.rda",FacCompDenseSp20="model-2019-04-21-11-10-34.rda"
                      , FacCompSparseSp5="model-2019-04-21-11-39-52.rda",FacCompSparseSp10="model-2019-04-21-22-58-58.rda",FacCompSparseSp20="model-2019-04-23-09-42-18.rda")


gjam_files_list<-list(EnvEvenSp5 ="./gjam_models/gjam5env.rda",EnvEvenSp10="./gjam_models/gjam10env.rda",EnvEvenSp20="./gjam_models/gjam20env.rda",
                      FacDenseSp5="./gjam_models/gjam5f.rda",FacDenseSp10="./gjam_models/gjam10fd.rda",FacDenseSp20="./gjam_models/gjam20fd.rda",
                      FacSparseSp5="./gjam_models/gjam5fs.rda",FacSparseSp10="./gjam_models/gjam10fs.rda",FacSparseSp20="./gjam_models/gjam20fs.rda",
                      CompDenseSp5="./gjam_models/gjam5cmpd.rda",CompDenseSp10="./gjam_models/gjam10cmpd.rda",CompDenseSp20="./gjam_models/gjam20cmpd.rda",
                      CompSparseSp5="./gjam_models/gjam5cmps.rda",CompSparseSp10="./gjam_models/gjam10cmps.rda",CompSparseSp20="./gjam_models/gjam20cmps.rda",
                      FacCompDenseSp5="./gjam_models/gjam5cmp_facd2.rda",FacCompDenseSp10="./gjam_models/gjamcmp_facd10.rda",FacCompDenseSp20="./gjam_models/gjamcmp_facd20.rda",
                      FacCompSparseSp5="./gjam_models/gjamcmp_facs5.rda",FacCompSparseSp10="./gjam_models/gjamcmp_facs10.rda",FacCompSparseSp20="./gjam_models/gjamcmp_facs20.rda")



hmsc_files_list<-list(EnvEvenSp5 ="./new_HMmodels/hm5env.rda",EnvEvenSp10="./new_HMmodels/hm10env.rda" ,EnvEvenSp20="./new_HMmodels/hm20env.rda" 
                      ,FacDenseSp5 ="./new_HMmodels/hm5fd.rda" ,FacDenseSp10 ="./new_HMmodels/hm10fd.rda",FacDenseSp20="./new_HMmodels/hm20fd.rda",
                      FacSparseSp5="./new_HMmodels/hm5fs.rda",FacSparseSp10="./new_HMmodels/hm10fs.rda",FacSparseSp20 ="./new_HMmodels/hm20fs.rda" ,
                      CompDenseSp5="./new_HMmodels/hm5cmpd.rda",CompDenseSp10 ="./new_HMmodels/hm10cmpd.rda" ,CompDenseSp20="./new_HMmodels/hm20cmpd.rda" ,
                      CompSparseSp5="./new_HMmodels/hm5cmps.rda",CompSparseSp10  ="./new_HMmodels/hm10cmps.rda",CompSparseSp20="./new_HMmodels/hm20cmps.rda",
                      FacCompDenseSp5 ="./new_HMmodels/hmcmp_facd5.rda",FacCompDenseSp10="./new_HMmodels/hmcmp_facd10.rda" ,FacCompDenseSp20 ="./new_HMmodels/hmcmp_facd20.rda" ,
                      FacCompSparseSp5 ="./new_HMmodels/hmcmp_facs5.rda",FacCompSparseSp10="./new_HMmodels/hmcmp_facs10.rda" ,FacCompSparseSp20="./new_HMmodels/hmcmp_facs20.rda" )

################################################################################
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

################Convergence check################################################
 jsdm_list<- lapply(jsdm_files_list, load_object)
 
# jsdm_metrics<-mapply(metrics_jsdm,jsdm_list,fac_inter,comp_inter,SIMPLIFY = FALSE)
# jsdm_rho_list<- lapply(jsdm_list,jsdm_rho )
# jsdm_tau_list<- lapply(jsdm_list,jsdm_tau )
# jsdm_convergence_par<- lapply(jsdm_list, jsdm_conv)
# jsdm_conv_dataset<- data.frame()
# 
# for(j in 1:length(jsdm_list)){
#   tmp1<- as.data.frame(jsdm_convergence_par[[j]][1])
#   colnames(tmp1)<-c("parameter", "effective_size", "value")
#   tmp1<-tmp1[,-2]
#   tmp1$type<-"Effective Size"
#   tmp1$Filtering<-sim_names[[j]]
#   tmp2<- as.data.frame(jsdm_convergence_par[[j]][2])
#   colnames(tmp2)<-c("parameter", "Rhat", "value")
#   tmp2<-tmp2[,-2]
#   tmp2$type<-"Rhat"
#   tmp2$Filtering<-sim_names[[j]]
#   jsdm_conv_dataset<-rbind(jsdm_conv_dataset,tmp1,tmp2)  
# }
# 
# for(i in 1:dim(jsdm_conv_dataset)[1]){
#   if(length(grep("Fac",jsdm_conv_dataset$Filtering[i]))>0){ jsdm_conv_dataset$gentype[i]<-"Facilitation"}
#   if(length(grep("Comp",jsdm_conv_dataset$Filtering[i]))>0){ jsdm_conv_dataset$gentype[i]<-"Competition"}
#   if(length(grep("Env",jsdm_conv_dataset$Filtering[i]))>0){ jsdm_conv_dataset$gentype[i]<-"Environmental"}
#   if((length(grep("Fac",jsdm_conv_dataset$Filtering[i]))>0)&(length(grep("Comp",jsdm_conv_dataset$Filtering[i]))>0)){
#     jsdm_conv_dataset$gentype[i]<-"Comp+Fac"
#   }
# }
# #pdf("plot_conv_jsdm.pdf")
# p2<- ggplot(jsdm_conv_dataset, aes(x=value, color=parameter,fill=parameter)) +
#   geom_histogram( alpha=0.4, position="identity") +
#   scale_color_brewer(palette="Dark2")+scale_fill_brewer(palette="Dark2")+ xlab("Rhat") +facet_wrap(type~.,scales = "free")+xlab(" ")+
#   ggtitle("Convergence parameters for the CM model")+theme_bw()+ theme(plot.title = element_text(hjust = 0.5))
# plot(p2)
# #dev.off()

########################################Data################################################################################
data_env<-sim_data$EnvEvenSp5
data_fac<-sim_data$FacDenseSp5
data_cmp<-sim_data$CompDenseSp5

#####True interactions 
#pdf(paste0("plot_test_caseSp7S10.pdf"))
par(mfrow=c(2,2))
cols = colorRampPalette(c("dark blue","white","red"))
col2 <- colorRampPalette(c("#4393C3", "#2166AC", "#053061",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#67001F", "#B2182B", "#D6604D", "#F4A582"))

gcols = colorRampPalette(c( "White", "White", "Black"))
corrplot(diag(5), diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
title("True interactions for environment only")
corrplot(fac_inter[[4]], diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
title("True interactions sparse facilitation")
corrplot((-1)*comp_inter[[10]], diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
title("True interactions sparse competition")


############################Fundamental niche and occurence for 3 datasets#############################################################

run_plots<-function(run_data,jsdm,S1,S2,lab, lab2){
  
 # pdf(paste0("plot_",lab,".pdf"))
  table_fundamental<-data.frame()
  nspecies<-ncol(run_data)-1
  np<-nrow(run_data)
  niche_optima  = seq(2, 98, length.out=nspecies)
  niche_breadth = 20
  xx<-run_data$env
  
  
  #dataframe for fundamental niches
  table_fundamental<-data.frame()
  
  for(i in 1:nspecies){
    tmp<-data.frame(xx=0:100,niche=dnorm(0:100,mean=niche_optima[i],sd=niche_breadth)/dnorm(niche_optima[i],mean=niche_optima[i],sd=niche_breadth),species=i)
    table_fundamental<-rbind(table_fundamental,tmp)
  }
  
  ############data frame for observations
  #table for observations
  Y_data = subset(run_data, select = -env)
  table_obs<-data.frame()
  
  for(i in 1:nspecies){
    tmp<-data.frame(xx=xx,obs=Y_data[,i],species=rep(i,np))
    table_obs<-rbind(table_obs,tmp)
  }
  
  p_0<- list()
  for(i in 1:nspecies)
    local({
      i<-i
      tmp_obs<-table_obs[which(table_obs$species==i),]
      tmp_fund<-table_fundamental[which(table_fundamental$species==i),]
      g<<-ggplot()+
        geom_point(aes(x=tmp_obs$xx,y=tmp_obs$obs),col="#000066",size=0.5) +xlab("Environmental gradient")+ylab("Probability of presence")+
        geom_line(data=tmp_fund,aes(x=tmp_fund$xx,y=tmp_fund$niche,color = "Fundamental niche"),lwd=1)+
        labs(title=paste0("Fundamental niche and realized points, Species ",i, " for ",lab2))+
        scale_color_manual(name = c("Legend"), values = c("Predicted probability" = "#FF6666","Fundamental niche"="#9999FF"))+theme_bw()+ theme(plot.title = element_text(hjust = 0.5))
      p_0<<-list.append(p_0,g)
      assign(paste0("p",i), g, pos =1)
   })
  
  tmp_obsS7<-table_obs[which(table_obs$species==S1),]
  tmp_obsS10<-table_obs[which(table_obs$species==S2),]
  tmp_fundS7<-table_fundamental[which(table_fundamental$species==S1),]
  tmp_fundS10<-table_fundamental[which(table_fundamental$species==S2),]
  color1=paste0("Fundamental niche S",S1)
  color2=paste0("Fundamental niche S",S2)
  g<-ggplot()+geom_point(aes(x=tmp_obsS7$xx,y=tmp_obsS7$obs),col="#FF6666",size=2.8,shape=1)+
    geom_point(aes(x=tmp_obsS10$xx,y=tmp_obsS10$obs),col="#9999FF",size=2.8,shape=2) +xlab("Environmental gradient")+ylab("Probability of presence")+
    geom_line(data=tmp_fundS7,aes(x=tmp_fundS7$xx,y=tmp_fundS7$niche,color =color1),lwd=1)+
    geom_line(data=tmp_fundS10,aes(x=tmp_fundS10$xx,y=tmp_fundS10$niche,color =color2),lwd=1)+
    labs(title=paste0("Fundamental niche and realized points, Species,",S1," and ", S2," for ",lab2) )+
    scale_color_manual(name = c("Legend"), values = c("#FF6666","#9999FF"), labels=c(color1, color2))+theme_bw()+ theme(plot.title = element_text(hjust = 0.5))
  
  plot(g)
  
  #######################################################################################################################
  #######################################################################################################################
  
  table_jsdm<-data.frame()
  pred_jsdm<- predict_cm(run_data,jsdm)
  pred_j_mean<-pred_jsdm$pred_j_mean
  pred_j_05<-pred_jsdm$pred_j_05
  pred_j_95<-pred_jsdm$pred_j_95
  
  
  for(i in 1:nspecies) {
    tmp<-data.frame(xx=xx,mean=pred_j_mean[,i],q_95=pred_j_95[,i],q_05=pred_j_05[,i],type=rep("CM",np),species=rep(i,np))
    table_jsdm<-rbind(table_jsdm,tmp)
    
  }
  
  #jsdm 
  
  p_1<-list()
  for(i in 1:nspecies)
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
        labs(title=paste0("CM, Species ",i, " for ",lab2))+
        scale_color_manual(name = c("Legend"), values = c("Predicted probability" = "#FF6666","Fundamental niche"="#9999FF"))+theme_bw()+ theme(plot.title = element_text(hjust = 0.5))
      
      p_1<<-list.append(p_1,g)
      assign(paste0("p",i), g, pos =1)
    })
  
  
 plot(p_1[[S1]])
 plot(p_1[[S2]])
  
  tmpS710<-table_jsdm[which(table_jsdm$species %in%c(S1,S2)),]
  tmp_obsS710<-table_obs[which(table_obs$species %in%c(S1,S2)),]
  tmp_fundS710<-table_fundamental[which(table_fundamental$species%in%c(S1,S2)),]
  g<-ggplot()+geom_ribbon(aes(x=tmpS710$xx,ymin=tmpS710$q_05,ymax=tmpS710$q_95, col =as.factor(tmpS710$species)),alpha=0.5)+
    geom_line(aes(x=tmpS710$xx,y=tmpS710$mean,col = as.factor(tmpS710$species)),lwd=1.5)+
    geom_point(aes(x=tmp_obsS710$xx,y=tmp_obsS710$obs,col= as.factor(tmp_obsS710$species)),size=0.5) +xlab("Environmental gradient")+ylab("Probability of presence")+
    geom_line(data=tmp_fundS710,aes(x=tmp_fundS710$xx,y=tmp_fundS710$niche,color = as.factor(tmp_fundS710$species)),lwd=1)+
    labs(title=paste0("Fundamental niche and realized points for Species ",S1," and ",S2, " for ",lab2))+theme_bw()+ theme(plot.title = element_text(hjust = 0.5))+
    scale_color_manual(name = c("Legend"), values = c("#FF6666","#9999FF"), labels=c(paste0("Species ",S1),paste0("Species ",S2)))+theme_bw()+ theme(plot.title = element_text(hjust = 0.5))
  plot(g)
#  dev.off()
  
}
#pdf("plot_S1S2_all5.pdf")
#####True interactions 
#pdf(paste0("plot_test_caseSp7S10.pdf"))
par(mfrow=c(2,2))
cols = colorRampPalette(c("dark blue","white","red"))
col2 <- colorRampPalette(c("#4393C3", "#2166AC", "#053061",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#67001F", "#B2182B", "#D6604D", "#F4A582"))

gcols = colorRampPalette(c( "White", "White", "Black"))
corrplot(diag(5), diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
title("True interactions for environment only")
#pdf("plot_compdenseSp20.pdf")
corrplot(fac_inter[[4]], diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
title("True interactions dense competition + facilitation 20 species")
#dev.off()
corrplot((-1)*comp_inter[[10]], diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
title("True interactions sparse competition")

#dev.off()

run_plots(data_env,jsdm=jsdm_list$EnvEvenSp5,S1=1,S2=2,"envS1S2", lab2=" environment")
run_plots(data_fac,jsdm=jsdm_list$FacDenseSp5,S1=2,S2=3,"facS2S3", lab2=" facilitation")
run_plots(data_cmp,jsdm=jsdm_list$CompDenseSp5,S1=1,S2=2,"cmpS1S12", lab2=" competition")


#dev.off()


###########################Add SDM model################################################################


predict_sdm_glm<-function(datas){
  data_temp<-dataprep(datas) 
  df<-cbind(data_temp$X[,-1],data_temp$Y)
  names(df)[1:2]<- c("env","env2" )
  #predY<- as.data.frame(matrix(NA,nrow(data_temp$X),ncol(data_temp$Y)))
  predY<- data.frame()
  predBeta<- data.frame()
  #colnames(predBeta)<- c("low","high","mean","sp","beta")
  i<-1
  for(i in 1:ncol(data_temp$Y)){
    df_temp<-df[,c("env","env2",colnames(df)[i+2])]
    colnames(df_temp)<-c("env","env2","yval")
    probitMod <- glm(yval~ env + env2, data=df_temp, family=binomial(link="probit"))  # build the logit model
    #predicted <- predict(probitMod, df[,c("env","env2")], type="response")  # predict the probability scores
    predicted <- predict(probitMod, df[,c("env","env2")], type="link",se.fit = TRUE)  # predict the probability scores
    table_sdm<-as.data.frame(datas$env)
    colnames(table_sdm)<-c("xx")
    table_sdm$mean<-pnorm(predicted$fit)
    confidence <- .95
    score <- qnorm((confidence / 2) + .5)
    table_sdm$q_95<-pnorm(predicted$fit + score * predicted$se.fit)
    table_sdm$q_05<-pnorm(predicted$fit - score * predicted$se.fit)
    table_sdm$type<-"SDM"
    table_sdm$species<-i
    predY<- rbind(predY,table_sdm)
    tmp_ci<-confint(probitMod)
    predBeta_tmp<- as.data.frame(tmp_ci)
    names(predBeta_tmp)<- c("low", "high")
    predBeta_tmp$mean <- probitMod$coefficients
    predBeta_tmp$sp <- rep(i,3)
    predBeta_tmp$beta<- c("inter","env","env2")
    predBeta<- rbind(predBeta,predBeta_tmp)
  }
  
  predBeta$model<-"SDM"
  
  return(list(beta_sdm=predBeta,pred_sdm=predY))
}


predict_sdm_bayes<-function(datas){
  data_temp<-dataprep(datas) 
  predY<- data.frame()
  predBeta<- data.frame()
  for(s in 1:ncol(data_temp$Y)){
    y_temp<-data_temp$Y[,s]
   # sim.par<-uni_bayes_probit(datas,D=3, i_col=s, N_sim=20000)
    sim.par <- bayes.probit(y_temp, data_temp$X, 20000)
    betas<- apply(sim.par$beta[-(1:5000),], 2, quantile, c(0.05, 0.5, 0.95))
    pred_j<-array(NA,dim=c(data_temp$n,15000))
    for(j in 1:15000){
      pred_j[,j]<-pnorm(data_temp$X%*%sim.par$beta[j+5000,])
      
    }
    pred_j_mean <- apply(pred_j, 1, mean)
    pred_j_05 <- apply(pred_j, 1, quantile,0.05)
    pred_j_95 <- apply(pred_j, 1, quantile,0.95)
    table_sdm<-as.data.frame(datas$env)
    colnames(table_sdm)<-c("xx")
    table_sdm$mean<-pred_j_mean
    table_sdm$q_95<-pred_j_05
    table_sdm$q_05<-pred_j_95
    table_sdm$type<-"SDM"
    table_sdm$species<-s
    predY<- rbind(predY,table_sdm)
    predBeta_tmp<- as.data.frame( t(betas))
    names(predBeta_tmp)<- c("low","mean", "high")
    predBeta_tmp$sp <- rep(s,3)
    predBeta_tmp$beta<- c("inter","env","env2")
    predBeta<- rbind(predBeta,predBeta_tmp)
  }
  
  predBeta$model<-"SDM"
  
  return(list(beta_sdm=predBeta,pred_sdm=predY))
}

#library(LearnBayes)

#datab<- dataprep(sim_data$EnvEvenSp5)
#sim.par <- bayes.probit(datab$Y[,1], datab$X, 20000)


sdm_b_env5 <- predict_sdm_bayes(sim_data$EnvEvenSp5)
sdm_env5 <- predict_sdm_glm(sim_data$EnvEvenSp5)



###########################Add other models#############################################################
plot_response<- function(datab,pred_gjam, pred_jsdm,pred_hmsc,sdmpred){
  np<-nrow(datab)
  nspecies<-ncol(datab)-1
  nsp<- nspecies
  niche_optima  = seq(2, 98, length.out=nspecies)
  niche_breadth = 20
  
  #dataframe for fundamental niches
  table_fundamental<-data.frame()
  
  for(i in 1:nspecies){
    tmp<-data.frame(xx=0:100,niche=dnorm(0:100,mean=niche_optima[i],sd=niche_breadth)/dnorm(niche_optima[i],mean=niche_optima[i],sd=niche_breadth),species=i)
    table_fundamental<-rbind(table_fundamental,tmp)
  }
  
  ######## JSDM
  xx<-datab$env
  
  table_jsdm<-data.frame()
  pred_j_mean<-pred_jsdm$pred_j_mean
  pred_j_05<-pred_jsdm$pred_j_05
  pred_j_95<-pred_jsdm$pred_j_95
  
  
  for(i in 1:nspecies) {
    tmp<-data.frame(xx=xx,mean=pred_j_mean[,i],q_95=pred_j_95[,i],q_05=pred_j_05[,i],type=rep("CM",np),species=rep(i,np))
    table_jsdm<-rbind(table_jsdm,tmp)
    
  }
  
  
  
  ######## GJAM
  table_gjam<-data.frame()
  
  pred_gj_mean<-pred_gjam$pred_gj_mean
  pred_gj_95<-pred_gjam$pred_gj_95
  pred_gj_05<-pred_gjam$pred_gj_05
  
  for(i in 1:nspecies) {
    tmp<-data.frame(xx=xx,mean=pred_gj_mean[,i],q_95=pred_gj_95[,i],q_05=pred_gj_05[,i],type=rep("GJAM",np),species=rep(i,np))
    table_gjam<-rbind(table_gjam,tmp)
  }
  
  
  # ######## GJAM_DR
  # table_gjam_dr<-data.frame()
  # 
  # pred_gj_dr_mean<-pred_gjam_dr$pred_gj_dr_mean
  # pred_gj_dr_95<-pred_gjam_dr$pred_gj_dr_95
  # pred_gj_dr_05<-pred_gjam_dr$pred_gj_dr_05
  # 
  # for(i in 1:nspecies) {
  #   tmp<-data.frame(xx=xx,mean=pred_gj_dr_mean[,i],q_95=pred_gj_dr_95[,i],q_05=pred_gj_dr_05[,i],type=rep("DR-GJAM",np),species=rep(i,np))
  #   table_gjam_dr<-rbind(table_gjam_dr,tmp)
  # }
  # 
  
  #HMSC
  table_hmsc<-data.frame()
  pred_hm_mean<-pred_hmsc$pred_hm_mean
  pred_hm_95<-pred_hmsc$pred_hm_95
  pred_hm_05<-pred_hmsc$pred_hm_05
  
  for(i in 1:nspecies) {
    tmp<-data.frame(xx=xx,mean=pred_hm_mean[,i],q_95=pred_hm_95[,i],q_05=pred_hm_05[,i],type=rep("HMSC",np),species=rep(i,np))
    table_hmsc<-rbind(table_hmsc,tmp)
  }
  
  #SDM 
  table_sdm<-sdmpred
  
  #table for predictions
  table<-rbind(table_jsdm,table_gjam,table_hmsc,table_sdm)
  
  #table for observations
  Y_data = subset(datab, select = -env)
  table_obs<-data.frame()
  
  for(i in 1:nspecies){
    tmp<-data.frame(xx=xx,obs=Y_data[,i],species=rep(i,np))
    table_obs<-rbind(table_obs,tmp)
  }
  tjur_list<-list()
  ######Tjur
  for(k in (1:nspecies)) {
    datap<- dataprep(datab)
    indx <- datap$Y[,k]==1
    Tjur_j <- mean(pred_jsdm$pred_j_mean[indx,k]) - mean(pred_jsdm$pred_j_mean[!indx,k])
    Tjur_gj <- mean(pred_gjam$pred_gj_mean[indx,k]) - mean(pred_gjam$pred_gj_mean[!indx,k])
    Tjur_hm <- mean(pred_hmsc$pred_hm_mean[indx,k]) - mean(pred_hmsc$pred_hm_mean[!indx,k])
    Tjur_sdm<- mean(table_sdm$mean[which(table_sdm$species==k)][indx]) - mean(table_sdm$mean[which(table_sdm$species==k)][!indx])
    Tjur<- sum(Tjur_j,Tjur_gj,Tjur_hm,Tjur_sdm)/4
    tjur_list<-c(tjur_list,Tjur)
  }
  
  #######Auc
  Auc_list<- list()
  for(k in (1:nspecies)) {
    datap<- dataprep(datab)
    indx <- datap$Y[,k]==1
    AUC_j <- auc(roc(pred_jsdm$pred_j_mean[,k],factor(datap$Y[,k]))) 
    AUC_gj <- auc(roc(pred_gjam$pred_gj_mean[,k],factor(datap$Y[,k])))
    AUC_hm <-  auc(roc(pred_hmsc$pred_hm_mean[,k],factor(datap$Y[,k])))  
    AUC_sdm<-  auc(roc(table_sdm$mean[which(table_sdm$species==k)],factor(datap$Y[,k])))
    AUC_mean<- sum(AUC_j,AUC_gj,AUC_hm,AUC_sdm)/4
    Auc_list<-c(Auc_list,AUC_mean)
  }

  ########SECOND TYPE OF PLOT
  p_4<-list()
  
  
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
        scale_color_manual(name = c("Legend"), values = c("SDM" = "#56B4E9","CM" = "#FF6666","GJAM" = "#66CC99","DR-GJAM" = "#FFFF10","HMSC" = "#FFB266","Fundamental niche"="#9999FF")) +theme_bw()+ theme(plot.title = element_text(hjust = 0.5))
      p_4<<-list.append(p_4,g)
      assign(paste0("p",i), g, pos =1)
    })
  
  if (nsp==5){
    grid.arrange(p1,p2,p3,p4,p5,nrow=ceiling(nspecies/2))
  }else{ 
    p_4
  }
  
 return(p_4)
}


#####################################HMSC##############
###############################################HMSC########################################################

#hm_conv(hm_mod)
#hm_mod_R<-hm_inter(hm_mod)
# 
# h_metric_e20<-metrics_hmsc(hm_mod_R$Rho_sign)
# h_metric_e20_p<-metrics_hmsc(hm_mod_R$Tau_sign)
# 
# mod_list_Rho<-list.append(mod_list_Rho,hmsc=hm_mod_R$Rho_sign)
# mod_list_Tau<-list.append(mod_list_Tau,hmsc=hm_mod_R$Tau_sign)


###############################################GJAM########################################################

data_test<- sim_data[c(1,4,11)]

gjam_list<- lapply(gjam_files_list[c(1,4,11)],gj_conv)
gjam_conv_dataset<- data.frame()
# 
# for(j in 1:length(gjam_list)){
#   tmp1<- as.data.frame(gjam_list[[j]][1])
#   colnames(tmp1)<-c("value","parameter")
#   tmp1$type<-"Effective Size"
#   tmp1$Filtering<-sim_names[[j]]
#   tmp2<- as.data.frame(gjam_list[[j]][2])
#   colnames(tmp2)<-c("value", "value2","parameter")
#   tmp2<-tmp2[,-2]
#   tmp2$type<-"Rhat"
#   tmp2$Filtering<-sim_names[[j]]
#   gjam_conv_dataset<-rbind(gjam_conv_dataset,tmp1,tmp2)  
# }
# 
# for(i in 1:dim(gjam_conv_dataset)[1]){
#   if(length(grep("Fac",gjam_conv_dataset$Filtering[i]))>0){ gjam_conv_dataset$gentype[i]<-"Facilitation"}
#   if(length(grep("Comp",gjam_conv_dataset$Filtering[i]))>0){ gjam_conv_dataset$gentype[i]<-"Competition"}
#   if(length(grep("Env",gjam_conv_dataset$Filtering[i]))>0){ gjam_conv_dataset$gentype[i]<-"Environmental"}
#   if((length(grep("Fac",gjam_conv_dataset$Filtering[i]))>0)&(length(grep("Comp",gjam_conv_dataset$Filtering[i]))>0)){
#     gjam_conv_dataset$gentype[i]<-"Comp+Fac"
#   }
# }
# 
# #pdf("plot_conv_gjam.pdf")
# p2<- ggplot(gjam_conv_dataset, aes(x=value, color=parameter,fill=parameter)) +
#   geom_histogram( alpha=0.4, position="identity") +
#   scale_color_brewer(palette="Dark2")+scale_fill_brewer(palette="Dark2") +facet_wrap(type~.,scales = "free_x") +xlab(" ")+
#   ggtitle("Convergence parameters for the GJAM model")+theme_bw()+ theme(plot.title = element_text(hjust = 0.5))
# p2
# #dev.off()

###########prediction############################################################################################
#ind<-c(1,4,10,16)
ind<-c(5)

pdf("Plots_CAM_10dense.pdf")
data_test<- sim_data[ind]
hmsc_files_l<- hmsc_files_list[ind]
gjam_files_l<-gjam_files_list[ind]
jsdm_mod_l<-jsdm_list[ind]

p_list<-list()
for(l in 1:length(data_test)){
  pred_gjam<-load_gjam(data_test[[l]],gjam_files_l[[l]])
  pred_hmsc<- hm_pred(data_test[[l]],hmsc_files_l[[l]])
  #sdm_pr<-predict_sdm_glm(data_test[[l]])
  sdm_b_pr<-predict_sdm_bayes(data_test[[l]])
  pred_jsdm<- predict_cm(data_test[[l]],jsdm_mod_l[[l]])
  p_plot<- plot_response(datab=data_test[[l]],pred_gjam$predict,pred_jsdm,pred_hmsc,sdm_b_pr$pred_sdm)
  p_plot
  p_list<- list.append(p_list, p_plot)
}
p_list
dev.off()


pred_jsdm$pred_j_mean[,2]
tjur<- c()
for(k in (1:ncol(pred_jsdm$pred_j_mean))) {
  datap<- dataprep(data_test[[l]])
  indx <- datap$Y[,k]==1
  Tjur <- mean(pred_jsdm$pred_j_mean[indx,k]) - mean(pred_jsdm$pred_j_mean[!indx,k])
  tjur<-c(tjur,Tjur)
}

tjur2<- c()
for(k in (1:ncol(pred_jsdm$pred_j_mean))) {
  datap<- dataprep(data_test[[l]])
  indx <- datap$Y[,k]==1
  Tjur <- mean(pred_hmsc$pred_hm_mean[indx,k]) - mean(pred_hmsc$pred_hm_mean[!indx,k])
  tjur2<-c(tjur2,Tjur)
}
#dev.off()
################################################################################################################

