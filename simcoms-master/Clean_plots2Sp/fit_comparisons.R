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

############Functions HMSC###################################################################################################


################################################Models for this test#############
jsdm_files_list<-list(EnvEvenSp10="model-2019-04-10-08-26-20.rda",
                       FacDenseSp10="model-2019-04-12-06-59-42.rda",
                      CompDenseSp10="model-2019-04-16-03-39-17.rda")

gjam_files_list<-list(EnvEvenSp10="./gjam_models/gjam10env.rda",
                      FacDenseSp10="./gjam_models/gjam10fd.rda",
                      CompDenseSp10="./gjam_models/gjam10cmpd.rda")


hmsc_files_list<-list(EnvEvenSp10="./new_HMmodels/hm10env.rda",
                      FacDenseSp10="./new_HMmodels/hm10fd.rda",
                      CompDenseSp10="./new_HMmodels/hm10cmpd.rda")

################################################################################


################Convergence check################################################
jsdm_list<- lapply(jsdm_files_list, load_object)

jsdm_metrics<-mapply(metrics_jsdm,jsdm_list,fac_inter[c(2,5,11)],comp_inter[c(2,5,11)],SIMPLIFY = FALSE)
jsdm_rho_list<- lapply(jsdm_list,jsdm_rho )
jsdm_tau_list<- lapply(jsdm_list,jsdm_tau )
jsdm_convergence_par<- lapply(jsdm_list, jsdm_conv)
jsdm_conv_dataset<- data.frame()

for(j in 1:length(jsdm_list)){
  tmp1<- as.data.frame(jsdm_convergence_par[[j]][1])
  colnames(tmp1)<-c("parameter", "effective_size", "value")
  tmp1<-tmp1[,-2]
  tmp1$type<-"Effective Size"
  tmp1$Filtering<-sim_names[[j]]
  tmp2<- as.data.frame(jsdm_convergence_par[[j]][2])
  colnames(tmp2)<-c("parameter", "Rhat", "value")
  tmp2<-tmp2[,-2]
  tmp2$type<-"Rhat"
  tmp2$Filtering<-sim_names[[j]]
  jsdm_conv_dataset<-rbind(jsdm_conv_dataset,tmp1,tmp2)  
}

for(i in 1:dim(jsdm_conv_dataset)[1]){
  if(length(grep("Fac",jsdm_conv_dataset$Filtering[i]))>0){ jsdm_conv_dataset$gentype[i]<-"Facilitation"}
  if(length(grep("Comp",jsdm_conv_dataset$Filtering[i]))>0){ jsdm_conv_dataset$gentype[i]<-"Competition"}
  if(length(grep("Env",jsdm_conv_dataset$Filtering[i]))>0){ jsdm_conv_dataset$gentype[i]<-"Environmental"}
  if((length(grep("Fac",jsdm_conv_dataset$Filtering[i]))>0)&(length(grep("Comp",jsdm_conv_dataset$Filtering[i]))>0)){
    jsdm_conv_dataset$gentype[i]<-"Comp+Fac"
  }
}
#pdf("plot_conv_jsdm.pdf")
p2<- ggplot(jsdm_conv_dataset, aes(x=value, color=parameter,fill=parameter)) +
  geom_histogram( alpha=0.4, position="identity") +
  scale_color_brewer(palette="Dark2")+scale_fill_brewer(palette="Dark2")+ xlab("Rhat") +facet_wrap(type~.,scales = "free")+xlab(" ")+
  ggtitle("Convergence parameters for the CM model")+theme_bw()+ theme(plot.title = element_text(hjust = 0.5))
plot(p2)
#dev.off()

########################################Data################################################################################
data_env<-sim_data$EnvEvenSp10
data_fac<-sim_data$FacDenseSp10
data_cmp<-sim_data$CompDenseSp10

#####True interactions 
#pdf(paste0("plot_test_caseSp7S10.pdf"))
par(mfrow=c(2,2))
cols = colorRampPalette(c("dark blue","white","red"))
col2 <- colorRampPalette(c("#4393C3", "#2166AC", "#053061",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#67001F", "#B2182B", "#D6604D", "#F4A582"))

gcols = colorRampPalette(c( "White", "White", "Black"))
corrplot(diag(10), diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
title("True interactions for environment only")
corrplot(fac_inter[[5]], diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
title("True interactions facilitation")
corrplot((-1)*comp_inter[[11]], diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
title("True interactions competition")


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
#pdf("plot_S4S6.pdf")
run_plots(data_env,jsdm=jsdm_list$EnvEvenSp10,S1=4,S2=6,"envS7S10", lab2=" environment")
run_plots(data_fac,jsdm=jsdm_list$FacDenseSp10,S1=4,S2=6,"facS7S10", lab2=" facilitation")
#run_plots(data_cmp,jsdm=jsdm_list$CompDenseSp10,S1=7,S2=10,"cmpS7S10", lab2=" competition")


#dev.off()


###########################Add SDM model################################################################



###########################Add other models#############################################################


