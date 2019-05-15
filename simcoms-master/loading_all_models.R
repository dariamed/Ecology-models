####Load all models
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



############Functions JSDM###################################################################################################
jsdm_conv<-function(mod) {
  n_eff_beta<- as.data.frame(t(mod$n.eff$B))
  n_eff_beta$parameter<- "beta"
  n_eff_sigma<- as.data.frame(mod$n.eff$Rho)
  n_eff_sigma$parameter<- "correlation"
  neff<- rbind(n_eff_beta,n_eff_sigma)
  neff_mod<- melt(neff, id=c("parameter"), variable_name = "effective size")
  # p<- ggplot(neff_mod, aes(x=value, color=parameter,fill=parameter)) +
  #   geom_histogram( alpha=0.4, position="identity",binwidth=100) + xlab("effective size") +
  #   scale_color_brewer(palette="Dark2")+scale_fill_brewer(palette="Dark2")+
  #   ggtitle("Effective size for the parameters beta and coefficients of correlation matrix")
  # plot(p)
  Rhat_beta<- as.data.frame(t(mod$Rhat$B))
  Rhat_beta$parameter<- "beta"
  Rhat_sigma<- as.data.frame(mod$Rhat$Rho)
  Rhat_sigma$parameter<- "correlation"
  Rhat<- rbind(Rhat_beta,Rhat_sigma)
  Rhat_mod<- melt(Rhat, id=c("parameter"), variable_name = "Rhat")
  Rhat_mod_nonNA<- Rhat_mod[!is.na(Rhat_mod$value),]
  # p2<- ggplot(Rhat_mod_nonNA, aes(x=value, color=parameter,fill=parameter)) +
  #   geom_histogram( alpha=0.4, position="identity",binwidth = 0.01) +
  #   scale_color_brewer(palette="Dark2")+scale_fill_brewer(palette="Dark2")+ xlab("Rhat") +
  #   ggtitle("Rhat for the parameters beta and the coefficients of correlation matrix")
  #plot(p2)
  return(list(neff=neff_mod,Rhat=Rhat_mod_nonNA))
  
}

############Functions JSDM###################################################################################################



##SET UP ##########################################################################################################

jsdm_files_list<-list(EnvEvenSp5="model-2019-04-09-19-02-16.rda",EnvEvenSp10="model-2019-04-10-08-26-20.rda",EnvEvenSp20 ="model-2019-04-11-19-06-02.rda",
                      FacDenseSp5="model-jsdm-2019-05-06-18-41-20.rda",FacDenseSp10="model-2019-04-12-06-59-42.rda",FacDenseSp20="model-2019-04-13-17-22-12.rda",
                      FacSparseSp5 ="model-2019-04-13-17-51-16.rda",FacSparseSp10="model-2019-04-14-05-17-58.rda",FacSparseSp20 ="model-2019-04-15-15-34-20.rda",
                      CompDenseSp5 ="model-2019-04-15-16-03-39.rda",CompDenseSp10="model-2019-04-16-03-39-17.rda",CompDenseSp20="model-2019-04-17-14-00-45.rda",
                      CompSparseSp5="model-2019-04-17-14-30-01.rda",CompSparseSp10 ="model-2019-04-18-01-55-45.rda",CompSparseSp20="model-2019-04-19-12-40-30.rda"
                      , FacCompDenseSp5="model-2019-04-19-13-08-47.rda",FacCompDenseSp10="model-2019-04-20-00-37-25.rda",FacCompDenseSp20="model-2019-04-21-11-10-34.rda"
                      , FacCompSparseSp5="model-2019-04-21-11-39-52.rda",FacCompSparseSp10="model-2019-04-21-22-58-58.rda",FacCompSparseSp20="model-2019-04-23-09-42-18.rda")

jsdm_list<- lapply(jsdm_files_list, load_object)
k<-load_object("model-2019-04-11-19-06-02.rda")

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

pdf("plot_conv_jsdm.pdf")
p2<- ggplot(jsdm_conv_dataset, aes(x=value, color=parameter,fill=parameter)) +
  geom_histogram( alpha=0.4, position="identity") +
  scale_color_brewer(palette="Dark2")+scale_fill_brewer(palette="Dark2")+ xlab("Rhat") +facet_wrap(type~.,scales = "free")+xlab(" ")+
   ggtitle("Convergence parameters for the CM model")+theme_bw()+ theme(plot.title = element_text(hjust = 0.5))
plot(p2)
dev.off()
############################################HMSC functions#####################################################################################
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
  # p<- ggplot(neff, aes(x=value, color=parameter,fill=parameter)) +
  #   geom_histogram( alpha=0.4, position="identity",binwidth =50) + xlab("effective size") +
  #   scale_color_brewer(palette="Dark2")+scale_fill_brewer(palette="Dark2")+
  #   ggtitle("Effective size for the parameters beta and coefficients of correlation matrix")
  # plot(p)
  # 
  Rhat_beta<-as.data.frame(gelman.diag(codaList$Beta,multivariate=FALSE)$psrf)
  Rhat_beta$parameter<- "beta"
  Rhat_sigma<- as.data.frame(gelman.diag(codaList$Omega[[1]], multivariate=FALSE)$psrf)
  Rhat_sigma$parameter<- "correlation"
  Rhat<- rbind(Rhat_beta,Rhat_sigma)
  # p2<- ggplot(Rhat, aes(x= Rhat$`Point est.`, color=parameter,fill=parameter)) +
  #   geom_histogram( alpha=0.4, position="identity",binwidth =0.01) +
  #   scale_color_brewer(palette="Dark2")+scale_fill_brewer(palette="Dark2")+ xlab("Rhat") +
  #   ggtitle("Rhat for the parameters beta and the coefficients of correlation matrix")
  # plot(p2)
  return(list(neff=neff,Rhat=Rhat))
}

###########################################HMSC functions#######################################################################

hmsc_files_list<-list(EnvEvenSp5 ="./new_HMmodels/hm5env.rda",EnvEvenSp10="./new_HMmodels/hm10env.rda" ,EnvEvenSp20="./new_HMmodels/hm20env.rda" 
                      ,FacDenseSp5 ="./new_HMmodels/hm5fd.rda" ,FacDenseSp10 ="./new_HMmodels/hm10fd.rda",FacDenseSp20="./new_HMmodels/hm20fd.rda",
                      FacSparseSp5="./new_HMmodels/hm5fs.rda",FacSparseSp10="./new_HMmodels/hm10fs.rda",FacSparseSp20 ="./new_HMmodels/hm20fs.rda" ,
                      CompDenseSp5="./new_HMmodels/hm5cmpd.rda",CompDenseSp10 ="./new_HMmodels/hm10cmpd.rda" ,CompDenseSp20="./new_HMmodels/hm20cmpd.rda" ,
                      CompSparseSp5="./new_HMmodels/hm5cmps.rda",CompSparseSp10  ="./new_HMmodels/hm10cmps.rda",CompSparseSp20="./new_HMmodels/hm20cmps.rda",
                      FacCompDenseSp5 ="./new_HMmodels/hmcmp_facd5.rda",FacCompDenseSp10="./new_HMmodels/hmcmp_facd10.rda" ,FacCompDenseSp20 ="./new_HMmodels/hmcmp_facd20.rda" ,
                      FacCompSparseSp5 ="./new_HMmodels/hmcmp_facs5.rda",FacCompSparseSp10="./new_HMmodels/hmcmp_facs10.rda" ,FacCompSparseSp20="./new_HMmodels/hmcmp_facs20.rda" )


hmsc_list<- lapply(hmsc_files_list, load_object)
hmsc_convergence_par<- lapply(hmsc_list, hm_conv)

hmsc_conv_dataset<- data.frame()

for(j in 1:length(hmsc_list)){
  tmp1<- as.data.frame(hmsc_convergence_par[[j]][1])
  colnames(tmp1)<-c("value","parameter")
  #tmp1<-tmp1[,-2]
  tmp1$type<-"Effective Size"
  tmp1$Filtering<-sim_names[[j]]
  tmp2<- as.data.frame(hmsc_convergence_par[[j]][2])
  colnames(tmp2)<-c("value", "value2","parameter")
  tmp2<-tmp2[,-2]
  tmp2$type<-"Rhat"
  tmp2$Filtering<-sim_names[[j]]
  hmsc_conv_dataset<-rbind(hmsc_conv_dataset,tmp1,tmp2)  
}

for(i in 1:dim(hmsc_conv_dataset)[1]){
  if(length(grep("Fac",hmsc_conv_dataset$Filtering[i]))>0){ hmsc_conv_dataset$gentype[i]<-"Facilitation"}
  if(length(grep("Comp",hmsc_conv_dataset$Filtering[i]))>0){ hmsc_conv_dataset$gentype[i]<-"Competition"}
  if(length(grep("Env",hmsc_conv_dataset$Filtering[i]))>0){ hmsc_conv_dataset$gentype[i]<-"Environmental"}
  if((length(grep("Fac",hmsc_conv_dataset$Filtering[i]))>0)&(length(grep("Comp",hmsc_conv_dataset$Filtering[i]))>0)){
    hmsc_conv_dataset$gentype[i]<-"Comp+Fac"
  }
}
pdf("plot_conv_hmsc.pdf")
p2<- ggplot(hmsc_conv_dataset, aes(x=value, color=parameter,fill=parameter)) +
  geom_histogram( alpha=0.4, position="identity") +
  scale_color_brewer(palette="Dark2")+scale_fill_brewer(palette="Dark2") +facet_wrap(type~.,scales = "free_x") +xlab(" ")+
  ggtitle("Convergence parameters for the HMSC model")+theme_bw()+ theme(plot.title = element_text(hjust = 0.5))

p2
dev.off()

###########################################GJAM functions#######################################################################


gj_conv<-function(name){
  gj_mod<-load_object(name)
  burn<-gj$m1$modelList$burnin
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



###########################################GJAM functions #######################################################################




gjam_files_list<-list(EnvEvenSp5 ="./gjam_models/gjam5env.rda",EnvEvenSp10="./gjam_models/gjam10env.rda",EnvEvenSp20="./gjam_models/gjam20env.rda",
                      FacDenseSp5="./gjam_models/gjam5f.rda",FacDenseSp10="./gjam_models/gjam10fd.rda",FacDenseSp20="./gjam_models/gjam20fd.rda",
                      FacSparseSp5="./gjam_models/gjam5fs.rda",FacSparseSp10="./gjam_models/gjam10fs.rda",FacSparseSp20="./gjam_models/gjam20fs.rda",
                      CompDenseSp5="./gjam_models/gjam5cmpd.rda",CompDenseSp10="./gjam_models/gjam10cmpd.rda",CompDenseSp20="./gjam_models/gjam20cmpd.rda",
                      CompSparseSp5="./gjam_models/gjam5cmps.rda",CompSparseSp10="./gjam_models/gjam10cmps.rda",CompSparseSp20="./gjam_models/gjam20cmps.rda",
                      FacCompDenseSp5="./gjam_models/gjam5cmp_facd2.rda",FacCompDenseSp10="./gjam_models/gjamcmp_facd10.rda",FacCompDenseSp20="./gjam_models/gjamcmp_facd20.rda",
                      FacCompSparseSp5="./gjam_models/gjamcmp_facs5.rda",FacCompSparseSp10="./gjam_models/gjamcmp_facs10.rda",FacCompSparseSp20="./gjam_models/gjamcmp_facs20.rda")



gjam_list<- lapply(gjam_files_list,gj_conv)
gjam_conv_dataset<- data.frame()

for(j in 1:length(gjam_list)){
  tmp1<- as.data.frame(gjam_list[[j]][1])
  colnames(tmp1)<-c("value","parameter")
  tmp1$type<-"Effective Size"
  tmp1$Filtering<-sim_names[[j]]
  tmp2<- as.data.frame(gjam_list[[j]][2])
  colnames(tmp2)<-c("value", "value2","parameter")
  tmp2<-tmp2[,-2]
  tmp2$type<-"Rhat"
  tmp2$Filtering<-sim_names[[j]]
  gjam_conv_dataset<-rbind(gjam_conv_dataset,tmp1,tmp2)  
}

for(i in 1:dim(gjam_conv_dataset)[1]){
  if(length(grep("Fac",gjam_conv_dataset$Filtering[i]))>0){ gjam_conv_dataset$gentype[i]<-"Facilitation"}
  if(length(grep("Comp",gjam_conv_dataset$Filtering[i]))>0){ gjam_conv_dataset$gentype[i]<-"Competition"}
  if(length(grep("Env",gjam_conv_dataset$Filtering[i]))>0){ gjam_conv_dataset$gentype[i]<-"Environmental"}
  if((length(grep("Fac",gjam_conv_dataset$Filtering[i]))>0)&(length(grep("Comp",gjam_conv_dataset$Filtering[i]))>0)){
    gjam_conv_dataset$gentype[i]<-"Comp+Fac"
  }
}


pdf("plot_conv_gjam.pdf")
p2<- ggplot(gjam_conv_dataset, aes(x=value, color=parameter,fill=parameter)) +
  geom_histogram( alpha=0.4, position="identity") +
  scale_color_brewer(palette="Dark2")+scale_fill_brewer(palette="Dark2") +facet_wrap(type~.,scales = "free_x") +xlab(" ")+
  ggtitle("Convergence parameters for the GJAM model")+theme_bw()+ theme(plot.title = element_text(hjust = 0.5))
p2
dev.off()


#########GJAM DIMENSION REDUCTION########################################################################




#########GJAM DIMENSION REDUCTION########################################################################
gjam_dr_files_list<-list(env5="./gjam_models/gjamDR5env.rda",env10="./gjam_models/gjamDR10env.rda",env20="./gjam_models/gjamDR20env.rda",
                         facd5="./gjam_models/gjamDR5facd.rda",facd10="./gjam_models/gjamDR10facd.rda",facd20="./gjam_models/gjamDR20facd.rda",
                         facs5="./gjam_models/gjamDR5facs.rda",facs10="./gjam_models/gjamDR10facs.rda",facs20="./gjam_models/gjamDR20facs.rda",
                      compd5="./gjam_models/gjamDR5compd.rda",compd10="./gjam_models/gjamDR10compd.rda",compd20="./gjam_models/gjamDR20compd.rda",
                      comps5="./gjam_models/gjamDR5comps.rda",comps10="./gjam_models/gjamDR10comps.rda",comps20="./gjam_models/gjamDR20comps.rda",
                      compfacd5="./gjam_models/gjamDR5compfacd.rda",compfacd10="./gjam_models/gjamDR10compfacd.rda",compfacd20="./gjam_models/gjamDR20compfacd.rda",
                      compfacs5="./gjam_models/gjamDR5compfacs.rda",compfacs10="./gjam_models/gjamDR10compfacs.rda",compfacs20="./gjam_models/gjamDR20compfacs.rda")



gjam_dr_list<- lapply(gjam_dr_files_list,gj_conv)
gjam_dr_conv_dataset<- data.frame()

for(j in 1:length(gjam_dr_list)){
  tmp1<- as.data.frame(gjam_dr_list[[j]][1])
  colnames(tmp1)<-c("value","parameter")
  tmp1$type<-"Effective Size"
  tmp1$Filtering<-sim_names[[j]]
  tmp2<- as.data.frame(gjam_dr_list[[j]][2])
  colnames(tmp2)<-c("value", "value2","parameter")
  tmp2<-tmp2[,-2]
  tmp2$type<-"Rhat"
  tmp2$Filtering<-sim_names[[j]]
  gjam_dr_conv_dataset<-rbind(gjam_dr_conv_dataset,tmp1,tmp2)  
}

for(i in 1:dim(gjam_dr_conv_dataset)[1]){
  if(length(grep("Fac",gjam_dr_conv_dataset$Filtering[i]))>0){ gjam_dr_conv_dataset$gentype[i]<-"Facilitation"}
  if(length(grep("Comp",gjam_dr_conv_dataset$Filtering[i]))>0){ gjam_dr_conv_dataset$gentype[i]<-"Competition"}
  if(length(grep("Env",gjam_dr_conv_dataset$Filtering[i]))>0){ gjam_dr_conv_dataset$gentype[i]<-"Environmental"}
  if((length(grep("Fac",gjam_dr_conv_dataset$Filtering[i]))>0)&(length(grep("Comp",gjam_dr_conv_dataset$Filtering[i]))>0)){
    gjam_dr_conv_dataset$gentype[i]<-"Comp+Fac"
  }
}

pdf("plot_conv_gjamdr.pdf")
p2<- ggplot(gjam_dr_conv_dataset, aes(x=value, color=parameter,fill=parameter)) +
  geom_histogram( alpha=0.4, position="identity") +
  scale_color_brewer(palette="Dark2")+scale_fill_brewer(palette="Dark2") +facet_wrap(type~.,scales = "free_x") +xlab(" ")+
  ggtitle("Convergence parameters for the DR_GJAM model")+theme_bw()+ theme(plot.title = element_text(hjust = 0.5))
p2
dev.off()



##############Mean correlations JSDM####################################################################################

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


mean_cor<-function(mod){
  mean_correlations <- do.call(
    rbind,
    lapply(
      seq_along(mod),
      function(i) {
        x <- mod[[i]]
        nm <- strsplit(
          names(mod)[[i]], "(?<=[a-z])(?=[A-Z])", perl = TRUE
        )[[1]]
        nsp <- ncol(x$mean$Rho)
        ut <- upper.tri(x$mean$Rho)
        sp <- arrayInd(which(ut), c(nsp, nsp))
        ans <- data.frame(
          model = i,
          sp1 = sp[, 1],
          sp2 = sp[, 2],
          rho = c(prob_cooccur_es(x$model$cluster1$data()$Y)[ut], x$mean$Rho[ut]),
          rho_type = rep(c("Effect-Size", "Residual CM"), each = sum(ut)),
          sgn = sign(x$mean$Rho)[ut],
          significant = x$overlap$Rho[ut],
          #cint = simulation_parameters$comp_inter[[i]][ut],
          #fint = simulation_parameters$fac_inter[[i]][ut],
          cint = comp_inter[[i]][ut],
          fint = fac_inter[[i]][ut],
          density = tail(nm, 2)[1],
          type = paste0(head(nm, -2), collapse = ""),
          nsp = nsp,
          stringsAsFactors = FALSE
        )
        ans$cint[is.na(ans$cint)] <- 0
        ans$fint[is.na(ans$fint)] <- 0
        ans$density[ans$density == "Even"] <- "None"
        ans$density <- factor(ans$density, c("None", "Sparse", "Dense"))
        ans$type[ans$type == "Env"] <- "Environmental\nFiltering Only"
        ans$type[ans$type == "Fac"] <- "Facilitation"
        ans$type[ans$type == "Comp"] <- "Competition"
        ans$type[ans$type == "FacComp"] <- "Facililation +\nCompetition"
        ans$type <- factor(
          ans$type,
          c(
            "Environmental\nFiltering Only", "Facilitation", "Competition",
            "Facililation +\nCompetition"
          )
        )
        ans$interaction <- "None"
        ans$interaction <- ifelse(ans$cint, "Competition", ans$interaction)
        ans$interaction <- ifelse(ans$fint, "Facilitation", ans$interaction)
        ans$status <- ifelse(
          ans$significant,
          ifelse(ans$sgn * -ans$cint == 1 | ans$sgn * ans$fint == 1, "TP", "FP"),
          ifelse(ans$cint == 0 & ans$fint == 0, "TN", "FN")
        )
        ans$interaction <- factor(
          ans$interaction, c("None", "Facilitation", "Competition")
        )
        ans
      }
    )
  )
  return(mean_correlations)
}


mean_corr_val<- mean_cor(jsdm_list)


# mean_correlations <- mean_cor(models)
# x <- subset(mean_correlations, type != "Environmental\nFiltering Only")
# acc <- by(x, x$model, function(x) sum(x$status == "TP" | x$status == "TN") / nrow(x))

#' Plot correlation parameter means
#+ plot-correlations
m<- ggplot(mean_corr_val) +
  aes(factor(nsp), rho, fill = interaction) +
  geom_hline(yintercept = 0) +
  geom_boxplot(
    outlier.size = .2, size = .1, position = position_dodge(preserve = "single")
  ) +
  scale_fill_manual(values = c("grey", "blue", "red")) +
  facet_grid(rho_type ~ type + density) +
  xlab("Number of species") +
  ylab("Correlation") +
  theme_bw() +
  theme(legend.position = "top")
m

######################################### HMSC

mean_cor_other<-function(mod,lab){
  mean_correlations <- do.call(
    rbind,
    lapply(
      seq_along(mod),
      function(i) {
        x <- mod[[i]]
        nm <- strsplit(
          names(mod)[[i]], "(?<=[a-z])(?=[A-Z])", perl = TRUE
        )[[1]]
        nsp <- ncol(x$Rho_mean)
        ut <- upper.tri(x$Rho_mean)
        sp <- arrayInd(which(ut), c(nsp, nsp))
        ans <- data.frame(
          model = i,
          sp1 = sp[, 1],
          sp2 = sp[, 2],
          rho = x$Rho_mean[ut],
          rho_type = rep(paste0("Residual ",lab), each = sum(ut)),
          sgn = sign(x$Rho_mean)[ut],
          significant = !(x$Rho_sign[ut]==0),

          cint = comp_inter[[i]][ut],
          fint = fac_inter[[i]][ut],
          density = tail(nm, 2)[1],
          type = paste0(head(nm, -2), collapse = ""),
          nsp = nsp,
          stringsAsFactors = FALSE
        )
        ans$cint[is.na(ans$cint)] <- 0
        ans$fint[is.na(ans$fint)] <- 0
        ans$density[ans$density == "Even"] <- "None"
        ans$density <- factor(ans$density, c("None", "Sparse", "Dense"))
        ans$type[ans$type == "Env"] <- "Environmental\nFiltering Only"
        ans$type[ans$type == "Fac"] <- "Facilitation"
        ans$type[ans$type == "Comp"] <- "Competition"
        ans$type[ans$type == "FacComp"] <- "Facililation +\nCompetition"
        ans$type <- factor(
          ans$type,
          c(
            "Environmental\nFiltering Only", "Facilitation", "Competition",
            "Facililation +\nCompetition"
          )
        )
        ans$interaction <- "None"
        ans$interaction <- ifelse(ans$cint, "Competition", ans$interaction)
        ans$interaction <- ifelse(ans$fint, "Facilitation", ans$interaction)
        ans$status <- ifelse(
          ans$significant,
          ifelse(ans$sgn * -ans$cint == 1 | ans$sgn * ans$fint == 1, "TP", "FP"),
          ifelse(ans$cint == 0 & ans$fint == 0, "TN", "FN")
        )
        ans$interaction <- factor(
          ans$interaction, c("None", "Facilitation", "Competition")
        )
        ans
      }
    )
  )
  return(mean_correlations)
}


hm_inter<-function(mod){
  getOmega = function(a,r=1)
    return(crossprod(a$Lambda[[r]]))
  ns<-mod$ns
  nsamples<-mod$samples
  nchains<-2
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

  return(list(Rho_mean=postRMean,Rho_sign=Toplot_R, Tau=postTMean,Tau_sign=Toplot_T)) 
}

R_list<-lapply(hmsc_list,hm_inter)
Tau_list<-lapply(R_list, function(x) list(Rho_mean=x$Tau,Rho_sign=x$Tau_sign))


hmsc_cor<-mean_cor_other(R_list,lab="HMSC")
hmsc_pcor<-mean_cor_other(Tau_list,lab="HMSC")



#####GJAM#########################################################################################################


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


gjam_Rho<-function(name){
  gj_mod<-load_object(name)
  burn<-gj_mod$m1$modelList$burnin
  J<-ncol(gj_mod$m1$inputs$y)
  gjam_bs<- mcmc.list(mcmc(gj_mod$m1$chains$bgibbsUn[-(1:burn),]),mcmc(gj_mod$m2$chains$bgibbsUn[-(1:burn),]))
  gjam_sigma<- mcmc.list(mcmc(gj_mod$m1$chains$sgibbs[-(1:burn),]),mcmc(gj_mod$m2$chains$sgibbs[-(1:burn),]))
   gjam_mods_2sgibbs<-abind(gj_mod$m1$chains$sgibbs,gj_mod$m2$chains$sgibbs, along = 3)
   postH<- apply(gjam_mods_2sgibbs, 2, quantile,0.95)
   postL<-apply(gjam_mods_2sgibbs, 2, quantile,0.05)
   post_mean<-apply(gjam_mods_2sgibbs, 2, mean)
   pH<-convert_to_m(postH)
   pL<-convert_to_m(postL)
   post_mean_s<-convert_to_m(post_mean)
   S_mean<-cov2cor(post_mean_s)
   R_sign<-cov2cor(post_mean_s)*(!(pH>0 & pL<0))
   sgibbs<-abind(gj_mod$m1$chains$sgibbs[-(1:burn),],gj_mod$m2$chains$sgibbs[-(1:burn),],along=1)
   
   tau<-array(NA,dim=c(J,J,dim(sgibbs)[1]))
   for(j in 1:dim(sgibbs)[1]){
     ss <- expandSigma_rmd(sgibbs[j,], S = J)
     si <- solve(ss)
     tau[,,j] <- -cov2cor(si)
   }
   
   tau_mean<-apply(tau,c(1,2), mean)
   tau_HI<-apply(tau,c(1,2),quantile,0.95)
   tau_LO<-apply(tau,c(1,2),quantile,0.05)
   Tau_sign<-tau_mean*(!(tau_HI>0 & tau_LO<0))
 
 return(list(Rho_mean=S_mean,Rho_sign=R_sign, Tau=tau_mean,Tau_sign=Tau_sign))
}


gjam_list_R<- lapply(gjam_files_list,gjam_Rho)
names(gjam_list_R)<- sim_names
Tau_list_gjam<-lapply(gjam_list_R, function(x) list(Rho_mean=x$Tau,Rho_sign=x$Tau_sign))
gjam_mean_cor<-mean_cor_other(gjam_list_R,lab="GJAM")
gjam_mean_cor_p<-mean_cor_other(Tau_list_gjam,lab="GJAM")



#####GJAM DR#########################################################################################################
gjam_dr_Rho<-function(name){
    gjam_dr_mods<-load_object(name)
    mod_gjam_red_1<-gjam_dr_mods$m1
    mod_gjam_red_2<-gjam_dr_mods$m2
    ns<-ncol(mod_gjam_red_1$inputs$y)
    ng<-mod_gjam_red_2$modelList$ng
    burnin<-mod_gjam_red_2$modelList$burnin  
    gjam_sigma<- mcmc.list(mcmc(mod_gjam_red_1$chains$sgibbs[-(1:burnin),]),mcmc(mod_gjam_red_2$chains$sgibbs[-(1:burnin),]))
    sgibbs<-abind(mod_gjam_red_1$chains$sgibbs[-(1:burnin),],mod_gjam_red_2$chains$sgibbs[-(1:burnin),],along=1)
    sigErrGibbs<-abind(mod_gjam_red_1$chains$sigErrGibbs[-(1:burnin)],mod_gjam_red_2$chains$sigErrGibbs[-(1:burnin)],along=1)
    kgibbs<-abind(mod_gjam_red_1$chains$kgibbs[-(1:burnin),],mod_gjam_red_2$chains$kgibbs[-(1:burnin),],along=1)
    sigma<-invsigma<-array(NA,dim=c(ns,ns,2*(ng-burnin))) 
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
    S_mean<-cov2cor(sigma_mean)
    Sigma_sign<--cov2cor(sigma_mean*(!(sigma_q95>0 & sigma_q05<0)))
    invsigma_mean<-apply(invsigma,c(1,2),mean) 
    invsigma_q05<-apply(invsigma,c(1,2),quantile,0.05) 
    invsigma_q95<-apply(invsigma,c(1,2),quantile,0.95) 
    INVSigma_sign<--cov2cor(sigma_mean*(!(invsigma_q95>0 & invsigma_q05<0)))
    
  return(list(Rho_mean=S_mean,Rho_sign=Sigma_sign, Tau=invsigma_mean,Tau_sign=INVSigma_sign))
}



gjam_dr_list_R<- lapply(gjam_dr_files_list,gjam_dr_Rho)
names(gjam_dr_list_R)<- sim_names
Tau_list_gj_dr<-lapply(gjam_dr_list_R, function(x) list(Rho_mean=x$Tau,Rho_sign=x$Tau_sign))

gjam_dr_mean_cor<-mean_cor_other(gjam_dr_list_R,lab="DR-GJAM")
# 
gjam_dr_mean_cor_p<-mean_cor_other(Tau_list_gj_dr,lab="DR-GJAM")


total_mean_cor<- rbind(mean_corr_val,hmsc_cor,gjam_mean_cor,gjam_dr_mean_cor)

#pdf("plot_mean_corr.pdf")
#png("plot_mean_corr.png")
m<- ggplot(total_mean_cor) +
  aes(factor(nsp), rho, fill = interaction) +
  geom_hline(yintercept = 0) +
  geom_boxplot(
    outlier.size = .2, size = .1, position = position_dodge(preserve = "single")
  ) +
  scale_fill_manual(values = c("grey", "blue", "red")) +
  facet_grid(rho_type ~ type + density) +
  xlab("Number of species") +
  ylab("Correlation") +
  theme_bw() +
  theme(legend.position = "top")
m
#dev.off()
#####Partial correlation

total_mean_part_cor<- rbind(hmsc_pcor,gjam_mean_cor_p,gjam_dr_mean_cor_p)

pdf("plot_mean_part_corr.pdf")
#png("plot_mean_part_corr.png")
m<- ggplot(total_mean_part_cor) +
  aes(factor(nsp), rho, fill = interaction) +
  geom_hline(yintercept = 0) +
  geom_boxplot(
    outlier.size = .2, size = .1, position = position_dodge(preserve = "single")
  ) +
  scale_fill_manual(values = c("grey", "blue", "red")) +
  facet_grid(rho_type ~ type + density) +
  xlab("Number of species") +
  ylab("Partial correlation") +
  theme_bw() +
  theme(legend.position = "top")
m
dev.off()

#############################################################################################################################################################
### Posterior predictive check
### Conditiobal Predictive Ordinate (CPO)

cor_relation[abs(cor_relation) < 0.6] <- NA


