

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

jsdm_files_list<-list(env5="model-2019-04-09-19-02-16.rda",env10="model-2019-04-10-08-26-20.rda",env20="model-2019-04-11-19-06-02.rda",
                      facd5="model-jsdm-2019-05-06-18-41-20.rda",facd10="model-2019-04-12-06-59-42.rda",facd20="model-2019-04-13-17-22-12.rda",
                      facs5="model-2019-04-13-17-51-16.rda",facs10="model-2019-04-14-05-17-58.rda",facs20="model-2019-04-15-15-34-20.rda",
                      compd5="model-2019-04-15-16-03-39.rda",compd10="model-2019-04-16-03-39-17.rda",compd20="model-2019-04-17-14-00-45.rda",
                      comps5="model-2019-04-17-14-30-01.rda",comps10="model-2019-04-18-01-55-45.rda",comps20="model-2019-04-19-12-40-30.rda"
                      ,compfacd5="model-2019-04-19-13-08-47.rda",compfacd10="model-2019-04-20-00-37-25.rda",compfacd20="model-2019-04-21-11-10-34.rda"
                      , compfacs5="model-2019-04-21-11-39-52.rda",compfacs10="model-2019-04-21-22-58-58.rda",compfacs20="model-2019-04-23-09-42-18.rda")

jsdm_list<- lapply(jsdm_files_list, load_object)


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
  return(list(neff=neff_mod,Rhat=Rhat))
}

###########################################HMSC functions#######################################################################

hmsc_files_list<-list(env5="./HMmodels/hm5env.rda",env10="./HMmodels/hm10env.rda" ,env20="./HMmodels/hm20env.rda" 
                      ,facd5="./HMmodels/hm5fd.rda" ,facd10="./HMmodels/hm10fd.rda",facd20="./HMmodels/hm20fd.rda",
                      facs5="./HMmodels/hm5fs.rda",facs10="./HMmodels/hm10fs.rda",facs20="./HMmodels/hm20fs.rda" ,
                      compd5="./HMmodels/hm5cmpd.rda",compd10="./HMmodels/hm10cmpd.rda" ,compd20="./HMmodels/hm20cmpd.rda" ,
                      comps5="./HMmodels/hm5cmps.rda",comps10="./HMmodels/hm10cmps.rda",comps20="./HMmodels/hm20cmps.rda",
                      compfacd5="./HMmodels/hmcmp_facd5.rda",compfacd10="./HMmodels/hmcmp_facd10.rda" ,compfacd20="./HMmodels/hmcmp_facd20.rda" ,
                      compfacs5="./HMmodels/hmcmp_facs5.rda",compfacs10="./HMmodels/hmcmp_facs10.rda" ,compfacs20="./HMmodels/hmcmp_facs20.rda" )


hmsc_list<- lapply(hmsc_files_list, load_object)


hmsc_convergence_par<- lapply(hmsc_list, hm_conv)

hmsc_conv_dataset<- data.frame()

for(j in 1:length(hmsc_list)){
  tmp1<- as.data.frame(hmsc_convergence_par[[j]][1])
  colnames(tmp1)<-c("parameter", "effective_size", "value")
  tmp1<-tmp1[,-2]
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




gjam_files_list<-list(env5="./gjam_models/gjam5env.rda",env10="./gjam_models/gjam10env.rda",env20="./gjam_models/gjam20env.rda",
                      facd5="./gjam_models/gjam5f.rda",facd10="./gjam_models/gjam10fd.rda",facd20="./gjam_models/gjam20fd.rda",
                      facs5="./gjam_models/gjam5fs.rda",facs10="./gjam_models/gjam10fs.rda",facs20="./gjam_models/gjam20fs.rda",
                      compd5="./gjam_models/gjam5cmpd.rda",compd10="./gjam_models/gjam10cmpd.rda",compd20="./gjam_models/gjam20cmpd.rda",
                      comps5="./gjam_models/gjam5cmps.rda",comps10="./gjam_models/gjam10cmps.rda",comps20="./gjam_models/gjam20cmps.rda",
                      compfacd5="./gjam_models/gjam5cmp_facd2.rda",compfacd10="./gjam_models/gjamcmp_facd10.rda",compfacd20="./gjam_models/gjamcmp_facd20.rda",
                      compfacs5="./gjam_models/gjamcmp_facs5.rda",compfacs10="./gjam_models/gjamcmp_facs10.rda",compfacs20="./gjam_models/gjamcmp_facs20.rda")



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




