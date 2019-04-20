
##########
library(knitr)
#library(arm)
#library(jagsUI)
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

setwd("~/Desktop/VirtualCommunity/simcoms-master/ExampleFiles")
# setwd("~/Documents/GitHub/Ecology-models/simcoms-master/ExampleFiles")
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


setwd("~/Desktop/VirtualCommunity/simcoms-master/")
sim_data<-readRDS("sim_data.rds")

gjam_fit<-function(data, it=2500,burn=500 , name="./gjam_models/temp.rda",interact=diag(ncol(data$Y))){
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
  ml   <- list(ng = it, burnin = burn, typeNames = 'PA')
  ####fit
  mod_gjam1  <- gjam(formula, xdata = xdata, ydata = ydata, modelList = ml)
  save(mod_gjam1, file = name)
}

hmsc_fit<-function(data,nsamples = 1000,nchains=2,name="./HMmodels/hmtemp.rda" ){
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
}
###########################################################################Environment
setwd("~/Desktop/VirtualCommunity/simcoms-master/new_models")
#5
data5<-sim_data$EnvEvenSp5
gjam_fit(data5,5000,500,"./gjam_models/gjam5env.rda")
hmsc_fit(data5,name="./HMmodels/hm5env.rda" )  
#10
data10<-sim_data$EnvEvenSp10
gjam_fit(data10,5000,500,"./gjam_models/gjam10env.rda")
hmsc_fit(data10,name="./HMmodels/hm10env.rda" )  
#20
data20<-sim_data$EnvEvenSp20
gjam_fit(data20,5000,500,"./gjam_models/gjam20env.rda")
hmsc_fit(data20,name="./HMmodels/hm20env.rda" )  
####################################################################################### 
###########################################################################Environment + Facilitation Dense
#5
data5<-sim_data$FacDenseSp5
gjam_fit(data5,5000,500,"./gjam_models/gjam5facd.rda")
hmsc_fit(data5,name="./HMmodels/hm5facd.rda" )  
#10
data10<-sim_data$FacDenseSp10
gjam_fit(data10,5000,500,"./gjam_models/gjam10facd.rda")
hmsc_fit(data10,name="./HMmodels/hm10facd.rda" )  
#20
data20<-sim_data$FacDenseSp20
gjam_fit(data20,5000,500,"./gjam_models/gjam20facd.rda")
hmsc_fit(data20,name="./HMmodels/hm20facd.rda" )  

####################################################################################### 
###########################################################################Environment + Facilitation Sparse
#5
data5<-sim_data$FacSparseSp5
gjam_fit(data5,5000,500,"./gjam_models/gjam5facs.rda")
hmsc_fit(data5,name="./HMmodels/hm5facs.rda" )  
#10
data10<-sim_data$FacSparseSp10
gjam_fit(data10,5000,500,"./gjam_models/gjam10facs.rda")
hmsc_fit(data10,name="./HMmodels/hm10facs.rda" )  
#20
data20<-sim_data$FacSparseSp20
gjam_fit(data20,5000,500,"./gjam_models/gjam20facs.rda")
hmsc_fit(data20,name="./HMmodels/hm20facs.rda" )  

####################################################################################### 
###########################################################################Environment + Facilitation Sparse
#5
data5<-sim_data$FacSparseSp5
gjam_fit(data5,5000,500,"./gjam_models/gjam5facs.rda")
hmsc_fit(data5,name="./HMmodels/hm5facs.rda" )  
#10
data10<-sim_data$FacSparseSp10
gjam_fit(data10,5000,500,"./gjam_models/gjam10facs.rda")
hmsc_fit(data10,name="./HMmodels/hm10facs.rda" )  
#20
data20<-sim_data$FacSparseSp20
gjam_fit(data20,5000,500,"./gjam_models/gjam20facs.rda")
hmsc_fit(data20,name="./HMmodels/hm20facs.rda" )  

####################################################################################### 
###########################################################################Environment + Competition Dense
#5
data5<-sim_data$CompDenseSp5
gjam_fit(data5,5000,500,"./gjam_models/gjam5compd.rda")
hmsc_fit(data5,name="./HMmodels/hm5compd.rda" )  
#10
data10<-sim_data$CompDenseSp10
gjam_fit(data10,5000,500,"./gjam_models/gjam10compd.rda")
hmsc_fit(data10,name="./HMmodels/hm10compd.rda" )  
#20
data20<-sim_data$CompDenseSp20
gjam_fit(data20,5000,500,"./gjam_models/gjam20compd.rda")
hmsc_fit(data20,name="./HMmodels/hm20compd.rda" )  

####################################################################################### 
###########################################################################Environment + Competition Sparse
#5
data5<-sim_data$CompSparseSp5
gjam_fit(data5,5000,500,"./gjam_models/gjam5comps.rda")
hmsc_fit(data5,name="./HMmodels/hm5comps.rda" )  
#10
data10<-sim_data$CompSparseSp10
gjam_fit(data10,5000,500,"./gjam_models/gjam10comps.rda")
hmsc_fit(data10,name="./HMmodels/hm10comps.rda" )  
#20
data20<-sim_data$CompSparseSp20
gjam_fit(data20,5000,500,"./gjam_models/gjam20comps.rda")
hmsc_fit(data20,name="./HMmodels/hm20comps.rda" )  

####################################################################################### 
###########################################################################Environment + Competition  + Facilitation Dense
#5
data5<-sim_data$FacCompDenseSp5
gjam_fit(data5,5000,500,"./gjam_models/gjam5comp_facd.rda")
hmsc_fit(data5,name="./HMmodels/hm5comp_facd.rda" )  
#10
data10<-sim_data$FacCompDenseSp10
gjam_fit(data10,5000,500,"./gjam_models/gjam10comp_facd.rda")
hmsc_fit(data10,name="./HMmodels/hm10comp_facd.rda" )  
#20
data20<-sim_data$FacCompDenseSp20
gjam_fit(data20,5000,500,"./gjam_models/gjam20comp_facd.rda")
hmsc_fit(data20,name="./HMmodels/hm20comp_facd.rda" )  

###########################################################################Environment + Competition  + Facilitation Sparse
#5
data5<-sim_data$FacCompSparseSp5
gjam_fit(data5,5000,500,"./gjam_models/gjam5comp_facs.rda")
hmsc_fit(data5,name="./HMmodels/hm5comp_facs.rda" )  
#10
data10<-sim_data$FacCompSparseSp10
gjam_fit(data10,5000,500,"./gjam_models/gjam10comp_facs.rda")
hmsc_fit(data10,name="./HMmodels/hm10comp_facs.rda" )  
#20
data20<-sim_data$FacCompSparseSp20
gjam_fit(data20,5000,500,"./gjam_models/gjam20comp_facs.rda")
hmsc_fit(data20,name="./HMmodels/hm20comp_facs.rda" )  


####################################################################################### 
