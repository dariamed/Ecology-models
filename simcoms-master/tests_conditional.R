

library(Hmsc)


rm(list=ls())
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
new_sim_data<-readRDS("new_sim_data.rds")


######Conditional prediction test
data_inp<-sim_data$FacCompSparseSp20
#data_inp<-new_sim_data$FacCompSparseSp20
gjam_5<- "./gjam_models/gjam5cmpd.rda"
data_inp5<-sim_data$CompDenseSp5

########################################Functions##################################################################

load_object <- function(file) {
  tmp <- new.env()
  load(file = file, envir = tmp)
  tmp[[ls(tmp)[1]]]
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
gjam_predict_out<- function(datap, name){
  mod<-load_object(name)
  np<-200
  xdata<- as.data.frame(scale(poly(datap$env, 2))[1:np,])
  colnames(xdata)<- c("env", "env2")
  newdata <- list(xdata = xdata, nsim=200)
  predict<-  gjamPredict(mod$m1, newdata = newdata)
  y_full    <- predict$sdList$yMu
  AUC_g<-vector()
  for(i in 1:(ncol(datap)-1)) AUC_g<-c(AUC_g,auc(roc(y_full[,i],factor(datap[1:np,i]))))
  return(list(AUC=AUC_g,pred=y_full) )
}


#########################################################################################################################


data20FCSP<-dataprep(sim_data$FacCompSparseSp20) 
data5d<-dataprep(sim_data$CompDenseSp5) 


xdata20FC<-as.data.frame(data20FCSP$X[,-1])
colnames(xdata20FC)<- c("env","env2")
ydata20FC<-as.data.frame(data20FCSP$Y)
####fit



gjamfc20<-"./gjam_models/gjamcmp_facs20.rda"
gjam_5<-"./gjam_models/gjam5cmpd.rda"
mod5fcd<-load_object(gjam_5)
mod20cd<-load_object(gjamfc20)

# 
# newdata <- list(xdata = xdata, nsim=200)
# predict<-  gjamPredict(out$m1, newdata = newdata)
# plot(data_inp$env,predict$sdList$yMu[,3])


new20fc <- list(ydataCond = mod20cd$m1$inputs$y[,-3], nsim=500)   # cond on obs PA data
p1  <- gjamPredict(output =mod20cd$m1, newdata = new20fc,FULL = TRUE)
plot( mod20cd$m1$inputs$xdata[,1],p1$sdList$yMu[,3])


new5fd <- list(ydataCond = mod5fcd$m1$inputs$y[,-3], nsim=5000)   # cond on obs PA data
p2  <- gjamPredict(output =mod5fcd$m1, newdata = new5fd,FULL = TRUE)
plot( mod5fcd$m1$inputs$xdata[,1],p2$sdList$yMu[,3])


########################################HMSC###################################################################
hm20fc<-"./new_HMmodels/hmcmp_facs20.rda"
hm5cd<- "./new_HMmodels/hm5cmpd.rda"

hm_mod20fcd<-load_object(hm20fc)
hm_mod5cd<-load_object(hm5cd)

X20FS<-scale(poly(sim_data$FacCompSparseSp20$env, 2))
colnames(X20FS)<-c("env","env2")
XDataNew20FS<-as.data.frame(X20FS)
studyDesignNew20FS = data.frame(sample = as.factor(1:nrow(XDataNew20FS)))
rLNew20FS = HmscRandomLevel(units = studyDesignNew20FS$sample)
predY20FS = predict(hm_mod20fcd, XData=XDataNew20FS, studyDesign=studyDesignNew20FS, ranLevels=list(sample = rLNew20FS),Yc=hm_mod20fcd$Y[,1:5])
prY<- apply(simplify2array(predY20FS), 1:2, mean)

y<- predict(hm_mod, Yc=hm_mod$Y[,-(3:10)])






########################################GLM prediction###################################################################


predict_sdm<-function(datas){
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

sdm_env5<-predict_sdm(datas=sim_data$EnvEvenSp5)
















x<- data$X
it<- out$m1$modelList$ng
mu<-array(NA,dim=c(data$n,data$J,it))
for(k in 1:it){
  for(j in 1:data$J){
    
    mu[,j,k] <- pnorm(x%*%out$m1$chains$bgibbs[k,(3*(j-1)+1):(3*j)])
    
  }
}

np<-200
mod<-load_object(name)
X<-scale(poly(datap$env[1:np], 2))
colnames(X)<-c("env","env2")
XDataNew<-as.data.frame(X)
studyDesignNew = data.frame(sample = as.factor(1:np))
rLNew = HmscRandomLevel(units = studyDesignNew$sample)
predY = predict(mod, XData=XDataNew, studyDesign=studyDesignNew, ranLevels=list(sample = rLNew), expected=TRUE )
prY<- apply(simplify2array(predY), 1:2, mean)
AUC_hmcs<-vector()
for(i in 1:(ncol(datap)-1)) AUC_hmcs<-c(AUC_hmcs,auc(roc(prY[,i],factor(datap[1:np,i]))))
return(AUC_hmcs)


new <- list( nsim=5000)   # cond on obs PA data

newdata   <- list(xdata = xdata, nsim = 50 )


p2  <- gjamPredict(output = out$m1)
ydataCond <- out$m1$inputs$y[,1,drop=FALSE]*0
newdata   <- list(ydataCond = ydataCond, nsim=50)
p0        <- gjamPredict(output = out$m1, newdata = newdata)


newdata   <- list(ydataCond = sim_data$FacCompSparseSp20$env, nsim=50)

p0        <- gjamPredict(output = out$m1, newdata = newdata)


######### HMSC

hm_mod<-load_object("./new_HMmodels/hmcmp_facs20.rda")
X<-scale(poly(data_inp$env, 2))
colnames(X)<-c("env","env2")
XDataNew<-as.data.frame(X)
studyDesignNew = data.frame(sample = as.factor(1:nrow(XDataNew)))
rLNew = HmscRandomLevel(units = studyDesignNew$sample)
predY = predict(hm_mod, XData=XDataNew, studyDesign=studyDesignNew, ranLevels=list(sample = rLNew),   Yc=hm_mod$Y)
prY<- apply(simplify2array(predY), 1:2, mean)

y<- predict(hm_mod, Yc=hm_mod$Y[,-(3:10)])





