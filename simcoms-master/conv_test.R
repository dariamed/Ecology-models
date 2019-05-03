
library(gjam)
library(coda)
library(fitR)
ng <- 10000
burnin <- 1000
data_10<- sim_data$EnvEvenSp10
ydata<-data_10[,-11]
xdata<-scale(poly(data_10$env, 2))
ml  <- list(ng = ng, burnin = burnin, typeNames = 'PA')
colnames(xdata)<- c("env","env2")
formula<- ~env + env2
mod_gjam_low <- gjam(formula, xdata, ydata, modelList = ml)

ind<-seq(1,dim(mod_gjam_low$chains$sgibbs)[1], by=1)
gjam_mc<- mcmc(mod_gjam_low$chains$sgibbs[ind,])
s2s1<- mcmc(mod_gjam_low$chains$sgibbs[ind,2])
acfplot(gjam_mc)  ##autocor plot
ggs_autocorrelation(ggs(gjam_mc)) ###autocr

plotESSBurn(gjam_mc)  ### describe how much should be done the burnin


hist(effectiveSize(mod_gjam1$chains$sgibbs), main="ess(beta)",lwd=2,col=gray(.6))
hist(effectiveSize(mod_gjam1$chains$bgibbs), main="ess(beta)",lwd=2,col=gray(.6))
# 
#JSDM

# me_mcmc<- mcmc(me5$sims.list$Rho[,4,1])
# xyplot(me_mcmc)
# 

