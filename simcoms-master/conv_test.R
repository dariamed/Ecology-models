
library(gjam)
library(coda)
library(fitR)
ng <- 10000
burnin <- 1000
data_10<- sim_data$EnvEvenSp10
ns<-10
ydata<-data_10[,-11]
xdata<-scale(poly(data_10$env, 2))
ml  <- list(ng = ng, burnin = burnin, typeNames = 'PA')
colnames(xdata)<- c("env","env2")
formula<- ~env + env2
mod_gjam_low <- gjam(formula, xdata, ydata, modelList = ml)
## to set thininhg - choose by =..
ind<-seq(1,dim(mod_gjam_low$chains$sgibbs)[1], by=10)


gjam_mc<- mcmc(mod_gjam_low$chains$sgibbs[ind,]) #thinned sigma
s2s1<- mcmc(mod_gjam_low$chains$sgibbs[ind,2]) #s2s1 chain
ind_d<-vector()
ind_d[1]<-1
for (i in 1:9) ind_d[i+1]<- ind_d[i]+ns-i+1
gjam_mc_nd<- mcmc(mod_gjam_low$chains$sgibbs[ind,-ind_d]) #thinned sigma non diagonal elements

#acf
#acfplot(gjam_mc)  ##autocor plot
ggs_autocorrelation(ggs(gjam_mc)) ###autocr

#traceplots
x11()
plot(gjam_mc) #traceplots

#nESS crazy low
xyplot(s2s1)
hist(effectiveSize(gjam_mc), main="ess(sigma)",lwd=2,col=gray(.6),breaks=100)
hist(effectiveSize(gjam_mc_nd), main="ess(sigma) only non diagonal terms",lwd=2,col=gray(.6))

#cumuplot
cumuplot(s2s1)
x11()
cumuplot(gjam_mc)
dev.off()

#plotESSBurn(gjam_mc)  ### describe how much should be done the burnin


# BETA

beta_mcmc<-mcmc(mod_gjam_low$chains$bgibbs)
ggs_autocorrelation(ggs(beta_mcmc)) ###autocr

hist(effectiveSize(beta_mcmc), main="ess(beta)",lwd=2,col=gray(.6))
x11()
plot(beta_mcmc)
dev.off()

x11()
cumuplot(beta_mcmc)
dev.off()

