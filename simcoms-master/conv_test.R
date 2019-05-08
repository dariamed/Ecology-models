
library(gjam)
library(coda)
library(fitR)
library(ggmcmc)
#to recreate sigma
Rcpp::sourceCpp('/Users/dariabystrova/Documents/GitHub/gjamed/src/cppFns.cpp')
source("/Users/dariabystrova/Documents/GitHub/gjamed/R/gjamHfunctions_mod.R")

#setwd("~/Tesi/Code/Ecology-models-master_1/simcoms-master")
setwd("~/Documents/GitHub/Ecology-models/simcoms-master/ExampleFiles")
# lapply(list.files(path = "."),load,.GlobalEnv)

#setwd("~/Tesi/Code/Ecology-models-master/simcoms-master")
load("params.rds")
load("sim_names.rds")
load("comp_inter.rds")
load("fac_inter.rds")
setwd("~/Documents/GitHub/Ecology-models/simcoms-master")
sim_data<-readRDS("sim_data.rds")

ng <- 10000
burnin <- 1000
data_10<- sim_data$EnvEvenSp10
ns<-10
ydata<-data_10[,-11]
xdata<-scale(poly(data_10$env, 2))
r<-9
rl <-list(N=ns-1, r=r)
ml  <- list(ng = ng, burnin = burnin, typeNames = 'PA',reductList=rl)

colnames(xdata)<- c("env","env2")
formula<- ~env + env2
mod_gjam_low <- gjam(formula, xdata, ydata, modelList = ml)


## to set thininhg - choose by =..
ind<-seq(1,dim(mod_gjam_low$chains$sgibbs)[1], by=1)


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


.expandSigmaChains(snames, sgibbs, otherpar, simIndex=simIndex,
                   sigErrGibbs, kgibbs, REDUCT)

sigma<-invsigma<-array(NA,dim=c(ns,ns,ng))

sgibbs<-mod_gjam_low$chains$sgibbs
sigErrGibbs<-mod_gjam_low$chains$sigErrGibbs
kgibbs<-mod_gjam_low$chains$kgibbs
N<-mod_gjam_low$modelList$reductList$N
r<-mod_gjam_low$modelList$reductList$r

for(j in 1:ng){
    Z  <- matrix(sgibbs[j,],N,r)
    sigma[,,j] <- .expandSigma(sigErrGibbs[j], ns, Z = Z, kgibbs[j,], REDUCT = T) #sigma
    invsigma[,,j] <- invWbyRcpp(sigErrGibbs[j], Z[kgibbs[j,],]) #inverse sigma
} 

sigma_mean<-apply(sigma,c(1,2),mean) 
sigma_q05<-apply(sigma,c(1,2),quantile,0.05) 
sigma_q95<-apply(sigma,c(1,2),quantile,0.95) 
Sigma_sign<--cov2cor(sigma_mean*(!(sigma_q95>0 & sigma_q05<0)))
corrplot(Sigma_sign, diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")

.expandSigma <- function(sigma, S, Z = NULL, K = NULL, REDUCT = F){
  
  if(REDUCT) return( sigma*diag(S) + tcrossprod(Z[K,]) )
  
  ss <- diag(S)
  ss[lower.tri(ss,diag=T)] <- sigma
  ss[upper.tri(ss)] <- t(ss)[upper.tri(ss)]
  ss
}


#########TEST for the convergence in real data

library(repmis)
d <- "https://github.com/jimclarkatduke/gjam/blob/master/forestTraits.RData?raw=True"
source_data(d)
xdata <- forestTraits$xdata[,c(1,2,8)]
y  <- gjamReZero(forestTraits$treesDeZero)  # extract y
treeYdata  <- gjamTrimY(y,10)$y             # at least 10 plots
dim(treeYdata)
#treeYdata[1:5,1:6]

rl   <- list(r = 8, N = 40)
ml   <- list(ng = 2500, burnin = 500, typeNames = 'DA', reductList = rl)
form <- as.formula( ~ temp*deficit + I(temp^2) + I(deficit^2) )
out  <- gjam(form, xdata = xdata, ydata = treeYdata, modelList = ml)
# 


########SIMPLE DATA SET
xdata_short<- xdata[1:500,1:2]
ydata_short<-treeYdata[1:500,1:5]
ml   <- list(ng = 25000, burnin = 5000, typeNames = 'DA')
form <- as.formula( ~ temp+ deficit)
mod_s  <- gjam(form, xdata = xdata_short, ydata = ydata_short, modelList = ml)

###############################Convergence#########################################

ind<-seq(1,dim(mod_s$chains$sgibbs)[1], by=1)
gjam_mc<- mcmc(mod_s$chains$sgibbs[ind,]) #thinned sigma
s2s1<- mcmc(mod_s$chains$sgibbs[ind,2]) #s2s1 chain
ind_d<-vector()
ind_d[1]<-1
for (i in 1:9) ind_d[i+1]<- ind_d[i]+ns-i+1
gjam_mc_nd<- mcmc(mod_s$chains$sgibbs[ind,-ind_d]) #thinned sigma non diagonal elements

#acf
#acfplot(gjam_mc)  ##autocor plot
ggs_autocorrelation(ggs(gjam_mc)) ###autocr


#traceplots
x11()


#nESS crazy low
xyplot(s2s1)
hist(effectiveSize(gjam_mc), main="ess(sigma)",lwd=2,col=gray(.6),breaks=100)
hist(effectiveSize(gjam_mc_nd), main="ess(sigma) only non diagonal terms",lwd=2,col=gray(.6))

library(mcclust.ext)

psm_k=comp.psm(mod_gjam_low$chains$kgibbs)

corrplot(psm_k, diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")


##########FOR ANY


######Clustering function
datab<-sim_data$EnvEvenSp5
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
datab<-data


gjam_dim_red<-function(datab,ng=2500, burnin=500,r=1){
  ydata<-subset(datab, select = -env)
  xdata<-scale(poly(datab$env, 2))
  ns<- ncol(ydata)
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
      ml  <- list(ng = ng, burnin = burnin, typeNames = 'PA',reductList=rl)
      mod_gjam_c <- gjam(formula, xdata, ydata, modelList = ml)
      DIC_array[j-1]<-mod_gjam_c$fit$DIC
      RSMPE_array[j-1]<-mod_gjam_c$fit$rmspeAll
    }
    r<- rseq[which.min(RSMPE_array)+1]
    print(DIC_array)
  }
  rl <-list(N=N, r=r)
  ml  <- list(ng = ng, burnin = burnin, typeNames = 'PA',reductList=rl)
  mod_gjam_red <- gjam(formula, xdata, ydata, modelList = ml)
  #.expandSigmaChains(snames, sgibbs, otherpar, simIndex=simIndex, sigErrGibbs, kgibbs, REDUCT)
  sigma<-invsigma<-array(NA,dim=c(ns,ns,ng-burnin))
  sgibbs<-mod_gjam_red$chains$sgibbs
  sigErrGibbs<-mod_gjam_red$chains$sigErrGibbs
  kgibbs<-mod_gjam_red$chains$kgibbs
  N<-mod_gjam_red$modelList$reductList$N
  r<-mod_gjam_red$modelList$reductList$r
  for(j in (1:(ng-burnin))){
   
    Z  <- matrix(sgibbs[j+burnin,],N,r)
    sigma[,,j] <- .expandSigma(sigErrGibbs[j+burnin], ns, Z = Z, kgibbs[j+burnin,], REDUCT = T) #sigma
    invsigma[,,j] <- invWbyRcpp(sigErrGibbs[j+burnin], Z[kgibbs[j+burnin,],]) #inverse sigma
  } 
  
  sigma_mean<-apply(sigma,c(1,2),mean) 
  sigma_q05<-apply(sigma,c(1,2),quantile,0.05) 
  sigma_q95<-apply(sigma,c(1,2),quantile,0.95) 
  Sigma_sign<--cov2cor(sigma_mean*(!(sigma_q95>0 & sigma_q05<0)))
  
  invsigma_mean<-apply(invsigma,c(1,2),mean) 
  invsigma_q05<-apply(invsigma,c(1,2),quantile,0.05) 
  invsigma_q95<-apply(invsigma,c(1,2),quantile,0.95) 
  INVSigma_sign<--cov2cor(sigma_mean*(!(invsigma_q95>0 & invsigma_q05<0)))
  gjam_mc_nd<- mcmc(mod_gjam_red$chains$sgibbs) #thinned sigma non diagonal elements
  hist(effectiveSize(gjam_mc), main="ess(sigma)",lwd=2,col=gray(.6),breaks=100)
  return(list(Rho_sign_d=Sigma_sign,Tau_sign_d=INVSigma_sign,r_opt=r, k=mod_gjam_red$chains$kgibbs ))
}

S<- gjam_dim_red(datab=sim_data$EnvEvenSp10,ng=2500, burnin=500,r=7)
corrplot(S$Rho_sign_d, diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")





##################FOR 5






gjam_dim_red_5<-function(datab,ng=2500, burnin=500){
  ydata<-subset(datab, select = -env)
  xdata<-scale(poly(datab$env, 2))
  ns<- ncol(ydata)
  N<-ncol(ydata)-1
  colnames(xdata)<- c("env","env2")
  formula<- ~env + env2
  r<-3
  rl <-list(N=N, r=r)
  ml  <- list(ng = ng, burnin = burnin, typeNames = 'PA',reductList=rl)
  mod_gjam_red <- gjam(formula, xdata, ydata, modelList = ml)
   sigma<-invsigma<-array(NA,dim=c(ns,ns,ng))
  sgibbs<-mod_gjam_red$chains$sgibbs
  sigErrGibbs<-mod_gjam_red$chains$sigErrGibbs
  kgibbs<-mod_gjam_red$chains$kgibbs
  N<-mod_gjam_red$modelList$reductList$N
  r<-mod_gjam_red$modelList$reductList$r
  for(j in 1:ng){
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
  gjam_mc_nd<- mcmc(mod_gjam_c$chains$sgibbs) #thinned sigma non diagonal elements
  hist(effectiveSize(gjam_mc), main="ess(sigma)",lwd=2,col=gray(.6),breaks=100)
  return(list(Rho_sign_d=Sigma_sign,Tau_sign_d=INVSigma_sign,r_opt=r))
}

S<- gjam_dim_red_5(sim_data$EnvEvenSp5,ng=2500, burnin=500)
corrplot(S$Rho_sign_d, diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")
