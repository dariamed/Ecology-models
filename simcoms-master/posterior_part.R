#TO INSTALL mcclust.ext
#devtools::install_url("http://wrap.warwick.ac.uk/71934/1/mcclust.ext_1.0.tar.gz")

library(knitr)
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
library(reshape)
library(mcclust.ext)

setwd("~/Documents/GitHub/Ecology-models/simcoms-master")
sim_data<-readRDS("sim_data.rds")
data_20<-sim_data$FacDenseSp20

ydata<-data_20[,-21]
xdata_20<-scale(poly(data_20$env, 2))
colnames(xdata_20)<- c("env","env2")
formula<- ~env + env2
rl   <- list(r = 5, N = 10)
ml   <- list(ng = 5000, burnin = 500, typeNames = 'PA', reductList = rl, PREDICTX = F )
out  <- gjam(formula, xdata = xdata_20, ydata = ydata, modelList = ml)

k_unique_clust<- apply(out$chains$kgibbs,1,unique)
k_unique_num<- lapply(k_unique_clust,length)
hist(unlist(k_unique_num[1:2500]))



###Demo
credibleball(c.star, cls.draw, c.dist = c("VI","Binder"), alpha = 0.05)
credibleball(c.star, out$chains$kgibbs, c.dist = c("VI","Binder"), alpha = 0.05)
data(galaxy.fit)
x=data.frame(x=galaxy.fit$x)
data(galaxy.pred)
data(galaxy.draw)
psm=comp.psm(galaxy.draw)
galaxy.VI=minVI(psm,galaxy.draw,method=("all"),include.greedy=TRUE)
summary(galaxy.VI)
plot(galaxy.VI,data=x,dx=galaxy.fit$fx,xgrid=galaxy.pred$x,dxgrid=galaxy.pred$fx)
VI(galaxy.VI$cl,galaxy.draw)
# Binder
galaxy.B=minbinder.ext(psm,galaxy.draw,method=("all"),include.greedy=TRUE)
summary(galaxy.B)

plot(galaxy.B,data=x,dx=galaxy.fit$fx,xgrid=galaxy.pred$x,dxgrid=galaxy.pred$fx)
# Uncertainty in partition estimate
galaxy.cb=credibleball(galaxy.VI$cl[1,],galaxy.draw)
summary(galaxy.cb)
plot(galaxy.cb,data=x,dx=galaxy.fit$fx,xgrid=galaxy.pred$x,dxgrid=galaxy.pred$fx)
# Compare with uncertainty in heat map of posterior similarity matrix
data(ex2.data)
x=data.frame(ex2.data[,c(1,2)])
cls.true=ex2.data$cls.true
plot(x[,1],x[,2],xlab="x1",ylab="x2")
k=max(cls.true)
for(l in 2:k){
  points(x[cls.true==l,1],x[cls.true==l,2],col=l)
  }
# Find representative partition of posterior
data(ex2.draw)
psm=comp.psm(ex2.draw)
ex2.VI=minVI(psm,ex2.draw,method=("all"),include.greedy=TRUE)
summary(ex2.VI)
plot(ex2.VI,data=x)

##########K
psm_k=comp.psm(out$chains$kgibbs)
k.VI=minVI(psm_k,out$chains$kgibbs[1:1000,],method=("all"))
summary(k.VI)

corrplot(psm_k, diag = FALSE, order = "original",tl.pos = "ld", tl.cex = 0.5, method = "color",col=cols(200), type = "lower")



######### Library I.Grpah
