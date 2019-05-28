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
new_sim_data<-readRDS("sim_data.rds")



######################Testing the variance in occurence of species########################################

run_data<- sim_data$EnvEvenSp10

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
      labs(title=paste0("Fundamental niche and realized points, Species ",i, " for "," for environment"))+
      scale_color_manual(name = c("Legend"), values = c("Predicted probability" = "#FF6666","Fundamental niche"="#9999FF"))+theme_bw()+ theme(plot.title = element_text(hjust = 0.5))
    p_0<<-list.append(p_0,g)
    assign(paste0("p",i), g, pos =1)
  })


p_0[[7]]



run_data<- sim_data$FacSparseSp10

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

p_f<- list()
for(i in 1:nspecies)
  local({
    i<-i
    tmp_obs<-table_obs[which(table_obs$species==i),]
    tmp_fund<-table_fundamental[which(table_fundamental$species==i),]
    g<<-ggplot()+
      geom_point(aes(x=tmp_obs$xx,y=tmp_obs$obs),col="#000066",size=0.5) +xlab("Environmental gradient")+ylab("Probability of presence")+
      geom_line(data=tmp_fund,aes(x=tmp_fund$xx,y=tmp_fund$niche,color = "Fundamental niche"),lwd=1)+
      labs(title=paste0("Fundamental niche and realized points, Species ",i, " for "," for facitilation"))+
      scale_color_manual(name = c("Legend"), values = c("Predicted probability" = "#FF6666","Fundamental niche"="#9999FF"))+theme_bw()+ theme(plot.title = element_text(hjust = 0.5))
    p_f<<-list.append(p_f,g)
    assign(paste0("p",i), g, pos =1)
  })


p_f[[7]]


#################new_sim_data################################################################


run_data<- new_sim_data$EnvEvenSp10

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

p_02<- list()
for(i in 1:nspecies)
  local({
    i<-i
    tmp_obs<-table_obs[which(table_obs$species==i),]
    tmp_fund<-table_fundamental[which(table_fundamental$species==i),]
    g<<-ggplot()+
      geom_point(aes(x=tmp_obs$xx,y=tmp_obs$obs),col="#000066",size=0.5) +xlab("Environmental gradient")+ylab("Probability of presence")+
      geom_line(data=tmp_fund,aes(x=tmp_fund$xx,y=tmp_fund$niche,color = "Fundamental niche"),lwd=1)+
      labs(title=paste0("Fundamental niche and realized points, Species ",i, " for "," for environment"))+
      scale_color_manual(name = c("Legend"), values = c("Predicted probability" = "#FF6666","Fundamental niche"="#9999FF"))+theme_bw()+ theme(plot.title = element_text(hjust = 0.5))
    p_02<<-list.append(p_02,g)
    assign(paste0("p",i), g, pos =1)
  })


p_02[[7]]



run_data<- new_sim_data$FacSparseSp10

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

p_f2<- list()
for(i in 1:nspecies)
  local({
    i<-i
    tmp_obs<-table_obs[which(table_obs$species==i),]
    tmp_fund<-table_fundamental[which(table_fundamental$species==i),]
    g<<-ggplot()+
      geom_point(aes(x=tmp_obs$xx,y=tmp_obs$obs),col="#000066",size=0.5) +xlab("Environmental gradient")+ylab("Probability of presence")+
      geom_line(data=tmp_fund,aes(x=tmp_fund$xx,y=tmp_fund$niche,color = "Fundamental niche"),lwd=1)+
      labs(title=paste0("Fundamental niche and realized points, Species ",i, " for "," for facitilation"))+
      scale_color_manual(name = c("Legend"), values = c("Predicted probability" = "#FF6666","Fundamental niche"="#9999FF"))+theme_bw()+ theme(plot.title = element_text(hjust = 0.5))
    p_f2<<-list.append(p_f2,g)
    assign(paste0("p",i), g, pos =1)
  })


p_f2[[7]]












