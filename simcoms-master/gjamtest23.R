

library(repmis)
d <- "https://github.com/jimclarkatduke/gjam/blob/master/forestTraits.RData?raw=True"
source_data(d)
xdata <- forestTraits$xdata[,c(1,2,8)]
y  <- gjamReZero(forestTraits$treesDeZero)  # extract y
treeYdata  <- gjamTrimY(y,10)$y             # at least 10 plots
dim(treeYdata)
#treeYdata[1:5,1:6]
burnin<- 2000
rl   <- list(r = 3, N = 4)
ml   <- list(ng = 10000, burnin = 2000, typeNames = 'DA', reductList = rl, PREDICTX = F )
form <- as.formula( ~ temp*deficit + I(temp^2) + I(deficit^2) )
out  <- gjam(form, xdata = xdata[1:500,], ydata = treeYdata[1:500,1:5], modelList = ml)

gjam_bs<- mcmc.list(mcmc(out$chains$bgibbsUn[-(1:burnin),]))
gjam_sigma<- mcmc.list(mcmc(out$chains$sgibbs[-(1:burnin),]))
n_eff_beta<- as.data.frame(effectiveSize(gjam_bs))
colnames(n_eff_beta)<- c("value")
n_eff_beta$parameter<- "beta"
n_eff_sigma<- as.data.frame(effectiveSize(gjam_sigma))
colnames(n_eff_sigma)<- c("value")
n_eff_sigma$parameter<- "correlation"
neff<- rbind(n_eff_beta,n_eff_sigma)
p<- ggplot(neff, aes(x=value, color=parameter,fill=parameter)) +
  geom_histogram( alpha=0.4, position="identity",binwidth = 100) + xlab("effective size") +
  scale_color_brewer(palette="Dark2")+scale_fill_brewer(palette="Dark2")+
  ggtitle("Effective size for the parameters beta and coefficients of correlation matrix")
plot(p)
Rhat_beta<-as.data.frame(gelman.diag(gjam_bs, multivariate=FALSE)$psrf)
Rhat_beta$parameter<- "beta"
Rhat_sigma<- as.data.frame(gelman.diag(gjam_sigma, multivariate=FALSE)$psrf)
Rhat_sigma$parameter<- "correlation"
Rhat<- rbind(Rhat_beta,Rhat_sigma)
p2<- ggplot(Rhat, aes(x= Rhat$`Point est.`, color=parameter,fill=parameter)) +
  geom_histogram( alpha=0.4, position="identity",binwidth =0.1) +
  scale_color_brewer(palette="Dark2")+scale_fill_brewer(palette="Dark2")+ xlab("Rhat") +
  ggtitle("Rhat for the parameters beta and the coefficients of correlation matrix")
plot(p2)



length(unique(out$chains$kgibbs[1000,]))

##########Clustering
########################################################"

dissimilarity <- 1 - cor(out$parameters$corMu)
distance <- as.dist(dissimilarity) 


plot(hclust(distance),
     main="Dissimilarity = 1 - Correlation", xlab="")

########################################################"

dissimilarity2 <- (1 - cor(out$parameters$corMu))/2
distance2 <- as.dist(dissimilarity2) 


plot(hclust(distance2),
     main="Dissimilarity = (1 - Correlation)/2", xlab="")


########################################################"

dissimilarity3 <- 1 - abs(cor(out$parameters$corMu))
distance3 <- as.dist(dissimilarity3) 


plot(hclust(distance3),
     main="Dissimilarity = 1 - abs(Correlation)", xlab="")


########################################################"

dissimilarity4 <- 1 - abs(cor(out$parameters$corMu))
distance4 <- as.dist(dissimilarity4) 


plot(hclust(distance4),
     main="Dissimilarity = sqrt(1 - (Correlation)^2)", xlab="")



########################################################"




