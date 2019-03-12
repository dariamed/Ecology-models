library(gjam)
library(ggplot2)
library(gridExtra)
#by default n=1000,S=10, Q=5
#f<- gjamSimData(n=500,S=10,Q=5, typeNames = "CA")
#summary(f)

#f$formula
#f$xdata

#x3_mean<-mean(f$xdata$x3)
#x4_mean<-mean(f$xdata$x4)
#x5_mean<-mean(f$xdata$x5)
#max_x2<-max(f$xdata$x2)
#min_x2<-min(f$xdata$x2)


## Example:
S <-5
f <- gjamSimData(n = 200, S = S, Q = 4, typeNames = 'PA')
ml <- list(ng = 50, burnin = 5, typeNames = f$typeNames, holdoutN = 10) 
out <- gjam(f$formula, f$xdata, f$ydata, modelList = ml)

##Predicted 5 Species and 

# predict data
#par(mfrow=c(1,3),bty='n')
#gjamPredict(out, y2plot = colnames(f$ydata)) #predict the data in-sample
#title('full sample')

# out-of-sample prediction
xdata     <- f$xdata[1:20,]
xdata[,3] <- mean(f$xdata[,3]) # mean for x[,3]
xdata[,4] <- mean(f$xdata[,4])# mean for x[,4]
xdata[,2] <- seq(-2.5,2.5,length=20)   # gradient x[,2]
newdata   <- list(xdata = xdata, nsim = 50 )
p1 <- gjamPredict(out, newdata = newdata)
# plus/minus 1 prediction SE, default effort = 1000
x2   <- p1$x[,2]
ylim <- c(0, max(p1$sdList$yMu[,1] + p1$sdList$yPe[,1]))

par(mfrow=c(2,3))
plot(x2, p1$sdList$yMu[,1],type='l',lwd=2, ylim=c(0, max(p1$sdList$yMu[,1] + p1$sdList$yPe[,1])), xlab='x2',
     ylab = 'Predicted')
mtext(paste(expression(beta),"=",round(out$parameters$betaMu[2,1],2)), side = 3)
lines(x2, p1$sdList$yMu[,1] + p1$sdList$yPe[,1], lty=2)
lines(x2, p1$sdList$yMu[,1] - p1$sdList$yPe[,1], lty=2)
# .95 prediction error
title('SE and prediction, Sp 1')

plot(x2, p1$sdList$yMu[,2],type='l',lwd=2, ylim=c(0, max(p1$sdList$yMu[,2] + p1$sdList$yPe[,2])), xlab='x2',
      ylab = 'Predicted')
mtext(paste(expression(beta),"=",round(out$parameters$betaMu[2,2],2)), side = 3)
lines(x2, p1$sdList$yMu[,2] + p1$sdList$yPe[,2], lty=2)
lines(x2, p1$sdList$yMu[,2] - p1$sdList$yPe[,2], lty=2)
title('SE and prediction, Sp 2')

plot(x2, p1$sdList$yMu[,3],type='l',lwd=2, ylim=c(0, max(p1$sdList$yMu[,3] + p1$sdList$yPe[,3])), xlab='x2',
      ylab = 'Predicted')
mtext(paste(expression(beta),"=",round(out$parameters$betaMu[2,3],3)), side = 3)
lines(x2, p1$sdList$yMu[,3] + p1$sdList$yPe[,3], lty=2)
lines(x2, p1$sdList$yMu[,3] - p1$sdList$yPe[,3], lty=2)
title('SE and prediction, Sp 3')

plot(x2, p1$sdList$yMu[,4],type='l',lwd=2, ylim=c(0, max(p1$sdList$yMu[,4] + p1$sdList$yPe[,4])), xlab='x2',
      ylab = 'Predicted')
mtext(paste(expression(beta),"=",round(out$parameters$betaMu[2,4],2)), side = 3)
lines(x2, p1$sdList$yMu[,4] + p1$sdList$yPe[,4], lty=2)
lines(x2, p1$sdList$yMu[,4] - p1$sdList$yPe[,4], lty=2)
title('SE and prediction, Sp 4')
plot(x2, p1$sdList$yMu[,5],type='l',lwd=2, ylim=c(0, max(p1$sdList$yMu[,5] + p1$sdList$yPe[,5])), xlab='x2',
      ylab = 'Predicted')
mtext(paste(expression(beta),"=",round(out$parameters$betaMu[2,5],5)), side = 3)
lines(x2, p1$sdList$yMu[,5] + p1$sdList$yPe[,5], lty=2)
lines(x2, p1$sdList$yMu[,5] - p1$sdList$yPe[,5], lty=2)

#lines(x2, p1$sdList$yMu[,1] + p1$sdList$yPe[,1], lty=2)
#lines(x2, p1$sdList$yMu[,1] - p1$sdList$yPe[,1], lty=2)
# .95 prediction error
#lines(x2, p1$piList$yLo[,1], lty=3)
#lines(x2, p1$piList$yHi[,1], lty=3)
title('SE and prediction, Sp 5')

#ind<-which((p1$sdList$yMu[,1]>0.5))
#densityplot(xdata[ind,2])

########## Partial Dependence plot 
n.sample<- 200
#xseq <- seq(-2,2,length=20)   # gradient x[,2]

presence_array<- matrix(nrow=n.sample, ncol =5)
for (i in 1:n.sample) {
  xtrain    <- f$xdata[1:200,]
  #xtrain[,2] <- xseq[i]    # value for ith predictor duplicate for all set 
  xtrain[,2] <-xtrain[i,2] 
  newdata   <- list(xdata = xtrain, nsim = 50)
  p_tr <- gjamPredict(out, newdata = newdata)
  presence_array[i,]<- colSums(p_tr$prPresent)/n.sample
}
par(mfrow=c(2,3))
plot(f$xdata[1:200,2],presence_array[,1],xlab=expression(x[2]),main="",ylab="Prediction of occurence S1",lwd=2)
rug(f$xdata[1:200,2], col="red")
plot(f$xdata[1:200,2],presence_array[,2],xlab=expression(x[2]),main="",ylab="Prediction of occurence S2",lwd=2)
rug(f$xdata[1:200,2], col="red")
plot(f$xdata[1:200,2],presence_array[,3],xlab=expression(x[2]),main="",ylab="Prediction of occurence S3",lwd=2)
rug(f$xdata[1:200,2], col="red")
plot(f$xdata[1:200,2],presence_array[,4],xlab=expression(x[2]),main="",ylab="Prediction of occurence S4",lwd=2)
rug(f$xdata[1:200,2], col="red")
plot(f$xdata[1:200,2],presence_array[,5],xlab=expression(x[2]),main="",ylab="Prediction of occurence S5",lwd=2)
rug(f$xdata[1:200,2], col="red")




n.sample<- 20
#xseq <- seq(-2,2,length=20)   # gradient x[,2]
presence_array<- matrix(nrow=n.sample, ncol =5)
pr_array<-array(0,dim=c(n.sample,n.sample,5))
for (i in 1:n.sample) {
  for (j in 1:n.sample){
    xtrain    <- f$xdata[1:20,]
    xtrain[,2] <-xtrain[i,2] 
    xtrain[,3] <-xtrain[j,3] 
    newdata   <- list(xdata = xtrain, nsim = 50)
    p_tr <- gjamPredict(out, newdata = newdata)
    pr_array[i,j,]<- colSums(p_tr$prPresent)/n.sample
  }
}
#library(plotly)



df <- data.frame(matrix(ncol = 3, nrow = 0))
x <- c("X2", "X3", "Prediction")
colnames(df) <- x

df<-data.frame()
for (i in 1:n.sample){
  for (j in 1:n.sample){
    df<- rbind(df, c(f$xdata[i,2],f$xdata[j,3],pr_array[i,j,1]))
  }
}

new_df <- df[order(df$X0.370312380933042, df$X.1.60280849294246),] 
K<-as.matrix(df)
p <- plot_ly(z=~K) %>% add_markers()
p

p <- plot_ly(x=df$X0.370312380933042,y=df$X.1.60280849294246,z=df$X0.737)
p




par(mfrow=c(2,3))
plot(f$xdata[1:20,2],presence_array[,1],xlab=expression(x[2]),main="",ylab="Prediction of occurence S1",lwd=2)
rug(f$xdata[1:20,2], col="red")
plot(f$xdata[1:20,2],presence_array[,2],xlab=expression(x[2]),main="",ylab="Prediction of occurence S2",lwd=2)
rug(f$xdata[1:20,2], col="red")
plot(f$xdata[1:20,2],presence_array[,3],xlab=expression(x[2]),main="",ylab="Prediction of occurence S3",lwd=2)
rug(f$xdata[1:20,2], col="red")
plot(f$xdata[1:20,2],presence_array[,4],xlab=expression(x[2]),main="",ylab="Prediction of occurence S4",lwd=2)
rug(f$xdata[1:20,2], col="red")
plot(f$xdata[1:20,2],presence_array[,5],xlab=expression(x[2]),main="",ylab="Prediction of occurence S5",lwd=2)
rug(f$xdata[1:20,2], col="red")




# conditional prediction
ydataCond <- out$inputs$y[,1,drop=FALSE]*0
newdata   <- list(ydataCond = ydataCond, nsim=50)
p0        <- gjamPredict(output = out, newdata = newdata)
ydataCond <- ydataCond + 20                  #first column is 20
newdata   <- list(ydataCond = ydataCond, nsim=50)
p1        <- gjamPredict(output = out, newdata = newdata)
plot(out$inputs$y[,4],p0$sdList$yMu[,4], cex=.4,col='orange'); abline(0,1,lty=2)

points(out$inputs$y[,4],p1$sdList$yMu[,4], cex=.4,col='blue')
title('Cond. on 1st Sp')







## Example:
S <-1
f <- gjamSimData(n = 200, S = S, Q = 3, typeNames = 'PA')
ml <- list(ng = 50, burnin = 5, typeNames = f$typeNames, holdoutN = 10) 
out <- gjam(f$formula, f$xdata, f$ydata, modelList = ml)

##Predicted 5 Species and 

# predict data
#par(mfrow=c(1,3),bty='n')
gjamPredict(out, y2plot = colnames(f$ydata)) #predict the data in-sample
title('full sample')

# out-of-sample prediction
xdata     <- f$xdata[1:20,]
xdata[,3] <- mean(f$xdata[,3])     # mean for x[,3]
xdata[,2] <- seq(-2,2,length=20)   # gradient x[,2]
newdata   <- list(xdata = xdata, nsim = 50 )
p1 <- gjamPredict(out, newdata = newdata)
# plus/minus 1 prediction SE, default effort = 1000
x2   <- p1$x[,2]
ylim <- c(0, max(p1$sdList$yMu[,1] + p1$sdList$yPe[,1]))
plot(x2, p1$sdList$yMu[,1],type='l',lwd=2, ylim=ylim, xlab='x2',
     ylab = 'Predicted')
lines(x2, p1$sdList$yMu[,2],type='l',lwd=3, ylim=ylim, xlab='x2',
      ylab = 'Predicted')
lines(x2, p1$sdList$yMu[,3],type='l',lwd=3, ylim=ylim, xlab='x2',
      ylab = 'Predicted')
lines(x2, p1$sdList$yMu[,4],type='l',lwd=3, ylim=ylim, xlab='x2',
      ylab = 'Predicted')
lines(x2, p1$sdList$yMu[,5],type='l',lwd=3, ylim=ylim, xlab='x2',
      ylab = 'Predicted')
#lines(x2, p1$sdList$yMu[,1] + p1$sdList$yPe[,1], lty=2)
#lines(x2, p1$sdList$yMu[,1] - p1$sdList$yPe[,1], lty=2)
# .95 prediction error
lines(x2, p1$piList$yLo[,1], lty=3)
lines(x2, p1$piList$yHi[,1], lty=3)
title('SE and prediction, Sp 1')

ind<-which((p1$sdList$yMu[,1]>0.5))
densityplot(xdata[ind,2])

hist(xdata[ind,2])

ind<-which((p1$sdList$yMu[,1]>0.5))
hist(xdata[ind,2])
densityplot(xdata[ind,2])



hist(xdata[,2])
densityplot(xdata[,2])


x<- seq(-2,2,length=20) 
mu<-2*x
sigma<-1

p<-rnorm(1000,mu,sigma)


#########Using Probit regression
# Partial Dependence plot 




#plot(out$prediction$presence[,1],out$prediction$xpredMu[,2])


# conditional prediction
ydataCond <- out$inputs$y[,1,drop=FALSE]*0
newdata   <- list(ydataCond = ydataCond, nsim=50)
p0        <- gjamPredict(output = out, newdata = newdata)
ydataCond <- ydataCond + 20                  #first column is 20
newdata   <- list(ydataCond = ydataCond, nsim=50)
p1        <- gjamPredict(output = out, newdata = newdata)
plot(out$inputs$y[,4],p0$sdList$yMu[,4], cex=.4,col='orange'); abline(0,1,lty=2)

points(out$inputs$y[,4],p1$sdList$yMu[,4], cex=.4,col='blue')
title('Cond. on 1st Sp')



# conditional prediction
ydataCond <- out$inputs$y[,1,drop=FALSE]*0
newdata   <- list(ydataCond = ydataCond, xdata = xdata, nsim=50)
p0        <- gjamPredict(output = out, newdata = newdata)
ydataCond <- ydataCond + 20                  #first column is 20
newdata   <- list(ydataCond = ydataCond, nsim=50)
p1        <- gjamPredict(output = out, newdata = newdata)
plot(out$inputs$y[,4],p0$sdList$yMu[,4], cex=.4,col='orange'); abline(0,1,lty=2)
points(out$inputs$y[,4],p1$sdList$yMu[,4], cex=.4,col='blue')
title('Cond. on 1st Sp')



summary(f$xdata)
q1<-qplot(f$xdata[,2], geom="histogram", xlab="x1") 
q2<-qplot(f$xdata[,3], geom="histogram", xlab="x2") 
q3<-qplot(f$xdata[,4], geom="histogram", xlab="x3") 
grid.arrange(q1, q2,q3, nrow = 1)


b1<-qplot(f$trueValues$beta[1,], geom="histogram", xlab="intercept") 
b2<-qplot(f$trueValues$beta[2,], geom="histogram", xlab="x2") 
b3<-qplot(f$trueValues$beta[3,], geom="histogram", xlab="x3") 
b4<-qplot(f$trueValues$beta[4,], geom="histogram", xlab="x4") 
grid.arrange(b1,b2,b3,b4, nrow = 1)


library(plotly)
p_S <- plot_ly(z = f$trueValues$sigma, colors = "Greys", type = "heatmap")
p_S

library(plotly)
p_R <- plot_ly(z = f$trueValues$corSpec, colors = "Greys", type = "heatmap")
p_R


par(bty = 'n', mfrow = c(1,2), family='')
h <- hist(c(-1,f$y),nclass = 50,plot = F)
plot(h$counts,h$mids,type = 's')
plot(f$w,f$y,cex = .2)

qplot(f$y, geom="histogram", xlab="y") 

ml  <- list(ng = 1000, burnin = 100, typeNames = f$typeNames)
out <- gjam(f$formula, f$xdata, f$ydata, modelList = ml)
summary(out)
















f   <- gjamSimData(n = 500, S = 10, typeNames = 'CA')
ml  <- list(ng = 1000, burnin = 200, typeNames = f$typeNames)
out <- gjam(f$formula, f$xdata, f$ydata, modelList = ml)
pl  <- list(trueValues = f$trueValues, GRIDPLOTS = T)
gjamPlot(output = out, plotPars = pl)


b<-out$parameters$betaMu  
summary(out$chains)








summary(out$chains)


out$chains$sgibbs

summary(out$modelSummary)

sim <- gjamSimData(n = 500, S = 10, q = 3, typeNames = 'CA')
ml <- list(ng = 2000, burnin = 500, typeNames = sim$typeNames)
out <- gjamGibbs(sim$formula, sim$xdata, sim$ydata, modelList = ml) 
pl <- list(trueValues = sim$trueValues,width = 3,height = 2,GRIDPLOTS = T, SMALLPLOTS = F) 
gjamPlot(output = out, plotPars = pl)

par(bty = 'n', mfrow = c(1,3)) 
plot(sim$trueValues$beta, out$modelSummary$betaMu) 
plot(sim$trueValues$corSpec, out$modelSummary$corMu) 
plot(sim$y,out$modelSummary$yMu, cex = .2)



#missing data

sim <- gjamSimData(S = 5, Q = 3, typeNames = 'CA', nmiss = 20)
ml <- list(ng = 2000, burnin = 500, typeNames = sim$typeNames, holdoutN = 50) 
out <- gjam(sim$formula, sim$xdata, sim$ydata, modelList = ml)
par(mfrow=c(1,3))
plot(out$x[out$missingIndex], out$modelSummary$xpredMu[out$missingIndex]) 
title('missing in x'); abline(0,1)
plot(out$x[out$holdoutIndex,-1], out$modelSummary$xpredMu[out$holdoutIndex,-1]) 
title('holdouts in x');abline(0,1)
