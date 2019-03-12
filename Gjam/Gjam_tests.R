library(gjam)
g1<- gjamSimData(n=500,S=10,Q=4, typeNames = "CA")

summary(g1)

g1$formula

par(bty = 'n', mfrow = c(1,2), family='')
h <- hist(c(-1,g1$y),nclass = 50,plot = T)
plot(h$counts,h$mids,type = 's')
plot(g1$w,g1$y,cex = .2)

ml  <- list(ng = 1000, burnin = 100, typeNames = g1$typeNames)
out <- gjam(g1$formula, g1$xdata, g1$ydata, modelList = ml)
summary(out)




## Not run: 
## ordinal data, show true parameter values
sim <- gjamSimData(S = 10, typeNames = 'OC')  
sim$ydata[1:5,]                              # example data
sim$trueValues$cuts                          # simulated partition
sim$trueValues$beta                          # coefficient matrix


library(plotly)
p <- plot_ly(z = sim$trueValues$sigma, colors = "Greys", type = "heatmap")
p
#chart_link = api_create(p, filename="heatmap-simple")
#chart_link

sim <- gjamSimData(n = 5, S = 5, typeNames = 'CA')  
sim$w
sim$y
plot(sim$w,sim$y)

types <- c(rep('DA',5), rep('CA',4))
sim   <- gjamSimData(n = 10, S = length(types), Q = 4, typeNames = types)
sim$typeNames
sim$ydata

sim <- gjamSimData(n = 10, S = 8, typeNames = 'CC')
totalCount <- rowSums(sim$ydata)
cbind(sim$ydata, totalCount)  # data with sample effort



## multiple categorical responses - compare matrix y and data.frqme ydata
types <- rep('CAT',2)
sim   <- gjamSimData(S = length(types), typeNames = types)
head(sim$ydata)
head(sim$y)

## discrete abundance, heterogeneous effort 
S   <- 5
n   <- 1000
ef  <- list( columns = 1:S, values = round(runif(n,.5,5),1) )
sim <- gjamSimData(n, S, typeNames = 'DA', effort = ef)
sim$effort$values[1:20]

## combinations of scales, partition only for 'OC' columns
types <- c('OC','OC','OC','CC','CC','CC','CC','CC','CA','CA','PA','PA')
sim   <- gjamSimData(S = length(types), typeNames = types)
sim$typeNames                           
head(sim$ydata)
sim$trueValues$cuts



######Predictions
f   <- gjamSimData(n = 500, S = 10, typeNames = 'CA')
ml  <- list(ng = 1000, burnin = 200, typeNames = f$typeNames)
out <- gjam(f$formula, f$xdata, f$ydata, modelList = ml)
pl  <- list(trueValues = f$trueValues, GRIDPLOTS = T)
gjamPlot(output = out, plotPars = pl)


f   <- gjamSimData(n = 500, S = 10, typeNames = 'CA')
ml  <- list(ng = 1000, burnin = 200, typeNames = f$typeNames)
out <- gjam(f$formula, f$xdata, f$ydata, modelList = ml)
pl  <- list(trueValues = f$trueValues, GRIDPLOTS = T)
gjamPlot(output = out, plotPars = pl)


f$ydata



