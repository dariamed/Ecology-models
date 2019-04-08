library(R2jags)
# An example model file is given in:
model.file <- system.file(package = "R2jags", "model", "schools.txt") 
 # data
J <- 8.0
y <- c(28.4,7.9,-2.8,6.8,-0.6,0.6,18.0,12.2)
sd <- c(14.9,10.2,16.3,11.0,9.4,11.4,10.4,17.6)
jags.data <- list("y","sd","J")
jags.params <- c("mu","sigma","theta")
jags.inits <- function(){
  list("mu"=rnorm(1),"sigma"=runif(1),"theta"=rnorm(J))
}
# Fit the model
jagsfit <- jags(data=list("y","sd","J"), inits = jags.inits,
                jags.params, n.iter = 10, model.file = model.file)

n.sim<-100
set.seed(123)
x1<-rnorm(n.sim,mean=5,sd=2)
x2<-rbinom(n.sim, size=1, prob=0.3)
e<-rnorm(n.sim, mean=0, sd=1)

b1<- 1.2
b2<- -3.1
a<- 1.5
y<- a+ b1*x1 + b2*x2 +e

sim.dat<-data.frame(y,x1,x2)

freq.mod<- lm(y ~ x1 + x2, data= sim.dat)
summary(freq.mod)



bayes.mod<-function(){
  for(i in 1:N){
    y[i] ~ dnorm(mu[i], tau)
    mu[i]<-alpha + beta1 * x1[i]+ beta2 * x2[i] }
  alpha ~ dnorm(0, .01)
  beta1~ dunif(-100, 100)
  beta2~ dunif(-100, 100)
  tau ~dgamma(.01, .01)
 }

bayes.mod <- function() {
  for(i in 1:N){
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha + beta1 * x1[i] + beta2 * x2[i] }
  alpha ~ dnorm(0, .01)
  beta1 ~ dunif(-100, 100)
  beta2 ~ dunif(-100, 100)
  tau ~ dgamma(.01, .01)
  }


y<-sim.dat$y
x1<-sim.dat$x1
x2<-sim.dat$x2
N<-nrow(sim.dat)


sim.dat.jags<- list("y","x1","x2","N")
#sim.dat.jags<- as.list(sim.dat)
#sim.dat.jags$N<-nrow(sim.dat)

bayes.mod.params<- c("alpha", "beta1", "beta2")


bayes.mod.inits<- function(){
  list("alpha"=rnorm(1), "beta1"=rnorm(1), "beta2"=rnorm(1))
}

bayes.mod.fit <- jags(data = sim.dat.jags, inits = bayes.mod.inits,
                          parameters.to.save = bayes.mod.params, n.chains = 3, n.iter = 9000,
                          n.burnin = 1000, model.file = bayes.mod)

bayes.mod.fit.upd <- update(bayes.mod.fit, n.iter=1000)
bayes.mod.fit.upd <- autojags(bayes.mod.fit)

print(bayes.mod.fit)
plot(bayes.mod.fit)
traceplot(bayes.mod.fit)

#pdf("bayes_trace.pdf")
#traceplot(bayes.mod.fit)
#dev.off()


bayes.mod.fit.mcmc <- as.mcmc(bayes.mod.fit)

summary(bayes.mod.fit.mcmc)
#xyplot(bayes.mod.fit.mcmc)
#densityplot(bayes.mod.fit.mcmc)
#densityplot(bayes.mod.fit.mcmc, layout=c(2,2), aspect="fill")
pdf("bayes_fit_mcmc_autocorr.pdf")
autocorr.plot(bayes.mod.fit.mcmc)
dev.off()
