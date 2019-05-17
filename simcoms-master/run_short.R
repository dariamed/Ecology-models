library(arm)
library(jagsUI)
library(ggplot2)
library(gridExtra)
library(parallel)
library(Rcpp)


setwd("~/Documents/GitHub/Ecology-models/simcoms-master/ExampleFiles")
load("params.rds")
load("sim_names.rds")
load("comp_inter.rds")
load("fac_inter.rds")

setwd("~/Documents/GitHub/Ecology-models/simcoms-master")
sim_data<-readRDS("sim_data.rds")

###############################################################################################


run_model <- function(data) {
  jsdm_jags <- function() {
    model.file <- tempfile()
    cat(
      "model {
      for (i in 1:n) {
      Z[i, 1:J] ~ dmnorm(Mu[i, ], Tau)
      for (j in 1:J) {
      Mu[i, j] <- inprod(B_raw[j, ], X[i, ])
      Y[i, j] ~ dbern(step(Z[i, j]))
      }
      }
      for (j in 1:J) {
      sigma_[j] <- sqrt(Sigma[j, j])
      env_sigma_[j] <- sqrt(EnvSigma[j, j])
      for (k in 1:K) {
      B_raw[j, k] ~ dnorm(mu[k], tau[k])
      B[j, k] <- B_raw[j, k] / sigma_[j]
      }
      for (j_ in 1:J) {
      Rho[j, j_] <- Sigma[j, j_] / (sigma_[j] * sigma_[j_])
      EnvRho[j, j_] <- EnvSigma[j, j_] / (env_sigma_[j] * env_sigma_[j_])
      EnvSigma[j, j_] <- sum(EnvSigma1[, j, j_]) + sum(EnvSigma2[, , j, j_])
      for (k in 2:K) {
      EnvSigma1[k - 1, j, j_] <- B[j, k] * B[j_, k]
      for (k_ in 2:K) {
      EnvSigma2[k - 1, k_ - 1, j, j_] <-
      B[j, k] * B[j_, k_] * ifelse(k_ != k, covx[k, k_], 0)
      }
      }
      }
      }
      for (k in 1:K) {
      mu[k] ~ dnorm(0, 1)
      tau[k] <- pow(sigma[k], -2)
      sigma[k] ~ dnorm(0, 1)T(0,)
      }
      Tau ~ dwish(I, df)
      Sigma <- inverse(Tau)
  }",
      file = model.file
    )
    model.file
    }
  
  inits <- function(data) {
    Y <- as.matrix(data$Y)
    X <- as.matrix(data$X)[, -1]
    Tau <- rWishart(1, data$df, data$I)[, , 1]
    Sigma <- solve(Tau)
    Z <- rep(0, data$J)
    Z <- mvrnorm(1, Z, Sigma)
    Z <- replicate(data$n, Z)
    Z <- t(Z)
    Z <- abs(Z)
    Z <- ifelse(Y, Z, -Z)
    Sigma <- cov(Z)
    B <- sapply(
      seq_len(data$J),
      function(x) coef(bayesglm(Y[, x] ~ X, family = binomial(link = "probit")))
    )
    B <- t(B)
    B_raw <- B * sqrt(diag(Sigma))
    mu <- apply(B_raw, 2, mean)
    sigma <- pmin(99, apply(B_raw, 2, sd))
    Tau <- solve(Sigma)
    list(Tau = Tau, Z = Z, B_raw = B_raw, mu = mu, sigma = sigma)
  }
  
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
  
  model <- jags(
    data,
    function() inits(data), c("B", "Rho", "EnvRho","Tau"), jsdm_jags(), ###added Tau
    n.chains = 5, n.iter = 2e3, n.adapt = 25e3, n.thin = 10,  ####changed chains from 15 to 5
    parallel = TRUE, DIC = FALSE
  )
  
  save(
    model,
    file = paste0("jsdmmodel-", format(Sys.time(), "%Y-%m-%d-%H-%M-%S"), ".rda")
  )
  model
  }


# Fit models

data<-sim_data$FacCompSparseSp5
M1<- run_model(data)
models <- lapply(sim_data, run_model)


