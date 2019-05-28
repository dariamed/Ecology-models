


require(mvtnorm)
# Library for sampling from Truncated Normal distribution
require(truncnorm)

uni_bayes_probit<-function(datab,D=3, i_col=1, N_sim=10000){
  
  datapr<- dataprep(datab)
  N<- nrow(datapr$X)
  X <- as.matrix(datapr$X)
  # Binary observation data y
  y <- datapr$Y[,i_col]
  
  
  theta_0 <- rep(0, D)
  Q_0 <- diag(10, D)
  
  # Initialize parameters
  theta <- rep(0, D)
  z <- rep(0, N)
  
  # Number of simulations for Gibbs sampler
  N_sim <- N_sim
  # Burn in period
  burn_in <- 5000
  # Matrix storing samples of the \theta parameter
  theta_chain <- matrix(0, nrow = N_sim, ncol = D)
  
  theta_chain <- matrix(0, nrow = N_sim, ncol = D)
  # ---------------------------------
  # Gibbs sampling algorithm
  # ---------------------------------
  
  # Compute posterior variance of theta
  prec_0 <- solve(Q_0)
  V <- solve(prec_0 + crossprod(X, X))
  
  N1  <- sum(y)  # Number of successes
  N0  <- N - N1 
  
  for (t in 2:N_sim) {
    # Update Mean of z
    mu_z <- X %*% theta
    # Draw latent variable z from its full conditional: z | \theta, y, X
    z[y == 0] <- rtruncnorm(N0, mean = mu_z[y == 0], sd = 1, a = -Inf, b = 0)
    z[y == 1] <- rtruncnorm(N1, mean = mu_z[y == 1], sd = 1, a = 0, b = Inf)
    
    # Compute posterior mean of theta
    M <- V %*% (prec_0 %*% theta_0 + crossprod(X, z))
    # Draw variable \theta from its full conditional: \theta | z, X
    theta <- c(rmvnorm(1, M, V))
    
    # Store the \theta draws
    theta_chain[t, ] <- theta
  }
  
  # ---------------------------
  # Get posterior mean of \theta
  # ---------------------------
  post_theta_mean <- colMeans(theta_chain[-(1:burn_in), ])
  post_theta_05 <- apply(theta_chain[-(1:burn_in), ], 2, quantile, 0.05)
  post_theta_95 <- apply(theta_chain[-(1:burn_in), ], 2, quantile, 0.95)
  # Plot covariates x versus observations y
  # plot(X[,2], y, main = "Synthetic data")
  # # Show the fitted function using the posterior mean estimates
  # lines(x = X[,2], y = pnorm(X %*% post_theta_mean), col = "blue3", lwd = 2)
  # legend("bottomright", legend=c("MLE","Post. Mean"), col=c("red3","blue3"), 
  #        bty = 'n', lwd = 2, inset = c(0.02, 0.08), lty = 1, cex = 0.9)
  # predict_mean<- pnorm(X %*% post_theta)
  # predict_q05<- pnorm(X %*% post_theta)
  # predict_q95<- pnorm(X %*% post_theta)
  
  return(list(beta =theta_chain))
}
