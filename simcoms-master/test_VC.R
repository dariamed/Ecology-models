

simulate_community <- function(
  env = runif(500, 0, 100), niche_optima  = seq(2, 98, 5), niche_breadth = 20,
  type = "original", comp_inter = NA, fac_inter = NA, beta_env = 1,
  beta_comp = 5, beta_fac = 0, beta_abun = 0, years = 20, K = 40,
  competition = "facilitation", intra_sp_com  = 0
) {
  sim_com <- function(
    env, niche_breadth, niche_optima, type, comp_inter, fac_inter, beta_env,
    beta_comp, beta_fac, beta_abun, years, K, competition,intra_sp_com
  ) {
    n_sp = length(niche_optima)
    
    if (type == "original") {
      species_niche_overlap_sym <- outer(
        niche_optima,
        niche_optima,
        function(x, y) 2 * pnorm(-abs((x - y)) / 2, sd = niche_breadth)
      )
      diag(species_niche_overlap_sym) <- intra_sp_com
      species_fac_sym <- species_niche_overlap_sym
    } else {
      if (length(comp_inter) == 1) comp_inter = matrix(comp_inter, n_sp, n_sp)
      if (length(fac_inter)  == 1) fac_inter  = matrix(fac_inter, n_sp, n_sp)
      species_niche_overlap_sym <- comp_inter
      species_fac_sym <- fac_inter
    }
    
    species_niche_overlap_asym <- outer(
      niche_optima,
      niche_optima,
      function(x, y) {
        sign <- ifelse(x > y, 1, 0)
        overlap <- 2 * pnorm(-abs((x - y)) / 2, sd = niche_breadth)
        sign * overlap
      }
    )
    
    diag(species_niche_overlap_asym) <- intra_sp_com
    
    log_p_env <- sapply(
      niche_optima, dnorm, mean = env, sd = niche_breadth, log = TRUE
    )
    log_p_env <- log_p_env  - log(dnorm(0) / 10)
    
    community <- factor(
      x      = sample(seq_along(niche_optima), K, replace = TRUE),
      levels = seq_len(n_sp)
    )
    
    abund <- table(community)
    
    for (j in seq_len(years)) {
      for (k in seq_len(K)) {
        f_comp <- 1 - colSums(species_fac_sym[community, ]) / K
        p_comp <- 1 - colSums(species_niche_overlap_sym[community, ]) / K
        
        if (competition == "asymmetric") {
          p_comp <- 1 - colSums(species_niche_overlap_asym[community, ]) / K
        }
        
        if (competition == "facilitation") {
          p_all <- exp(
            beta_env * log_p_env - beta_fac * log(f_comp) + 
              log(1 + beta_abun * abund)
          )
        } else {
          p_all <- exp(
            beta_env * log_p_env + beta_comp * log(p_comp) +
              log(1 + beta_abun * abund)
          )
        }
        
        if (competition == "both") {
          p_all <- exp(
            beta_env * log_p_env + beta_comp * log(p_comp) - beta_fac *
              log(f_comp) + log(1 + beta_abun * abund)
          )
        }
        
        p_all <- ifelse(is.na(p_all), min(p_all, na.rm = TRUE), p_all)
        if (all(is.na(p_all)) || identical(min(p_all), max(p_all))) p_all = NULL
        if (any(is.infinite(p_all))) {
          community[sample(K, 1)] <- sample(seq_len(n_sp)[p_all == Inf], 1)
        } else {
          community[sample(K, 1)] <- sample(n_sp, 1, prob = p_all)
        }
        abund <- table(community)
      }
    }
    as.integer(abund) > 0
  }
  
  ans <- lapply(
    env, sim_com, niche_breadth, niche_optima, type, comp_inter, fac_inter,
    beta_env, beta_comp, beta_fac, beta_abun, years, K, competition,
    intra_sp_com
    #,mc.cores = detectCores()
  )
  ans <- do.call(rbind, ans)
  ans <- cbind(ans, env)
  sp_labs <- paste0(
    "sp_", gsub(" ", 0, format(seq_along(niche_optima), width = 2))
  )
  colnames(ans) <- c(sp_labs, "env")
  as.data.frame(ans)
}


test.coms <- simulate_community()


xdata<-as.data.frame(cbind(rep(1,500),test.coms$env))

colnames(xdata)<- c("int","env")
ydata<-as.data.frame(test.coms[1:20])
formula<-as.formula( ~env)

rl   <- list(r = 8, N = 20)
ml   <- list(ng = 2500, burnin = 500, typeNames = 'PA', reductList = rl)

out3  <- gjam(formula, xdata = xdata, ydata = ydata, modelList = ml)
summary(out3)


specNames <- colnames(ydata)
specColor <- rep('black',ncol(ydata))
specColor[ c(grep('quer',specNames),grep('cary',specNames)) ] <- 'brown'
specColor[ c(grep('acer',specNames),grep('frax',specNames)) ] <- 'darkgreen'
specColor[ c(grep('abie',specNames),grep('pice',specNames)) ] <- 'blue'

pl   <- list(GRIDPLOTS=T, specColor = specColor)
x<-gjamPlot(output = out3, plotPars = pl)





#x <- model.matrix(formula, xdata)
#qr(x)$rank
#dim(x)

#f<- gjamSimData(n=500, S=20,Q=2, typeNames = "PA")
#ml  <- list(ng = 1000, burnin = 100, typeNames = f$typeNames)
#out2 <- gjam(f$formula, f$xdata, f$ydata, modelList = ml)
#summary(out2)

Y_cor<-cor(ydata) 
R <- heatmap(out3$parameters$corMu, Rowv=NA, Colv=NA)
E <- heatmap(out3$parameters$ematrix, Rowv=NA, Colv=NA)
Yc <- heatmap(Y_cor, Rowv=NA, Colv=NA, col = cm.colors(256), scale="column", margins=c(5,10))


Prec<-100*solve(out3$parameters$sigMu)

Pr_mat <- heatmap(Prec,, Rowv=NA, Colv=NA)

P<- as.matrix(Pr_mat)


