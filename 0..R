
calc_posterior_predictions_gumbel6agg <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  OUTLEN <- length(mu)
  
  p <- matrix(NA_real_, nrow = OUTLEN, ncol = 2)
  
  G1<-function(x, mu){
    ordinal::dgumbel(x, max = FALSE)*(1-ordinal::pgumbel(x,mu, max = FALSE))**2
  }


  for (j in seq_len(OUTLEN)) {
    p[j, 1] <- integrate(G1,-Inf,Inf, mu = -mu[j],
                         rel.tol = .Machine$double.eps^0.5)$value    
  }
  p[, 2] <- 1-p[,1]
  return(p)
}

posterior_epred_gumbel6agg <- function(prep) {
  nobs <- prep$nobs
  out <- array(NA_real_, dim = c(prep$ndraws, prep$nobs, 2), 
               dimnames = list(seq(prep$ndraws), seq(prep$nobs), 
                               c("correct", "incorrect")))
  for (i in seq_len(nobs)) {
    tmp <- calc_posterior_predictions_gumbel6agg(i = i, prep = prep)
    out[,i,] <- tmp
  }
  return(out)
}

#############
calc_posterior_predictions_uvsdt6agg <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  discsignal <- brms::get_dpar(prep, "discsignal", i = i)

  OUTLEN <- length(mu)
  
  p <- matrix(NA_real_, nrow = OUTLEN, ncol = 2)
  
  G1<-function(x, mu, sd){
    dnorm(x)*(1-pnorm(x,mu,sd))**2
  }

  for (j in seq_len(OUTLEN)) {
    p[j, 1] <- integrate(G1,-Inf,Inf, mu = mu[j], sd = discsignal[j],
                         rel.tol = .Machine$double.eps^0.5)$value    
  }
  p[, 2] <- 1-p[,1]
  return(p)
}

posterior_epred_uvsdt6agg <- function(prep) {
  nobs <- prep$nobs
  out <- array(NA_real_, dim = c(prep$ndraws, prep$nobs, 2), 
               dimnames = list(seq(prep$ndraws), seq(prep$nobs), 
                               c("correct", "incorrect")))
  for (i in seq_len(nobs)) {
    tmp <- calc_posterior_predictions_uvsdt6agg(i = i, prep = prep)
    out[,i,] <- tmp
  }
  return(out)
}
