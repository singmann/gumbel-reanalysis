### (int y, real mu, int r, int maxrank)
gumbelafctrial_stanvars <- "
  real gumbelafctrial_lpmf(int y, real mu, int m) {
    real g = mu;
    real r = 1;
    
    real log_p;
    real e_neg_g = exp(-g);
    // p = exp(-g) * beta(exp(-g), m);
    log_p = -g + lbeta(e_neg_g, m);
    return bernoulli_logit_lpmf(y | log_p - log1m_exp(log_p));
    
    //real log_p;
    //real e_neg_g = exp(-g);
    //log_p = (-g + lgamma(m) + lgamma(r - 1 + e_neg_g) - lgamma(r) - lgamma(m + e_neg_g));
    //return bernoulli_logit_lpmf(y | log_p - log1m_exp(log_p));
    
    //real p;
    //p = (exp(-g)*tgamma(m)*tgamma(r-1 + exp(-g))) / (tgamma(r)*tgamma(m + exp(-g)));
    //return bernoulli_lpmf(y | p);
  }
"

gumbelafctrial_family <- custom_family(
  name = "gumbelafctrial", 
  dpars = c("mu"), 
  links = c("identity"), lb = c(NA),
  type = "int", vars = paste0("vint", 1, "[n]")
)
sv_gumbelafctrial <- stanvar(scode = gumbelafctrial_stanvars, block = "functions") 

calc_posterior_predictions_gumbelafctrial <- function(i, prep) {
  g <- brms::get_dpar(prep, "mu", i = i)
  OUTLEN <- length(g)
  m <- prep$data$vint1[i]
  r <- 1
  
  e_neg_g = exp(-g);
  log_p = (-g + lgamma(m) + lgamma(r - 1 + e_neg_g) - lgamma(r) - lgamma(m + e_neg_g));
  
  #p = (exp(-g)*gamma(m)*gamma(r-1 + exp(-g))) / (gamma(r)*gamma(m + exp(-g)));
  return(exp(log_p))
}

log_lik_gumbelafctrial <- function(i, prep) {
  p <- calc_posterior_predictions_gumbelafctrial(i = i, prep = prep)
  out <- p
  out[prep$data$Y[i] == 0] <- 1 - p
  log(out)
}

posterior_epred_gumbelafctrial <- function(prep) {
  nobs <- prep$nobs
  out <- matrix(NA_real_, nrow = prep$ndraws, ncol = prep$nobs)
  for (i in seq_len(nobs)) {
    tmp <- calc_posterior_predictions_gumbelafctrial(i = i, prep = prep)
    out[,i] <- tmp
  }
  return(out)
}

posterior_predict_gumbelafctrial <- function(i, prep, ...) {
  p <- calc_posterior_predictions_gumbelafctrial(i = i, prep = prep)
  rbinom(length(p), 1, p)
}

