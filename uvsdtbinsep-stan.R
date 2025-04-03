uvsdtbinsep_stanvars <- "
   real uvsdtbinsep_lpmf(int y, real mu, real discsignal,
                   real cr,
                   int N, int type) {
    real disc = 1/discsignal;
    
     real p;

    // calculate probabilities
    if (type == 0) {
      p = 1 - Phi(disc * (cr - mu));
    } else if (type == 1) {
      p = 1 - Phi(cr);
    }
    
    return binomial_lpmf(y | N, p);
   }
"

uvsdtbinsep_family <- custom_family(
  name = "uvsdtbinsep", 
  dpars = c("mu", "discsignal", "cr"), 
  links = c("identity", "log", "identity"), lb = c(NA, 0, NA),
  type = "int", vars = c("vint1[n]", "vint2[n]")
)
sv_uvsdtbinsep <- stanvar(scode = uvsdtbinsep_stanvars, block = "functions")

calc_posterior_predictions_uvsdtbinsep <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  discsignal <- 1/brms::get_dpar(prep, "discsignal", i = i)
  cr <- brms::get_dpar(prep, "cr", i = i)
  type <- prep$data$vint2[i]
  
  # calculate probabilities
  if (type == 0) {
    p = 1 - pnorm(discsignal * (cr - mu))
  } else {
    p = 1 - pnorm(cr, 0)
  }
  return(p)
}

log_lik_uvsdtbinsep <- function(i, prep) {
  use <- calc_posterior_predictions_uvsdtbinsep(i = i, prep = prep)
  y <- prep$data$Y[i]
  N <- prep$data$vint1[i]
  dbinom(y, size = N, prob = use, log = TRUE) 
}

posterior_epred_uvsdtbinsep <- function(prep) {
  nobs <- prep$nobs
  out <- matrix(NA_real_, nrow = prep$ndraws, ncol = prep$nobs, 
               dimnames = list(seq(prep$ndraws), seq(prep$nobs)))
  for (i in seq_len(nobs)) {
    tmp <- calc_posterior_predictions_uvsdtbinsep(i = i, prep = prep)
    out[,i] <- tmp
  }
  return(out)
}

posterior_predict_uvsdtbinsep <- function(i, prep, ...) {
  use <- calc_posterior_predictions_uvsdtbinsep(i = i, prep = prep)
  y <- prep$data$Y[i]
  N <- prep$data$vint1[i]
  
  lout <- length(use)
  out <- rbinom(n = rep(1, lout), size = N, prob = use)
  out
  #apply(out, 1, function(x) list(x))
  #out[,1]
}
