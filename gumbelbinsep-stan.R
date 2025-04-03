gumbelbinsep_stanvars <- "
   real gumbelmin(real x, real mu, real disc){
     //return 1- exp(-exp(-(-x-mu)/disc));
     //return exp(gumbel_lccdf(-x | mu,disc));
     return 1 - gumbel_cdf(-x|mu,disc);
   }
   real gumbelbinsep_lpmf(int y, real mu, 
                   real cr,
                   int N, int type) {

    real disc = 1;
    real p;

    // calculate probabilities
    if (type == 0) {
      p = 1 - gumbelmin(cr, -mu, disc);
    } else if (type == 1) {
      p = 1 - gumbelmin(cr, 0, disc);
    }
    
    return binomial_lpmf(y | N, p);
   }
"

gumbelbinsep_family <- custom_family(
  name = "gumbelbinsep", 
  dpars = c("mu", "cr"), 
  links = c("identity", "identity"), 
  type = "int", vars = c("vint1[n]", "vint2[n]")
)
sv_gumbelbinsep <- stanvar(scode = gumbelbinsep_stanvars, block = "functions")

calc_posterior_predictions_gumbelbinsep <- function(i, prep) {
  gumbelmin <- function(x, mu, disc) {
    return(1 - extraDistr::pgumbel(-x,mu,disc))
  }
  mu <- brms::get_dpar(prep, "mu", i = i)
  #discsignal <- brms::get_dpar(prep, "discsignal", i = i)
  cr <- brms::get_dpar(prep, "cr", i = i)
  type <- prep$data$vint2[i]

  OUTLEN <- length(mu)
  
  disc <- 1
  #p <- vector("numeric", OUTLEN)
  
  # calculate probabilities
  if (type == 0) {
    p <- 1 - gumbelmin(cr, -mu, disc)  
  } else {
    p <- 1 - gumbelmin(cr, 0, disc)  
  }
  return(p)
}

log_lik_gumbelbinsep <- function(i, prep) {
  use <- calc_posterior_predictions_gumbelbinsep(i = i, prep = prep)
  y <- prep$data$Y[i]
  N <- prep$data$vint1[i]
  
  dbinom(y, size = N, prob = use, log = TRUE) 
}

posterior_epred_gumbelbinsep <- function(prep) {
  nobs <- prep$nobs
  out <- matrix(NA_real_, nrow = prep$ndraws, ncol = prep$nobs, 
               dimnames = list(seq(prep$ndraws), seq(prep$nobs)))
  for (i in seq_len(nobs)) {
    tmp <- calc_posterior_predictions_gumbelbinsep(i = i, prep = prep)
    out[,i] <- tmp
  }
  return(out)
}

posterior_predict_gumbelbinsep <- function(i, prep, ...) {
  use <- calc_posterior_predictions_gumbelbinsep(i = i, prep = prep)
  y <- prep$data$Y[i]
  N <- prep$data$vint1[i]

  lout <- length(use)

  out <- rbinom(n = rep(1, lout), size = N, prob = use)
  out
  #browser()
  #lapply(seq_len(nrow(out)), function(i) out[i,])
  #apply(out, 1, function(x) list(x))
  #out[,1]
}
