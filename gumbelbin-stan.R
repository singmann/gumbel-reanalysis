gumbelbin_stanvars <- "
   real gumbelmin(real x, real mu, real disc){
     //return 1- exp(-exp(-(-x-mu)/disc));
     //return exp(gumbel_lccdf(-x | mu,disc));
     return 1 - gumbel_cdf(-x|mu,disc);
   }
   real gumbelbin_lpmf(int y, real mu, 
                   real cr,
                   int Nold, int ynew, int Nnew) {

    real disc = 1;
    real pold;
    real pnew;

    // calculate probabilities
    pold = 1 - gumbelmin(cr, -mu, disc);
    pnew = 1 - gumbelmin(cr, 0, disc);
    
    return binomial_lpmf(y | Nold, pold) + binomial_lpmf(ynew | Nnew, pnew);
   }
"

gumbelbin_family <- custom_family(
  name = "gumbelbin", 
  dpars = c("mu", "cr"), 
  links = c("identity", "identity"), 
  type = "int", vars = c("vint1[n]", "vint2[n]", "vint3[n]")
)
sv_gumbelbin <- stanvar(scode = gumbelbin_stanvars, block = "functions")

calc_posterior_predictions_gumbelbin <- function(i, prep) {
  gumbelmin <- function(x, mu, disc) {
    return(1 - extraDistr::pgumbel(-x,mu,disc))
  }
  mu <- brms::get_dpar(prep, "mu", i = i)
  #discsignal <- brms::get_dpar(prep, "discsignal", i = i)
  cr <- brms::get_dpar(prep, "cr", i = i)

  OUTLEN <- length(mu)
  
  disc <- 1
  pold <- vector("numeric", OUTLEN)
  pnew <- vector("numeric", OUTLEN)
  
  # calculate probabilities
  pold = 1 - gumbelmin(cr, -mu, disc)
  pnew = 1 - gumbelmin(cr, 0, disc)

    return(list(
    pold = pold, pnew = pnew
  ))
}

log_lik_gumbelbin <- function(i, prep) {
  use <- calc_posterior_predictions_gumbelbin(i = i, prep = prep)
  yold <- prep$data$Y[i]
  Nold <- prep$data$vint1[i]
  ynew <- prep$data$vint2[i]
  Nnew <- prep$data$vint3[i]
  
  dbinom(yold, size = Nold, prob = use$pold, log = TRUE) +
    dbinom(ynew, size = Nnew, prob = use$pnew, log = TRUE)
  }

posterior_epred_gumbelbin <- function(prep) {
  nobs <- prep$nobs
  out <- array(NA_real_, dim = c(prep$ndraws, prep$nobs, 2), 
               dimnames = list(seq(prep$ndraws), seq(prep$nobs), 
                               c("old", "new")))
  for (i in seq_len(nobs)) {
    tmp <- calc_posterior_predictions_gumbelbin(i = i, prep = prep)
    out[,i,] <- c(tmp$pold, tmp$pnew)
  }
  return(out)
}

posterior_predict_gumbelbin <- function(i, prep, ...) {
  use <- calc_posterior_predictions_gumbelbin(i = i, prep = prep)
  yold <- prep$data$Y[i]
  Nold <- prep$data$vint1[i]
  ynew <- prep$data$vint2[i]
  Nnew <- prep$data$vint3[i]
  
  lout <- nrow(use$pold)
  out <- cbind(old = rbinom(n = rep(1, lout), size = Nold, prob = use$pold), 
               new = rbinom(n = rep(1, lout), size = Nnew, prob = use$pnew))
  #browser()
  lapply(seq_len(nrow(out)), function(i) out[i,])
  #apply(out, 1, function(x) list(x))
  #out[,1]
}
