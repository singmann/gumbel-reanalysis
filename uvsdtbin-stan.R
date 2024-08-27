uvsdtbin_stanvars <- "
   real uvsdtbin_lpmf(int y, real mu, real discsignal,
                   real cr,
                   int Nold, int ynew, int Nnew) {

    real pold;
    real pnew;

    // calculate probabilities
    pold = 1 - Phi(discsignal * (cr - mu));
    pnew = 1 - Phi(cr);
    
    return binomial_lpmf(y | Nold, pold) + binomial_lpmf(ynew | Nnew, pnew);
   }
"

uvsdtbin_family <- custom_family(
  name = "uvsdtbin", 
  dpars = c("mu", "discsignal", "cr"), 
  links = c("identity", "log", "identity"), lb = c(NA, 0, NA),
  type = "int", vars = c("vint1[n]", "vint2[n]", "vint3[n]")
)
sv_uvsdtbin <- stanvar(scode = uvsdtbin_stanvars, block = "functions")

calc_posterior_predictions_uvsdtbin <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  discsignal <- brms::get_dpar(prep, "discsignal", i = i)
  cr <- brms::get_dpar(prep, "cr", i = i)

  OUTLEN <- length(mu)
  pold <- vector("numeric", OUTLEN)
  pnew <- vector("numeric", OUTLEN)
  
  # calculate probabilities
  pold = 1 - pnorm(discsignal * (cr - mu))
  pnew = 1 - pnorm(cr, 0)

    return(list(
    pold = pold, pnew = pnew
  ))
}

log_lik_uvsdtbin <- function(i, prep) {
  use <- calc_posterior_predictions_uvsdtbin(i = i, prep = prep)
  yold <- prep$data$Y[i]
  Nold <- prep$data$vint1[i]
  ynew <- prep$data$vint2[i]
  Nnew <- prep$data$vint3[i]
  
  dbinom(yold, size = Nold, prob = use$pold, log = TRUE) +
    dbinom(ynew, size = Nnew, prob = use$pnew, log = TRUE)
}

posterior_epred_uvsdtbin <- function(prep) {
  nobs <- prep$nobs
  out <- array(NA_real_, dim = c(prep$ndraws, prep$nobs, 2), 
               dimnames = list(seq(prep$ndraws), seq(prep$nobs), 
                               c("old", "new")))
  for (i in seq_len(nobs)) {
    tmp <- calc_posterior_predictions_uvsdtbin(i = i, prep = prep)
    out[,i,] <- c(tmp$pold, tmp$pnew)
  }
  return(out)
}

posterior_predict_uvsdtbin <- function(i, prep, ...) {
  use <- calc_posterior_predictions_uvsdtbin(i = i, prep = prep)
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
