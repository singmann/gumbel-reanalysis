gumbel6agg_stanvars <- "
   real gumbelmin(real x, real mu, real disc){
     //return 1- exp(-exp(-(-x-mu)/disc));
     //return exp(gumbel_lccdf(-x | mu,disc));
     return 1 - gumbel_cdf(-x|mu,disc);
   }
   real gumbel6agg_lpmf(int y, real mu, 
                   real crc, real crlm, real crll, real crhm, real crhh, 
                   array[] int oldvec, array[] int newvec) {
    int nthres = 5;
    vector[nthres+1] pold;
    vector[nthres+1] pnew;
    
    real disc = 1;
    vector[nthres] thres;
    
    // calculate thresholds
    thres[1] = crc - (exp(crlm) + exp(crll));
    thres[2] = crc - (exp(crlm));
    thres[3] = crc;
    thres[4] = crc + (exp(crhm));
    thres[5] = crc + (exp(crhm) + exp(crhh));
     
    // calculate probabilities
    pold[1] = gumbelmin(thres[1], mu, disc);
    for (i in 2:nthres) {
      pold[i] = gumbelmin(thres[i], mu, disc) - gumbelmin(thres[i-1], mu, disc);
    }
    pold[6] = 1 - gumbelmin(thres[nthres], mu, disc);
    pnew[1] = gumbelmin(thres[1], 0, disc);
    for (i in 2:nthres) {
      pnew[i] = gumbelmin(thres[i], 0, disc) - gumbelmin(thres[i-1], 0, disc);
    }
    pnew[6] = 1 - gumbelmin(thres[nthres], 0, disc);
    return multinomial_lpmf(oldvec | pold) + multinomial_lpmf(newvec | pnew);
   }
"

gumbel6agg_family <- custom_family(
  name = "gumbel6agg", 
  dpars = c("mu", "crc", "crlm", "crll", "crhm", "crhh"), 
  links = c("identity", rep("identity", 5)), lb = c(NA, rep(NA, 5)),
  type = "int", vars = c("oldmat2[n]", "newmat2[n]")
)
sv_gumbel6agg <- stanvar(scode = gumbel6agg_stanvars, block = "functions")

calc_posterior_predictions_gumbel6agg <- function(i, prep) {
  gumbelmin <- function(x, mu, disc) {
    return(1 - extraDistr::pgumbel(-x,mu,disc))
  }
  mu <- brms::get_dpar(prep, "mu", i = i)
  #discsignal <- brms::get_dpar(prep, "discsignal", i = i)
  crc <- brms::get_dpar(prep, "crc", i = i)
  crlm <- brms::get_dpar(prep, "crlm", i = i)
  crll <- brms::get_dpar(prep, "crll", i = i)
  crhm <- brms::get_dpar(prep, "crhm", i = i)
  crhh <- brms::get_dpar(prep, "crhh", i = i)

  oldvec <- prep$data$oldmat[i,]
  newvec <- prep$data$newmat[i,]
  OUTLEN <- length(mu)
  
  disc <- 1
  nthres <- 5
  thres <- matrix(NA_real_, nrow = OUTLEN, ncol = nthres)
  pold <- matrix(NA_real_, nrow = OUTLEN, ncol = nthres+1)
  pnew <- matrix(NA_real_, nrow = OUTLEN, ncol = nthres+1)
  
  thres[,1] = crc - (exp(crlm) + exp(crll));
  thres[,2] = crc - (exp(crlm));
  thres[,3] = crc;
  thres[,4] = crc + (exp(crhm));
  thres[,5] = crc + (exp(crhm) + exp(crhh));
  
  # calculate probabilities
  pold[,1] = gumbelmin(thres[,1], mu, disc)
  #ordinal::pgumbel(thres[,1], mu, disc, max = FALSE) ## does NOT 
  for (j in 2:nthres) {
    pold[,j] = gumbelmin(thres[,j], mu, disc) - gumbelmin(thres[,j-1], mu, disc);
  }
  pold[,6] = 1 - gumbelmin(thres[,nthres], mu, disc);

  pnew[,1] = gumbelmin(thres[,1], 0, disc);
  for (j in 2:nthres) {
    pnew[,j] = gumbelmin(thres[,j], 0, disc) - gumbelmin(thres[,j-1], 0, disc);
  }
  pnew[,6] = 1 - gumbelmin(thres[,nthres], 0, disc);
  return(list(
    pold = pold, pnew = pnew
  ))
}

log_lik_gumbel6agg <- function(i, prep) {
  use <- calc_posterior_predictions_gumbel6agg(i = i, prep = prep)
  oldvec <- prep$data$oldmat[i,]
  newvec <- prep$data$newmat[i,]
  extraDistr::dmnom(x = oldvec, size = sum(oldvec), prob = use$pold, log = TRUE) + 
     extraDistr::dmnom(x = newvec, size = sum(newvec), prob = use$pnew, log = TRUE)
}

posterior_epred_gumbel6agg <- function(prep) {
  nobs <- prep$nobs
  out <- array(NA_real_, dim = c(prep$ndraws, prep$nobs, 12), 
               dimnames = list(seq(prep$ndraws), seq(prep$nobs), 
                               c(colnames(prep$data$oldmat), colnames(prep$data$newmat))))
  for (i in seq_len(nobs)) {
    tmp <- calc_posterior_predictions_gumbel6agg(i = i, prep = prep)
    out[,i,] <- c(tmp$pold, tmp$pnew)
  }
  return(out)
}

posterior_predict_gumbel6agg <- function(i, prep, ...) {
  use <- calc_posterior_predictions_gumbel6agg(i = i, prep = prep)
  oldvec <- prep$data$oldmat[i,]
  newvec <- prep$data$newmat[i,]
  
  lout <- nrow(use$pold)
  out <- cbind(extraDistr::rmnom(n = rep(1, lout), size = sum(oldvec), prob = use$pold), 
               extraDistr::rmnom(n = rep(1, lout), size = sum(oldvec), prob = use$pnew))
  colnames(out) <- c(colnames(prep$data$oldmat), colnames(prep$data$newmat))
  #browser()
  lapply(seq_len(nrow(out)), function(i) out[i,])
  #apply(out, 1, function(x) list(x))
  #out[,1]
}

# log_lik_gumbel6agg <- function(i, prep) {
#   gumbelmin <- function(x, mu, disc) {
#     return(1 - extraDistr::pgumbel(-x,mu,disc))
#   }
#   mu <- brms::get_dpar(prep, "mu", i = i)
#   #discsignal <- brms::get_dpar(prep, "discsignal", i = i)
#   crc <- brms::get_dpar(prep, "crc", i = i)
#   crlm <- brms::get_dpar(prep, "crlm", i = i)
#   crll <- brms::get_dpar(prep, "crll", i = i)
#   crhm <- brms::get_dpar(prep, "crhm", i = i)
#   crhh <- brms::get_dpar(prep, "crhh", i = i)
# 
#   oldvec <- prep$data$oldmat[i,]
#   newvec <- prep$data$newmat[i,]
#   OUTLEN <- length(mu)
#   
#   disc <- 1
#   nthres <- 5
#   thres <- matrix(NA_real_, nrow = OUTLEN, ncol = nthres)
#   pold <- matrix(NA_real_, nrow = OUTLEN, ncol = nthres+1)
#   pnew <- matrix(NA_real_, nrow = OUTLEN, ncol = nthres+1)
#   
#   thres[,1] = crc - (exp(crlm) + exp(crll));
#   thres[,2] = crc - (exp(crlm));
#   thres[,3] = crc;
#   thres[,4] = crc + (exp(crhm));
#   thres[,5] = crc + (exp(crhm) + exp(crhh));
#   
#   # calculate probabilities
#   pold[,1] = gumbelmin(thres[,1], mu, disc)
#   #ordinal::pgumbel(thres[,1], mu, disc, max = FALSE) ## does NOT 
#   for (j in 2:nthres) {
#     pold[,j] = gumbelmin(thres[,j], mu, disc) - gumbelmin(thres[,j-1], mu, disc);
#   }
#   pold[,6] = 1 - gumbelmin(thres[,nthres], mu, disc);
# 
#   pnew[,1] = gumbelmin(thres[,1], 0, disc);
#   for (j in 2:nthres) {
#     pnew[,j] = gumbelmin(thres[,j], 0, disc) - gumbelmin(thres[,j-1], 0, disc);
#   }
#   pnew[,6] = 1 - gumbelmin(thres[,nthres], 0, disc);
#   
#   #browser()
#   extraDistr::dmnom(x = oldvec, size = sum(oldvec), prob = pold, log = TRUE) + 
#     extraDistr::dmnom(x = newvec, size = sum(newvec), prob = pnew, log = TRUE)
# }

# posterior_predict_gumbel6agg <- function(i, prep, ...) {
#   mu <- brms::get_dpar(prep, "mu", i = i)
#   #discsignal <- brms::get_dpar(prep, "discsignal", i = i)
#   crc <- brms::get_dpar(prep, "crc", i = i)
#   crlm <- brms::get_dpar(prep, "crlm", i = i)
#   crll <- brms::get_dpar(prep, "crll", i = i)
#   crhm <- brms::get_dpar(prep, "crhm", i = i)
#   crhh <- brms::get_dpar(prep, "crhh", i = i)
#   itemtype <- prep$data$vint1[i]
#   nsamples <- length(mu)
#   nthres <- 5
#   thres <- matrix(NA_real_, nrow = nsamples, ncol = nthres)
#   p <- matrix(NA_real_, nrow = nsamples, ncol = nthres+1)
#   disc = rep(1, nsamples)
#   thres[,1] = crc - (exp(crlm) + exp(crll));
#   thres[,2] = crc - (exp(crlm));
#   thres[,3] = crc;
#   thres[,4] = crc + (exp(crhm));
#   thres[,5] = crc + (exp(crhm) + exp(crhh));
#   
#   for (y in 1:(nthres+1)) {
#     if (y == 1) {
#       p[,y] = ordinal::pgumbel(thres[,1], mu, disc, max = FALSE)
#     } else if (y == nthres + 1) {
#       p[,y] = 1 - ordinal::pgumbel(thres[,nthres], mu, disc, max = FALSE)
#     } else {
#       p[,y] = ordinal::pgumbel(thres[,y], mu, disc, max = FALSE) -
#           ordinal::pgumbel(thres[,y-1], mu, disc, max = FALSE)
#     }
#   }
#   apply(p, 1, extraDistr::rcat, n = 1)
# }
