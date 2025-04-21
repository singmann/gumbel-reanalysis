source("gumbelmin_dist-stan.R")
## cat(gumbelmin_dist)
gumbel8agglog_stanvars <- "
   real gumbel8agglog_lpmf(int y, real mu, real crc, 
                    real crlm, real crll, real crlx, 
                    real crhm, real crhh, real crhx,
                    int y1, int y2, int y3, int y4, int y5, int y6, int y7, 
                    int y8, int y9, int y10, int y11, int y12, int y13, int y14, int y15) {
    int nthres = 7;
    vector[nthres+1] pold;
    vector[nthres+1] pnew;
    array[nthres+1] int oldvec = { y, y1, y2, y3, y4, y5, y6, y7 };
    array[nthres+1] int newvec = { y8, y9, y10, y11, y12, y13, y14, y15 };
    real tmp;
    
    real disc = 1;
    vector[nthres] thres;
    
    // calculate thresholds
    thres[1] = crc - (crlm + crll + crlx);
    thres[2] = crc - (crlm + crll);
    thres[3] = crc - (crlm);
    thres[4] = crc;
    thres[5] = crc + (crhm);
    thres[6] = crc + (crhm + crhh);
    thres[7] = crc + (crhm + crhh + crhx);
     
    // calculate probabilities
    pold[1] = gumbelmin_lcdf(thres[1]| mu);
    for (i in 2:nthres) {
      tmp = gumbelmin_lcdf(thres[i]| mu);
      pold[i] = tmp + log1m_exp(gumbelmin_lcdf(thres[i-1]| mu) - tmp);
    }
    pold[nthres+1] = gumbelmin_lccdf(thres[nthres] | mu);
    pnew[1] = gumbelmin_lcdf(thres[1]| 0);
    for (i in 2:nthres) {
      tmp = gumbelmin_lcdf(thres[i]| 0);
      pnew[i] = tmp + log1m_exp(gumbelmin_lcdf(thres[i-1]| 0) - tmp);
    }
    pnew[nthres+1] = gumbelmin_lccdf(thres[nthres] | 0);
    return multinomial_logit_lpmf(oldvec | pold) + 
      multinomial_logit_lpmf(newvec | pnew);
   }
"

gumbel8agglog_family <- custom_family(
  name = "gumbel8agglog", 
  dpars = c("mu", "crc", "crlm", "crll", "crlx", "crhm", "crhh", "crhx"), 
  links = c("identity", "identity", rep("log", 6)), 
  lb = c(NA, NA, rep(0, 6)),
  type = "int", vars = paste0("vint", 1:15, "[n]")
)
sv_gumbel8agglog <- stanvar(scode = gumbelmin_dist, block = "functions") + 
  stanvar(scode = gumbel8agglog_stanvars, block = "functions")

calc_posterior_predictions_gumbel8agglog <- function(i, prep) {
  gumbelmin <- function(x, mu, disc) {
    return(1 - extraDistr::pgumbel(-x,-mu,disc))
  }
  mu <- brms::get_dpar(prep, "mu", i = i)
  #discsignal <- brms::get_dpar(prep, "discsignal", i = i)
  crc <- brms::get_dpar(prep, "crc", i = i)
  crlm <- brms::get_dpar(prep, "crlm", i = i)
  crll <- brms::get_dpar(prep, "crll", i = i)
  crlx <- brms::get_dpar(prep, "crlx", i = i)
  crhm <- brms::get_dpar(prep, "crhm", i = i)
  crhh <- brms::get_dpar(prep, "crhh", i = i)
  crhx <- brms::get_dpar(prep, "crhx", i = i)

  oldvec <- prep$data$oldmat[i,]
  newvec <- prep$data$newmat[i,]
  OUTLEN <- length(mu)
  
  disc <- 1
  nthres <- 7
  thres <- matrix(NA_real_, nrow = OUTLEN, ncol = nthres)
  pold <- matrix(NA_real_, nrow = OUTLEN, ncol = nthres+1)
  pnew <- matrix(NA_real_, nrow = OUTLEN, ncol = nthres+1)
  
  thres[,1] = crc - (crlm + (crll) + (crlx));
  thres[,2] = crc - ((crlm) + (crll));
  thres[,3] = crc - ((crlm));
  thres[,4] = crc;
  thres[,5] = crc + ((crhm));
  thres[,6] = crc + ((crhm) + (crhh));
  thres[,7] = crc + ((crhm) + (crhh) + (crhx));
  
  # calculate probabilities
  pold[,1] = gumbelmin(thres[,1], mu, disc)
  #ordinal::pgumbel(thres[,1], mu, disc, max = FALSE) ## does NOT 
  for (j in 2:nthres) {
    pold[,j] = gumbelmin(thres[,j], mu, disc) - gumbelmin(thres[,j-1], mu, disc);
  }
  pold[,nthres+1] = 1 - gumbelmin(thres[,nthres], mu, disc);

  pnew[,1] = gumbelmin(thres[,1], 0, disc);
  for (j in 2:nthres) {
    pnew[,j] = gumbelmin(thres[,j], 0, disc) - gumbelmin(thres[,j-1], 0, disc);
  }
  pnew[,nthres+1] = 1 - gumbelmin(thres[,nthres], 0, disc);
  return(list(
    pold = pold, pnew = pnew
  ))
}

log_lik_gumbel8agglog <- function(i, prep) {
  use <- calc_posterior_predictions_gumbel8agglog(i = i, prep = prep)
  oldvec <- c(prep$data$Y[i], prep$data$vint1[i], prep$data$vint2[i], prep$data$vint3[i], 
              prep$data$vint4[i], prep$data$vint5[i], prep$data$vint6[i], prep$data$vint7[i])
  newvec <- c(prep$data$vint8[i], prep$data$vint9[i], prep$data$vint10[i], prep$data$vint11[i], 
              prep$data$vint12[i], prep$data$vint13[i], prep$data$vint14[i], prep$data$vint15[i])
  extraDistr::dmnom(x = oldvec, size = sum(oldvec), prob = use$pold, log = TRUE) + 
     extraDistr::dmnom(x = newvec, size = sum(newvec), prob = use$pnew, log = TRUE)
}

posterior_epred_gumbel8agglog <- function(prep) {
  nobs <- prep$nobs
  out <- array(NA_real_, dim = c(prep$ndraws, prep$nobs, 16), 
               dimnames = list(seq(prep$ndraws), seq(prep$nobs), 
                               c(colnames(prep$data$oldmat), colnames(prep$data$newmat))))
  for (i in seq_len(nobs)) {
    tmp <- calc_posterior_predictions_gumbel8agglog(i = i, prep = prep)
    out[,i,] <- c(tmp$pold, tmp$pnew)
  }
  return(out)
}

posterior_predict_gumbel8agglog <- function(i, prep, ...) {
  use <- calc_posterior_predictions_gumbel8agglog(i = i, prep = prep)
  oldvec <- c(prep$data$Y[i], prep$data$vint1[i], prep$data$vint2[i], prep$data$vint3[i], 
              prep$data$vint4[i], prep$data$vint5[i], prep$data$vint6[i], prep$data$vint7[i])
  newvec <- c(prep$data$vint8[i], prep$data$vint9[i], prep$data$vint10[i], prep$data$vint11[i], 
              prep$data$vint12[i], prep$data$vint13[i], prep$data$vint14[i], prep$data$vint15[i])
  
  lout <- nrow(use$pold)
  out <- cbind(extraDistr::rmnom(n = rep(1, lout), size = sum(oldvec), prob = use$pold), 
               extraDistr::rmnom(n = rep(1, lout), size = sum(oldvec), prob = use$pnew))
  colnames(out) <- c(colnames(prep$data$oldmat), colnames(prep$data$newmat))
  #browser()
  lapply(seq_len(nrow(out)), function(i) out[i,])
  #apply(out, 1, function(x) list(x))
  #out[,1]
}
