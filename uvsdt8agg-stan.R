uvsdt8agg_stanvars <- "
  real uvsdt8agg_lpmf(int y, real mu, real discsignal, real crc, 
                    real crlm, real crll, real crlx, 
                    real crhm, real crhh, real crhx,
                    int y1, int y2, int y3, int y4, int y5, int y6, int y7, 
                    int y8, int y9, int y10, int y11, int y12, int y13, int y14, int y15) {
  int nthres = 7;
  real disc = 1/discsignal;
  vector[nthres+1] pold;
  vector[nthres+1] pnew;
  array[nthres+1] int oldvec = { y, y1, y2, y3, y4, y5, y6, y7 };
  array[nthres+1] int newvec = { y8, y9, y10, y11, y12, y13, y14, y15 };
  
  vector[nthres] thres;
  
  // calculate thresholds
  thres[1] = crc - (exp(crlm) + exp(crll) + exp(crlx));
  thres[2] = crc - (exp(crlm) + exp(crll));
  thres[3] = crc - (exp(crlm));
  thres[4] = crc;
  thres[5] = crc + (exp(crhm));
  thres[6] = crc + (exp(crhm) + exp(crhh));
  thres[7] = crc + (exp(crhm) + exp(crhh) + exp(crhx));
   
  // calculate probabilities
  pold[1] = Phi(disc * (thres[1] - mu));
  for (i in 2:nthres) {
    pold[i] = Phi(disc * (thres[i] - mu)) - Phi(disc * (thres[i - 1] - mu));
  }
  pold[nthres+1] = 1 - Phi(disc * (thres[nthres] - mu));
  
  pnew[1] = Phi((thres[1]));
  for (i in 2:nthres) {
    pnew[i] = Phi((thres[i])) - Phi((thres[i - 1]));
  }
  pnew[nthres+1] = 1 - Phi((thres[nthres]));
  return multinomial_lpmf(oldvec | pold) + multinomial_lpmf(newvec | pnew);
  }
"

uvsdt8agg_family <- custom_family(
  name = "uvsdt8agg", 
  dpars = c("mu", "discsignal", "crc", "crlm", "crll", "crlx", "crhm", "crhh", "crhx"), 
  links = c("identity", "log", rep("identity", 7)), lb = c(NA, 0, rep(NA, 7)),
  type = "int", vars = paste0("vint", 1:15, "[n]")
)
sv_uvsdt8agg <- stanvar(scode = uvsdt8agg_stanvars, block = "functions")

calc_posterior_predictions_uvsdt8agg <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  discsignal <- 1/brms::get_dpar(prep, "discsignal", i = i)
  crc <- brms::get_dpar(prep, "crc", i = i)
  crlm <- brms::get_dpar(prep, "crlm", i = i)
  crll <- brms::get_dpar(prep, "crll", i = i)
  crlx <- brms::get_dpar(prep, "crlx", i = i)
  crhm <- brms::get_dpar(prep, "crhm", i = i)
  crhh <- brms::get_dpar(prep, "crhh", i = i)
  crhx <- brms::get_dpar(prep, "crhx", i = i)

  OUTLEN <- length(mu)
  
  nthres <- 7
  thres <- matrix(NA_real_, nrow = OUTLEN, ncol = nthres)
  pold <- matrix(NA_real_, nrow = OUTLEN, ncol = nthres+1)
  pnew <- matrix(NA_real_, nrow = OUTLEN, ncol = nthres+1)
  
  thres[,1] = crc - (exp(crlm) + exp(crll) + exp(crlx));
  thres[,2] = crc - (exp(crlm) + exp(crll));
  thres[,3] = crc - (exp(crlm));
  thres[,4] = crc;
  thres[,5] = crc + (exp(crhm));
  thres[,6] = crc + (exp(crhm) + exp(crhh));
  thres[,7] = crc + (exp(crhm) + exp(crhh) + exp(crhx));
  
  # calculate probabilities
  pold[,1] = pnorm(discsignal * (thres[,1] - mu))
  for (j in 2:nthres) {
    pold[,j] = pnorm(discsignal * (thres[,j] - mu)) -
        pnorm(discsignal * (thres[,j - 1] - mu))
  }
  pold[,nthres+1] = 1 - pnorm(discsignal * (thres[,nthres] - mu));

  pnew[,1] = pnorm(thres[,1])
  for (j in 2:nthres) {
    pnew[,j] = pnorm(thres[,j]) - pnorm(thres[,j - 1])
  }
  pnew[,nthres+1] = 1 - pnorm(thres[,nthres]);
  return(list(
    pold = pold, pnew = pnew
  ))
}

log_lik_uvsdt8agg <- function(i, prep) {
  use <- calc_posterior_predictions_uvsdt8agg(i = i, prep = prep)
  oldvec <- c(prep$data$Y[i], prep$data$vint1[i], prep$data$vint2[i], prep$data$vint3[i], 
              prep$data$vint4[i], prep$data$vint5[i], prep$data$vint6[i], prep$data$vint7[i])
  newvec <- c(prep$data$vint8[i], prep$data$vint9[i], prep$data$vint10[i], prep$data$vint11[i], 
              prep$data$vint12[i], prep$data$vint13[i], prep$data$vint14[i], prep$data$vint15[i])
  extraDistr::dmnom(x = oldvec, size = sum(oldvec), prob = use$pold, log = TRUE) + 
     extraDistr::dmnom(x = newvec, size = sum(newvec), prob = use$pnew, log = TRUE)
}

posterior_epred_uvsdt8agg <- function(prep) {
  nobs <- prep$nobs
  out <- array(NA_real_, dim = c(prep$ndraws, prep$nobs, 16), 
               dimnames = list(seq(prep$ndraws), seq(prep$nobs), 
                               c(colnames(prep$data$oldmat), colnames(prep$data$newmat))))
  for (i in seq_len(nobs)) {
    tmp <- calc_posterior_predictions_uvsdt8agg(i = i, prep = prep)
    out[,i,] <- c(tmp$pold, tmp$pnew)
  }
  return(out)
}

posterior_predict_uvsdt8agg <- function(i, prep, ...) {
  use <- calc_posterior_predictions_uvsdt8agg(i = i, prep = prep)
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

