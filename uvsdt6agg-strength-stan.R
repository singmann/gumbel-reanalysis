uvsdt6aggreg_stanvars <- "
  real uvsdt6aggreg_lpmf(int y, real mu, real discsignal, 
                   real crc, real crlm, real crll, real crhm, real crhh, 
                   int y1, int y2, int y3, int y4, int y5, 
                   int type) {
  int nthres = 5;
  real disc = 1/discsignal;
  vector[nthres+1] p;
  array[6] int resvec = { y, y1, y2, y3, y4, y5 };
  
  vector[nthres] thres;
  
    // calculate thresholds
    thres[1] = crc - (crlm + crll);
    thres[2] = crc - (crlm);
    thres[3] = crc;
    thres[4] = crc + (crhm);
    thres[5] = crc + (crhm + crhh);
   
  // calculate probabilities
  if (type == 1) {
    p[1] = Phi(disc * (thres[1] - mu));
    for (i in 2:nthres) {
      p[i] = Phi(disc * (thres[i] - mu)) - Phi(disc * (thres[i - 1] - mu));
    }
    p[6] = 1 - Phi(disc * (thres[nthres] - mu));
  } else if (type == 0) {
    p[1] = Phi((thres[1]));
    for (i in 2:nthres) {
      p[i] = Phi((thres[i])) - Phi((thres[i - 1]));
    }
    p[6] = 1 - Phi((thres[nthres]));
  }
  
  return multinomial_lpmf(resvec | p);
  }
"

uvsdt6aggreg_family <- custom_family(
  name = "uvsdt6aggreg", 
  dpars = c("mu", "discsignal", "crc", "crlm", "crll", "crhm", "crhh"), 
  links = c("identity", "log", "identity", rep("log", 4)), 
  lb = c(NA, 0, NA, rep(0, 4)),
  type = "int", vars = paste0("vint", 1:6, "[n]")
)
sv_uvsdt6aggreg <- stanvar(scode = uvsdt6aggreg_stanvars, block = "functions")

calc_posterior_predictions_uvsdt6aggreg <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  discsignal <- 1/brms::get_dpar(prep, "discsignal", i = i)
  crc <- brms::get_dpar(prep, "crc", i = i)
  crlm <- brms::get_dpar(prep, "crlm", i = i)
  crll <- brms::get_dpar(prep, "crll", i = i)
  crhm <- brms::get_dpar(prep, "crhm", i = i)
  crhh <- brms::get_dpar(prep, "crhh", i = i)
  type <- prep$data$vint6[i]

  OUTLEN <- length(mu)
  
  nthres <- 5
  thres <- matrix(NA_real_, nrow = OUTLEN, ncol = nthres)
  pold <- matrix(NA_real_, nrow = OUTLEN, ncol = nthres+1)
  pnew <- matrix(NA_real_, nrow = OUTLEN, ncol = nthres+1)
  p <- matrix(NA_real_, nrow = OUTLEN, ncol = nthres+1)
  
  thres[,1] = crc - ((crlm) + (crll));
  thres[,2] = crc - ((crlm));
  thres[,3] = crc;
  thres[,4] = crc + ((crhm));
  thres[,5] = crc + ((crhm) + (crhh));
  
  # calculate probabilities
  if (type == 1) {
    p[,1] = pnorm(discsignal * (thres[,1] - mu))
    for (j in 2:nthres) {
      p[,j] = pnorm(discsignal * (thres[,j] - mu)) -
          pnorm(discsignal * (thres[,j - 1] - mu))
    }
    p[,6] = 1 - pnorm(discsignal * (thres[,nthres] - mu));
  } else if (type == 0) {
    p[,1] = pnorm(thres[,1])
    for (j in 2:nthres) {
      p[,j] = pnorm(thres[,j]) - pnorm(thres[,j - 1])
    }
    p[,6] = 1 - pnorm(thres[,nthres]);
  }
  return(p)
}

log_lik_uvsdt6aggreg <- function(i, prep) {
  use <- calc_posterior_predictions_uvsdt6aggreg(i = i, prep = prep)
  resvec <- c(prep$data$Y[i], prep$data$vint1[i], prep$data$vint2[i], 
              prep$data$vint3[i], prep$data$vint4[i], prep$data$vint5[i])
  extraDistr::dmnom(x = resvec, size = sum(resvec), prob = use, log = TRUE) 
}

posterior_epred_uvsdt6aggreg <- function(prep) {
  nobs <- prep$nobs
  out <- array(NA_real_, dim = c(prep$ndraws, prep$nobs, 6), 
               dimnames = list(seq(prep$ndraws), seq(prep$nobs), 
                               paste0("r", 1:6)))
  for (i in seq_len(nobs)) {
    tmp <- calc_posterior_predictions_uvsdt6aggreg(i = i, prep = prep)
    out[,i,] <- tmp
  }
  return(out)
}

