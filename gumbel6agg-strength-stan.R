source("gumbelmin_dist-stan.R")
## cat(gumbelmin_dist)
gumbel6aggreg_stanvars <- "
   real gumbel6aggreg_lpmf(int y, real mu, 
                   real crc, real crlm, real crll, real crhm, real crhh, 
                   int y1, int y2, int y3, int y4, int y5, 
                   int type) {
    int nthres = 5;
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
      p[1] = gumbelmin_cdf(thres[1], mu);
      for (i in 2:nthres) {
        p[i] = gumbelmin_cdf(thres[i], mu) - gumbelmin_cdf(thres[i-1], mu);
      }
      p[6] = 1 - gumbelmin_cdf(thres[nthres], mu);
    } else if (type == 0) {
      p[1] = gumbelmin_cdf(thres[1], 0);
      for (i in 2:nthres) {
        p[i] = gumbelmin_cdf(thres[i], 0) - gumbelmin_cdf(thres[i-1], 0);
      }
      p[6] = 1 - gumbelmin_cdf(thres[nthres], 0);
    }
    
    return multinomial_lpmf(resvec | p);
   }
"

gumbel6aggreg_family <- custom_family(
  name = "gumbel6aggreg", 
  dpars = c("mu", "crc", "crlm", "crll", "crhm", "crhh"), 
  links = c("identity", "identity", rep("log", 4)), 
  lb = c(NA, NA, rep(0, 4)),
  type = "int", vars = paste0("vint", 1:6, "[n]")
)

sv_gumbel6aggreg <- stanvar(scode = gumbelmin_dist, block = "functions") + 
  stanvar(scode = gumbel6aggreg_stanvars, block = "functions")

calc_posterior_predictions_gumbel6aggreg <- function(i, prep) {
  gumbelmin <- function(x, mu, disc) {
    return(1 - extraDistr::pgumbel(-x,-mu,disc))
  }
  mu <- brms::get_dpar(prep, "mu", i = i)
  #discsignal <- brms::get_dpar(prep, "discsignal", i = i)
  crc <- brms::get_dpar(prep, "crc", i = i)
  crlm <- brms::get_dpar(prep, "crlm", i = i)
  crll <- brms::get_dpar(prep, "crll", i = i)
  crhm <- brms::get_dpar(prep, "crhm", i = i)
  crhh <- brms::get_dpar(prep, "crhh", i = i)
  type <- prep$data$vint6[i]

  OUTLEN <- length(mu)
  
  disc <- 1
  nthres <- 5
  thres <- matrix(NA_real_, nrow = OUTLEN, ncol = nthres)
  p <- matrix(NA_real_, nrow = OUTLEN, ncol = nthres+1)
  
  thres[,1] = crc - ((crlm) + (crll));
  thres[,2] = crc - ((crlm));
  thres[,3] = crc;
  thres[,4] = crc + ((crhm));
  thres[,5] = crc + ((crhm) + (crhh));
  
  # calculate probabilities
  if (type == 1) {
    p[,1] = gumbelmin(thres[,1], mu, disc)
    #ordinal::pgumbel(thres[,1], mu, disc, max = FALSE) ## does NOT 
    for (j in 2:nthres) {
      p[,j] = gumbelmin(thres[,j], mu, disc) - gumbelmin(thres[,j-1], mu, disc);
    }
    p[,6] = 1 - gumbelmin(thres[,nthres], mu, disc);
  } else if (type == 0) {
    p[,1] = gumbelmin(thres[,1], 0, disc);
    for (j in 2:nthres) {
      p[,j] = gumbelmin(thres[,j], 0, disc) - gumbelmin(thres[,j-1], 0, disc);
    }
    p[,6] = 1 - gumbelmin(thres[,nthres], 0, disc);
  }



  return(p)
}

log_lik_gumbel6aggreg <- function(i, prep) {
  use <- calc_posterior_predictions_gumbel6aggreg(i = i, prep = prep)
  resvec <- c(prep$data$Y[i], prep$data$vint1[i], prep$data$vint2[i], 
              prep$data$vint3[i], prep$data$vint4[i], prep$data$vint5[i])
  extraDistr::dmnom(x = resvec, size = sum(resvec), prob = use, log = TRUE) 
}

posterior_epred_gumbel6aggreg <- function(prep) {
  nobs <- prep$nobs
  out <- array(NA_real_, dim = c(prep$ndraws, prep$nobs, 6), 
               dimnames = list(seq(prep$ndraws), seq(prep$nobs), 
                               paste0("r", 1:6)))
  for (i in seq_len(nobs)) {
    tmp <- calc_posterior_predictions_gumbel6aggreg(i = i, prep = prep)
    out[,i,] <- tmp
  }
  return(out)
}
