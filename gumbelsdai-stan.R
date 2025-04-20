gumbelsdai_stanvars <- '
  real pinternal1(real a, real b, real g1) {
    return exp(-exp(-g1 + a) + b);
  }
  real pinternal2(real a, real b, real g1, real g2) {
    return exp(-exp(-g1 + a) - exp(-g2 + a) + b);
  }
  real pfun(real l, real u, real g1, real g2) {
    real denominator;
    
    denominator = exp(g1) + exp(g2);
    
    return (pinternal1(l, g1, g1) - pinternal1(u, g1, g1) + pinternal1(l, g2, g1) - pinternal1(u, g2, g1) + pinternal2(u, g2, g1, g2) - pinternal2(l, g2, g1, g2)) / denominator;
  }
  real p_hit_gumbel(real l, real u) {
    return gumbelmin_cdf(u, 0)^2 - gumbelmin_cdf(l, 0)^2;
  }
  real gumbelsdai_lpmf(int y, real mu,
                   real crc, real crl, real crh, 
                   int y1, int y2, int y3, int y4, int y5, 
                   int y6, int y7, int y8, int y9, int y10, int y11) {
  int nthres = 3;
  vector[4] p_tabs;
  vector[8] p_tpres;
  array[4] int res_tabs = { y, y1, y2, y3 };
  array[8] int res_tpres = { y4, y5, y6, y7, y8, y9, y10, y11 };
  int problimit = 100;
  
  vector[nthres] thres;
  
  // calculate thresholds
  thres[1] = crc - crl;
  thres[2] = crc;
  thres[3] = crc + crh;
  
  p_tabs[1] = p_hit_gumbel(-problimit, thres[1]);
  p_tabs[2] = p_hit_gumbel(thres[1], thres[2]);
  p_tabs[3] = p_hit_gumbel(thres[2], thres[3]);
  p_tabs[4] = p_hit_gumbel(thres[3], problimit);
  
  p_tpres[1] = pfun(-problimit, thres[1], 0, mu);
  p_tpres[2] = pfun(thres[1], thres[2], 0, mu);
  p_tpres[3] = pfun(thres[2], thres[3],  0, mu);
  p_tpres[4] = pfun(thres[3], problimit, 0, mu);
  
  p_tpres[5] = pfun(-problimit, thres[1], mu, 0);
  p_tpres[6] = pfun(thres[1], thres[2], mu, 0);
  p_tpres[7] = pfun(thres[2], thres[3],  mu, 0);
  p_tpres[8] = pfun(thres[3], problimit, mu, 0);
  
  return multinomial_lpmf(res_tabs | p_tabs) + multinomial_lpmf(res_tpres | p_tpres);
  }
'
gumbelsdai_family <- custom_family(
  name = "gumbelsdai", 
  dpars = c("mu", "crc", "crl", "crh"), 
  links = c("identity", "identity", rep("log", 2)), 
  lb = c(NA, NA, rep(0, 2)), 
  type = "int", vars = paste0("vint", 1:11, "[n]")
)

source("gumbelmin_dist-stan.R")
sv_gumbelsdai <- stanvar(scode = gumbelmin_dist, block = "functions") +
  stanvar(scode = gumbelsdai_stanvars, block = "functions") 

calc_posterior_predictions_gumbelsdai <- function(i, prep) {
  pgumbmin = function(p, g=0){
    1 - exp(-exp(p-g))
  }
  dgumbmin = function(x, g=0){
    exp(x-g)*exp(-exp(x-g))
  }
  p_hit_gumb <- function(l, u, g1 = 0, g2=0){
    pgumbmin(u, g1)*pgumbmin(u, g2) - pgumbmin(l, g1)*pgumbmin(l, g2)
  }
  int_inst_gumb <- function(x, m, g1, g2){
    y = (pgumbmin(x, g2)^(m-1)) * dgumbmin(x, g1)
    return(ifelse(is.na(y) | is.nan(y) | is.infinite(y), 0, y))
  }
  m = 2
  
  mu <- brms::get_dpar(prep, "mu", i = i)
  crc <- brms::get_dpar(prep, "crc", i = i)
  crl <- brms::get_dpar(prep, "crl", i = i)
  crh <- brms::get_dpar(prep, "crh", i = i)
  
  OUTLEN <- length(mu)
  
  nthres <- 3
  thres <- matrix(NA_real_, nrow = OUTLEN, ncol = nthres)
  p_tabs <- matrix(NA_real_, nrow = OUTLEN, ncol = 4)
  p_tpres <- matrix(NA_real_, nrow = OUTLEN, ncol = 8)
  
  thres[,1] = crc - crl;
  thres[,2] = crc;
  thres[,3] = crc + crh;
  
  p_tabs[,1] = p_hit_gumb(-Inf, thres[,1]);
  p_tabs[,2] = p_hit_gumb(thres[,1], thres[,2]);
  p_tabs[,3] = p_hit_gumb(thres[,2], thres[,3]);
  p_tabs[,4] = p_hit_gumb(thres[,3], Inf);
  
  for (i in seq_len(OUTLEN)) {
    p_tpres[i, 1] = integrate(int_inst_gumb, -Inf, thres[i,1], 
                              m=m, g1=0, g2=mu[i], 
                              stop.on.error = F)$value
    p_tpres[i, 2] = integrate(int_inst_gumb, thres[i,1], thres[i,2], 
                              m=m, g1=0, g2=mu[i],
                              stop.on.error = F)$value
    p_tpres[i, 3] = integrate(int_inst_gumb, thres[i,2], thres[i,3], 
                              m=m, g1=0, g2=mu[i],
                              stop.on.error = F)$value
    p_tpres[i, 4] = integrate(int_inst_gumb, thres[i,3], Inf,
                              m=m, g1=0, g2=mu[i],
                              stop.on.error = F)$value
    p_tpres[i, 5] = integrate(int_inst_gumb, -Inf, thres[i,1], 
                              m=m, g2=0, g1=mu[i],
                              stop.on.error = F)$value
    p_tpres[i, 6] = integrate(int_inst_gumb, thres[i,1], thres[i,2], 
                              m=m, g2=0, g1=mu[i],
                              stop.on.error = F)$value
    p_tpres[i, 7] = integrate(int_inst_gumb, thres[i,2], thres[i,3], 
                              m=m, g2=0, g1=mu[i],
                              stop.on.error = F)$value
    p_tpres[i, 8] = integrate(int_inst_gumb, thres[i,3], Inf,
                              m=m, g2=0, g1=mu[i],
                              stop.on.error = F)$value
  }
  return(list(
    tabs = p_tabs, tpres = p_tpres
  ))
}


posterior_epred_gumbelsdai <- function(prep) {
  nobs <- prep$nobs
  out <- array(NA_real_, dim = c(prep$ndraws, prep$nobs, 12), 
               dimnames = list(seq(prep$ndraws), seq(prep$nobs), 
                               c(paste0("tabs", 1:4), paste0("tpres", 1:8))))
  for (i in seq_len(nobs)) {
    tmp <- calc_posterior_predictions_gumbelsdai(i = i, prep = prep)
    out[,i,] <- c(tmp$tabs, tmp$tpres)
  }
  return(out)
}


