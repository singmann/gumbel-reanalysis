uvsdtsdai_stanvars <- '
   //real int_inst_uvg(real x, real d1, real s1, real d2, real s2) {
    real int_inst_uvg(real x,             // Function argument
             real xc,            // Complement of function argument
                                //  on the domain (defined later)
             array[] real theta, // parameters
             array[] real x_r,   // data (real)
             array[] int x_i) {  // data (integer)
     real d1 = theta[1];
     real s1 = theta[2];
     real d2 = theta[3];
     real s2 = theta[4];
     real m = 2;
     return (normal_cdf(x | d2, s2)^(m-1)) * exp(normal_lpdf(x | d1, s1));
   }
  real p_hit_uvg(real l, real u) {
    // pnorm(u, d1, s1)*pnorm(u, d2, s2) - pnorm(l, d1, s1)*pnorm(l, d2, s2)
    return std_normal_cdf(u)^2 - std_normal_cdf(l)^2;
  }
  real uvsdtsdai_lpmf(int y, real mu, real discsignal, 
                   real crc, real crl, real crh, 
                   int y1, int y2, int y3, int y4, int y5, 
                   int y6, int y7, int y8, int y9, int y10, int y11,
                   data array[] real x_r, data array[] int x_i) {
  int nthres = 3;
  real disc = 1/discsignal;
  vector[4] p_tabs;
  vector[8] p_tpres;
  array[4] int res_tabs = { y, y1, y2, y3 };
  array[8] int res_tpres = { y4, y5, y6, y7, y8, y9, y10, y11 };
  
  vector[nthres] thres;
  
  // calculate thresholds
  thres[1] = crc - (exp(crl));
  thres[2] = crc;
  thres[3] = crc + (exp(crh));
  
  p_tabs[1] = p_hit_uvg(negative_infinity(), thres[1]);
  p_tabs[2] = p_hit_uvg(thres[1], thres[2]);
  p_tabs[3] = p_hit_uvg(thres[2], thres[3]);
  p_tabs[4] = p_hit_uvg(thres[3], positive_infinity());
  
  p_tpres[1] = integrate_1d(int_inst_uvg, negative_infinity(),
                         thres[1],
                         { 0, 1, mu, disc }, x_r, x_i);
  p_tpres[2] = integrate_1d(int_inst_uvg, thres[1],
                         thres[2],
                         { 0, 1, mu, disc }, x_r, x_i);
  p_tpres[3] = integrate_1d(int_inst_uvg, thres[2],
                         thres[3],
                         { 0, 1, mu, disc }, x_r, x_i);
  p_tpres[4] = integrate_1d(int_inst_uvg, thres[3],
                         positive_infinity(),
                         { 0, 1, mu, disc }, x_r, x_i);
  
  p_tpres[5] = integrate_1d(int_inst_uvg, negative_infinity(),
                         thres[1],
                         { mu, disc, 0, 1 }, x_r, x_i);
  p_tpres[6] = integrate_1d(int_inst_uvg, thres[1],
                         thres[2],
                         { mu, disc, 0, 1 }, x_r, x_i);
  p_tpres[7] = integrate_1d(int_inst_uvg, thres[2],
                         thres[3],
                         { mu, disc, 0, 1 }, x_r, x_i);
  p_tpres[8] = integrate_1d(int_inst_uvg, thres[3],
                         positive_infinity(),
                         { mu, disc, 0, 1 }, x_r, x_i);
  
  return multinomial_lpmf(res_tabs | p_tabs) + multinomial_lpmf(res_tpres | p_tpres);
  }
'
uvsdtsdai_family <- custom_family(
  name = "uvsdtsdai", 
  dpars = c("mu", "discsignal", "crc", "crl", "crh"), 
  links = c("identity", "log", rep("identity", 3)), lb = c(NA, 0, rep(NA, 3)),
  type = "int", vars = c(paste0("vint", 1:11, "[n]"), "x_r", "x_i")
)
uvsdtsdai_stanvars_tdata <- "
    array[0] real x_r;
    array[0] int x_i;
"

sv_uvsdtsdai <- stanvar(scode = uvsdtsdai_stanvars, block = "functions") +
  stanvar(scode = uvsdtsdai_stanvars_tdata, block = "tdata")

calc_posterior_predictions_uvsdtsdai <- function(i, prep) {
  p_hit_uvg <- function(l, u, d1 = 0, s1 = 1, d2=0, s2=1){
    pnorm(u, d1, s1)*pnorm(u, d2, s2) - pnorm(l, d1, s1)*pnorm(l, d2, s2)
  }
  int_inst_uvg <- function(x, m, d1, s1, d2, s2){
    y = (pnorm(x, d2, s2)^(m-1)) * dnorm(x, d1, s1)
    return(y)
  }
  m = 2
  
  mu <- brms::get_dpar(prep, "mu", i = i)
  discsignal <- 1/brms::get_dpar(prep, "discsignal", i = i)
  crc <- brms::get_dpar(prep, "crc", i = i)
  crl <- brms::get_dpar(prep, "crl", i = i)
  crh <- brms::get_dpar(prep, "crh", i = i)
  
  OUTLEN <- length(mu)
  
  nthres <- 3
  thres <- matrix(NA_real_, nrow = OUTLEN, ncol = nthres)
  p_tabs <- matrix(NA_real_, nrow = OUTLEN, ncol = 4)
  p_tpres <- matrix(NA_real_, nrow = OUTLEN, ncol = 8)
  
  thres[,1] = crc - (exp(crl));
  thres[,2] = crc;
  thres[,3] = crc + (exp(crh));
  
  p_tabs[,1] = p_hit_uvg(-Inf, thres[,1]);
  p_tabs[,2] = p_hit_uvg(thres[,1], thres[,2]);
  p_tabs[,3] = p_hit_uvg(thres[,2], thres[,3]);
  p_tabs[,4] = p_hit_uvg(thres[,3], Inf);
  
  for (i in seq_len(OUTLEN)) {
    p_tpres[i, 1] = integrate(int_inst_uvg, -Inf, thres[i,1], 
                              m=m, d1=0, s1=1, d2=mu[i], s2=discsignal[i],
                              stop.on.error = F)$value
    p_tpres[i, 2] = integrate(int_inst_uvg, thres[i,1], thres[i,2], 
                              m=m, d1=0, s1=1, d2=mu[i], s2=discsignal[i],
                              stop.on.error = F)$value
    p_tpres[i, 3] = integrate(int_inst_uvg, thres[i,2], thres[i,3], 
                              m=m, d1=0, s1=1, d2=mu[i], s2=discsignal[i],
                              stop.on.error = F)$value
    p_tpres[i, 4] = integrate(int_inst_uvg, thres[i,3], Inf,
                              m=m, d1=0, s1=1, d2=mu[i], s2=discsignal[i],
                              stop.on.error = F)$value
    p_tpres[i, 5] = integrate(int_inst_uvg, -Inf, thres[i,1], 
                              m=m, d2=0, s2=1, d1=mu[i], s1=discsignal[i],
                              stop.on.error = F)$value
    p_tpres[i, 6] = integrate(int_inst_uvg, thres[i,1], thres[i,2], 
                              m=m, d2=0, s2=1, d1=mu[i], s1=discsignal[i],
                              stop.on.error = F)$value
    p_tpres[i, 7] = integrate(int_inst_uvg, thres[i,2], thres[i,3], 
                              m=m, d2=0, s2=1, d1=mu[i], s1=discsignal[i],
                              stop.on.error = F)$value
    p_tpres[i, 8] = integrate(int_inst_uvg, thres[i,3], Inf,
                              m=m, d2=0, s2=1, d1=mu[i], s1=discsignal[i],
                              stop.on.error = F)$value
  }
  return(list(
    tabs = p_tabs, tpres = p_tpres
  ))
}

# log_lik_uvsdtsdai <- function(i, prep) {
#   use <- calc_posterior_predictions_uvsdtsdai(i = i, prep = prep)
#   oldvec <- c(prep$data$Y[i], prep$data$vint1[i], prep$data$vint2[i], 
#               prep$data$vint3[i], prep$data$vint4[i], prep$data$vint5[i])
#   newvec <- c(prep$data$vint6[i], prep$data$vint7[i], prep$data$vint8[i], 
#               prep$data$vint9[i], prep$data$vint10[i], prep$data$vint11[i])
#   extraDistr::dmnom(x = oldvec, size = sum(oldvec), prob = use$pold, log = TRUE) + 
#     extraDistr::dmnom(x = newvec, size = sum(newvec), prob = use$pnew, log = TRUE)
# }


posterior_epred_uvsdtsdai <- function(prep) {
  nobs <- prep$nobs
  out <- array(NA_real_, dim = c(prep$ndraws, prep$nobs, 12), 
               dimnames = list(seq(prep$ndraws), seq(prep$nobs), 
                               c(paste0("tabs", 1:4), paste0("tpres", 1:8))))
  for (i in seq_len(nobs)) {
    tmp <- calc_posterior_predictions_uvsdtsdai(i = i, prep = prep)
    out[,i,] <- c(tmp$tabs, tmp$tpres)
  }
  return(out)
}

# posterior_predict_uvsdtsdai <- function(i, prep, ...) {
#   use <- calc_posterior_predictions_uvsdtsdai(i = i, prep = prep)
#   oldvec <- c(prep$data$Y[i], prep$data$vint1[i], prep$data$vint2[i], 
#               prep$data$vint3[i], prep$data$vint4[i], prep$data$vint5[i])
#   newvec <- c(prep$data$vint6[i], prep$data$vint7[i], prep$data$vint8[i], 
#               prep$data$vint9[i], prep$data$vint10[i], prep$data$vint11[i])
#   
#   lout <- nrow(use$pold)
#   out <- cbind(extraDistr::rmnom(n = rep(1, lout), size = sum(oldvec), prob = use$pold), 
#                extraDistr::rmnom(n = rep(1, lout), size = sum(oldvec), prob = use$pnew))
#   colnames(out) <- c(colnames(prep$data$oldmat), colnames(prep$data$newmat))
#   #browser()
#   lapply(seq_len(nrow(out)), function(i) out[i,])
#   #apply(out, 1, function(x) list(x))
#   #out[,1]
# }

