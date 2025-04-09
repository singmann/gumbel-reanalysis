uvsdtrank_stanvars <- "
  real getp(real x,             // Function argument
             real xc,            // Complement of function argument
                                //  on the domain (defined later)
             array[] real theta, // parameters
             array[] real x_r,   // data (real)
             array[] int x_i) {  // data (integer)
    real mu = theta[1];
    real sigma = theta[2];
    int r = x_i[1];
    int m = x_i[2];
  
    return( Phi(x)^(m-r) * exp(normal_lpdf(x | mu, sigma)) *  (1 - Phi(x))^(r-1)  ) ;
  }
  
  real uvsdtrank_lpmf(int y, real mu, real discsignal,
                      int r, int maxrank, 
                      data array[] real x_r, data array[] int x_i) {
  real p;
  real g = mu;
  int m = maxrank;
  real disc = 1/discsignal;
  
  if (y == 0) {
    return 0;
  } else {
    p = choose(m-1, r-1) * integrate_1d(getp, negative_infinity(),
                             positive_infinity(),
                             { mu, disc }, x_r, {r, maxrank});
    return y * log(p);
  }
  
  }
"

uvsdtrank_stanvars_tdata <- "
    array[0] real x_r;
    array[0] int x_i;
"

uvsdtrank_family <- custom_family(
  name = "uvsdtrank", 
  dpars = c("mu", "discsignal"), 
  links = c("identity", "log"), lb = c(NA, 0),
  type = "int", vars = c(paste0("vint", 1:2, "[n]"), "x_r", "x_i")
)
sv_uvsdtrank <- stanvar(scode = uvsdtrank_stanvars, block = "functions") +
  stanvar(scode = uvsdtrank_stanvars_tdata, block = "tdata")

calc_posterior_predictions_uvsdtrank <- function(i, prep) {
  int_Rr_uvg = function(x, d, s, m, r){
    pnorm(x, lower.tail=F)^(r-1) * dnorm(x, d, s) * pnorm(x)^(m-r)
  }
  
  mu <- brms::get_dpar(prep, "mu", i = i)
  discsignal <- brms::get_dpar(prep, "discsignal", i = i)

  OUTLEN <- length(mu)
  
  r <- prep$data$vint1[i]
  m <- prep$data$vint2[i]
  multiplier <- choose(m-1, r-1)
  
  p <- vector(mode = "numeric", OUTLEN)
  
  for (j in seq_len(OUTLEN)) {
    p[j] <- multiplier * integrate(int_Rr_uvg, -Inf, Inf, 
                                   d=mu[j], s=1/discsignal[j], 
                                   m=m, r=r)$value
  }
  return(p)
}

### need to make sure log-likelihood has correct normalisation constant
### this will be inefficient as it needs to be calculated multiple times
### for each multinomial distribution
# log_lik_uvsdtrank <- function(i, prep) {
#   p <- calc_posterior_predictions_uvsdtrank(i = i, prep = prep)
#   dvec <- c(prep$data$Y[i], prep$data$vint1[i], prep$data$vint2[i], 
#               prep$data$vint3[i])
#   extraDistr::dmnom(x = dvec, size = sum(dvec), prob = p, log = TRUE)
# }

posterior_epred_uvsdtrank <- function(prep) {
  nobs <- prep$nobs
  out <- matrix(NA_real_, prep$ndraws, prep$nobs,
                dimnames = list(seq(prep$ndraws), seq(prep$nobs)))
  for (i in seq_len(nobs)) {
    tmp <- calc_posterior_predictions_uvsdtrank(i = i, prep = prep)
    out[,i] <- tmp
  }
  return(out)
}

### same as for log-lik, functions need adjustment. Unclear if it can be fixed
# posterior_predict_uvsdtrank <- function(i, prep, ...) {
#   p <- calc_posterior_predictions_uvsdtrank(i = i, prep = prep)
#   dvec <- c(prep$data$Y[i], prep$data$vint1[i], prep$data$vint2[i], 
#               prep$data$vint3[i])
#   
#   lout <- length(p)
#   out <- extraDistr::rmnom(n = rep(1, lout), size = sum(dvec), prob = p)  
#   colnames(out) <- c("R1", "R2", "R3", "R4")
#   #browser()
#   lapply(seq_len(nrow(out)), function(i) out[i,])
#   #apply(out, 1, function(x) list(x))
#   #out[,1]
# }

