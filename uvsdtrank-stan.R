uvsdtrank_stanvars <- "
  real getp1(real x,             // Function argument
             real xc,            // Complement of function argument
                                //  on the domain (defined later)
             array[] real theta, // parameters
             array[] real x_r,   // data (real)
             array[] int x_i) {  // data (integer)
    real mu = theta[1];
    real sigma = theta[2];
  
    return ( (Phi(x)^3) * exp(normal_lpdf(x | mu, sigma)) );
  }
  real getp2(real x,             // Function argument
             real xc,            // Complement of function argument
                                //  on the domain (defined later)
             array[] real theta, // parameters
             array[] real x_r,   // data (real)
             array[] int x_i) {  // data (integer)
    real mu = theta[1];
    real sigma = theta[2];
  
    return ( (Phi(x)^2) * exp(normal_lpdf(x | mu, sigma)) * (1 - Phi(x)) ) * 3;
  }
  real getp3(real x,             // Function argument
             real xc,            // Complement of function argument
                                //  on the domain (defined later)
             array[] real theta, // parameters
             array[] real x_r,   // data (real)
             array[] int x_i) {  // data (integer)
    real mu = theta[1];
    real sigma = theta[2];
  
    return ( (Phi(x)) * exp(normal_lpdf(x | mu, sigma)) * (1 - Phi(x))^2 ) * 3;
  }
  
  real uvsdtrank_lpmf(int y, real mu, real discsignal, 
                   int y1, int y2, int y3, 
                   data array[] real x_r, data array[] int x_i) {
  real disc = 1/discsignal;
  vector[4] p;
  array[4] int respvec = { y, y1, y2, y3 };
  
  // ((pnorm(x)^3)*dnorm(x,mean=mu,sd=ss))
  p[1] = integrate_1d(getp1, negative_infinity(),
                             positive_infinity(),
                             { mu, disc }, x_r, x_i);
  p[2] = integrate_1d(getp2, negative_infinity(),
                             positive_infinity(),
                             { mu, disc }, x_r, x_i);
  p[3] = integrate_1d(getp3, negative_infinity(),
                             positive_infinity(),
                             { mu, disc }, x_r, x_i);
  p[4] = 1-p[1]-p[2]-p[3];

  return multinomial_lpmf(respvec | p);
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
  type = "int", vars = c(paste0("vint", 1:3, "[n]"), "x_r", "x_i")
)
sv_uvsdtrank <- stanvar(scode = uvsdtrank_stanvars, block = "functions") +
  stanvar(scode = uvsdtrank_stanvars_tdata, block = "tdata")

calc_posterior_predictions_uvsdtrank <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  discsignal <- 1/brms::get_dpar(prep, "discsignal", i = i)

  OUTLEN <- length(mu)
  
  p <- matrix(NA_real_, nrow = OUTLEN, ncol = 4)
  
  G1<-function(x, mu, sd){
    ((pnorm(x)^3)*dnorm(x,mean=mu,sd=ss))
  }
  
  G2<-function(x, mu, sd){
    ((pnorm(x)^2)*dnorm(x,mean=mu,sd=ss)*(1-pnorm(x)))*3
  }
  
  G3<-function(x, mu, sd){
    (pnorm(x)*dnorm(x,mean=mu,sd=ss)*(1-pnorm(x))^2)*3
  }

  for (j in seq_len(OUTLEN)) {
    p[j, 1] <- integrate(G1,-Inf,Inf, mu = mu[j], sd = discsignal[j],
                         rel.tol = .Machine$double.eps^0.5)$value    
    p[j, 2] <- integrate(G2,-Inf,Inf, mu = mu[j], sd = discsignal[j],
                         rel.tol = .Machine$double.eps^0.5)$value
    p[j, 3] <- integrate(G3,-Inf,Inf, mu = mu[j], sd = discsignal[j],
                         rel.tol = .Machine$double.eps^0.5)$value
  }
  p[, 4] <- 1-p[,1]-p[,2]-p[,3]
  return(p)
}

log_lik_uvsdtrank <- function(i, prep) {
  p <- calc_posterior_predictions_uvsdtrank(i = i, prep = prep)
  dvec <- c(prep$data$Y[i], prep$data$vint1[i], prep$data$vint2[i], 
              prep$data$vint3[i])
  extraDistr::dmnom(x = dvec, size = sum(dvec), prob = p, log = TRUE)
}

posterior_epred_uvsdtrank <- function(prep) {
  nobs <- prep$nobs
  out <- array(NA_real_, dim = c(prep$ndraws, prep$nobs, 4), 
               dimnames = list(seq(prep$ndraws), seq(prep$nobs), 
                               c("R1", "R2", "R3", "R4")))
  for (i in seq_len(nobs)) {
    tmp <- calc_posterior_predictions_uvsdtrank(i = i, prep = prep)
    out[,i,] <- tmp
  }
  return(out)
}

posterior_predict_uvsdtrank <- function(i, prep, ...) {
  p <- calc_posterior_predictions_uvsdtrank(i = i, prep = prep)
  dvec <- c(prep$data$Y[i], prep$data$vint1[i], prep$data$vint2[i], 
              prep$data$vint3[i])
  
  lout <- length(p)
  out <- extraDistr::rmnom(n = rep(1, lout), size = sum(dvec), prob = p)  
  colnames(out) <- c("R1", "R2", "R3", "R4")
  #browser()
  lapply(seq_len(nrow(out)), function(i) out[i,])
  #apply(out, 1, function(x) list(x))
  #out[,1]
}

