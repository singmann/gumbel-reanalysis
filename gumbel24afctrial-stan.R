gumbel24afctrial_stanvars <- "
   real gumbelmin(real x, real mu, real disc){
     //return 1- exp(-exp(-(-x-mu)/disc));
     //return exp(gumbel_lccdf(-x | mu,disc));
     return 1 - gumbel_cdf(-x|mu,disc);
   }
    real getp4(real x,             // Function argument
             real xc,            // Complement of function argument
                                //  on the domain (defined later)
             array[] real theta, // parameters
             array[] real x_r,   // data (real)
             array[] int x_i) {  // data (integer)
    real mu = theta[1];
  
    return ( (gumbelmin(x, 0, 1)^3) * exp(gumbel_lpdf(-x | mu, 1)) );
  }
  real getp2(real x,             // Function argument
             real xc,            // Complement of function argument
                                //  on the domain (defined later)
             array[] real theta, // parameters
             array[] real x_r,   // data (real)
             array[] int x_i) {  // data (integer)
    real mu = theta[1];
  
    return ( gumbelmin(x, 0, 1) * exp(gumbel_lpdf(-x | mu, 1)) );
  }

  real gumbel24afctrial_lpmf(int y, real mu,
                   int y1, 
                   data array[] real x_r, data array[] int x_i) {
  real p;
  
  if (y1 == 2) {
    p = integrate_1d(getp2, negative_infinity(),
                               positive_infinity(),
                               { -mu }, x_r, x_i);
  } else if (y1 == 4) {
    p = integrate_1d(getp4, negative_infinity(),
                             positive_infinity(),
                             { -mu }, x_r, x_i);
  }
  
    

  return bernoulli_lpmf(y | p);
  }
"

gumbel24afctrial_stanvars_tdata <- "
    array[0] real x_r;
    array[0] int x_i;
"

gumbel24afctrial_family <- custom_family(
  name = "gumbel24afctrial", 
  dpars = c("mu"), 
  links = c("identity"), lb = c(NA),
  type = "int", vars = c(paste0("vint", 1, "[n]"), "x_r", "x_i")
)
sv_gumbel24afctrial <- stanvar(scode = gumbel24afctrial_stanvars, block = "functions") +
  stanvar(scode = gumbel24afctrial_stanvars_tdata, block = "tdata")

calc_posterior_predictions_gumbel24afctrial <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  OUTLEN <- length(mu)
  afc <- prep$data$vint1[i]
  
  p <- vector("numeric", OUTLEN)
  
  G4<-function(x, mu){
    ((ordinal::pgumbel(x, max = FALSE)^3)*ordinal::dgumbel(x, mu, max = FALSE))
  }
  G2<-function(x, mu){
    ((ordinal::pgumbel(x, max = FALSE))*ordinal::dgumbel(x, mu, max = FALSE))
  }
  
  for (j in seq_len(OUTLEN)) {
    if (afc == 2) {
      p[j] <- integrate(G2,-Inf,Inf, mu = -mu[j],
                           rel.tol = .Machine$double.eps^0.5)$value  
    } else if (afc == 4) {
      p[j] <- integrate(G4,-Inf,Inf, mu = -mu[j],
                           rel.tol = .Machine$double.eps^0.5)$value
    }
  }
  return(p)
}

log_lik_gumbel24afctrial <- function(i, prep) {
  p <- calc_posterior_predictions_gumbel24afctrial(i = i, prep = prep)
  out <- p
  out[prep$data$Y[i] == 0] <- 1 - p
  log(out)
}

posterior_epred_gumbel24afctrial <- function(prep) {
  nobs <- prep$nobs
  out <- matrix(NA_real_, nrow = prep$ndraws, ncol = prep$nobs)
  for (i in seq_len(nobs)) {
    tmp <- calc_posterior_predictions_gumbel24afctrial(i = i, prep = prep)
    out[,i] <- tmp
  }
  return(out)
}

posterior_predict_gumbel24afctrial <- function(i, prep, ...) {
  p <- calc_posterior_predictions_gumbel24afctrial(i = i, prep = prep)
  rbinom(length(p), 1, p)
}

