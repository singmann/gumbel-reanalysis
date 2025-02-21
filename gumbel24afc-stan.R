gumbel24afc_stanvars <- "
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

  real gumbel24afc_lpmf(int y, real mu,
                   int y1, int y2, int y3,
                   data array[] real x_r, data array[] int x_i) {
  real p2;
  real p4;
  
  p2 = integrate_1d(getp2, negative_infinity(),
                             positive_infinity(),
                             { -mu }, x_r, x_i);
  p4 = integrate_1d(getp4, negative_infinity(),
                             positive_infinity(),
                             { -mu }, x_r, x_i);

  return binomial_lpmf(y | y1, p2) + binomial_lpmf(y2 | y3, p4);
  }
"

gumbel24afc_stanvars_tdata <- "
    array[0] real x_r;
    array[0] int x_i;
"

gumbel24afc_family <- custom_family(
  name = "gumbel24afc", 
  dpars = c("mu"), 
  links = c("identity"), lb = c(NA),
  type = "int", vars = c(paste0("vint", 1:3, "[n]"), "x_r", "x_i")
)
sv_gumbel24afc <- stanvar(scode = gumbel24afc_stanvars, block = "functions") +
  stanvar(scode = gumbel24afc_stanvars_tdata, block = "tdata")

calc_posterior_predictions_gumbel24afc <- function(i, prep) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  OUTLEN <- length(mu)
  
  p <- matrix(NA_real_, nrow = OUTLEN, ncol = 2)
  
  G4<-function(x, mu){
    ((ordinal::pgumbel(x, max = FALSE)^3)*ordinal::dgumbel(x, mu, max = FALSE))
  }
  G2<-function(x, mu){
    ((ordinal::pgumbel(x, max = FALSE))*ordinal::dgumbel(x, mu, max = FALSE))
  }
  
  for (j in seq_len(OUTLEN)) {
    p[j, 1] <- integrate(G2,-Inf,Inf, mu = -mu[j],
                         rel.tol = .Machine$double.eps^0.5)$value
    p[j, 2] <- integrate(G4,-Inf,Inf, mu = -mu[j],
                         rel.tol = .Machine$double.eps^0.5)$value
  }
  return(p)
}

log_lik_gumbel24afc <- function(i, prep) {
  p <- calc_posterior_predictions_gumbel24afc(i = i, prep = prep)
  dvec <- c(prep$data$Y[i], prep$data$vint1[i], prep$data$vint2[i], 
              prep$data$vint3[i])
  dbinom(prep$data$Y[i], prep$data$vint1[i], p[,1], log = TRUE) + 
    dbinom(prep$data$vint2[i], prep$data$vint3[i], p[,2], log = TRUE)
}

posterior_epred_gumbel24afc <- function(prep) {
  nobs <- prep$nobs
  out <- array(NA_real_, dim = c(prep$ndraws, prep$nobs, 2), 
               dimnames = list(seq(prep$ndraws), seq(prep$nobs), 
                               c("AFC2", "AFC4")))
  for (i in seq_len(nobs)) {
    tmp <- calc_posterior_predictions_gumbel24afc(i = i, prep = prep)
    out[,i,] <- tmp
  }
  return(out)
}

## not yet implemented correctly
# posterior_predict_gumbel24afc <- function(i, prep, ...) {
#   p <- calc_posterior_predictions_gumbel24afc(i = i, prep = prep)
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

